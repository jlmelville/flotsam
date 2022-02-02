# flotsam

Finicky LOcal Tangent Space Alignment Method

<!-- badges: start -->
[![R-CMD-check](https://github.com/jlmelville/flotsam/workflows/R-CMD-check/badge.svg)](https://github.com/jlmelville/flotsam/actions)
[![Codecov test coverage](https://codecov.io/gh/jlmelville/flotsam/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jlmelville/flotsam?branch=main)
<!-- badges: end -->

An implementation of the [Local Tangent Space
Alignment](https://doi.org/10.1137/S1064827502419154) method (Zhang & Zha,
2004), a spectral method for dimensionality reduction. It is one of the few
methods that can unroll the "swiss roll" data set and its variants without
major distortion.

The F is for 'Finicky' (see 'Current Status'), but one day I would like it to
be Fast.

## Installation

``` r
devtools::install_github("jlmelville/flotsam")
```

This package relies on a package
([rnndescent](https://www.github.com/jlmelville/rnndescent)) which is not on
CRAN.

## Example

``` r
library(flotsam)

# Create a "swiss roll": 2D rectangle rolled up in 3D
n <- 1000
max_z <- 10

phi <- stats::runif(n, min = 1.5 * pi, max = 4.5 * pi)
x <- phi * cos(phi)
y <- phi * sin(phi)
z <- stats::runif(n, max = max_z)
swiss_roll <- data.frame(x, y, z)

# See side section to prove it's rolled
plot(swiss_roll$x, swiss_roll$y, col = phi)

# unroll it
swiss_unrolled <-
  ltsa(swiss_roll,
    verbose = TRUE,
  )
plot(swiss_unrolled, col = phi)
```

## A Brief Description of LTSA

The mathematics and some of the descriptions of LTSA are not always very clear,
but this is what happens:

1. You have a dataset `X` with `N` items that you would like to reduce to `d`
dimensions.
1. Create an empty `N` x `N` matrix B (filled with zeros).
1. For each item in the dataset, define a neighborhood, i.e. the `k`-nearest
neighbors.
1. Calculate the SVD for just that neighborhood.
1. Create the `k` x `k` matrix `W = I - V.V'` where `V` is the right singular
vectors of the SVD (with some centering).
1. Update the elements of the `N` x `N` matrix `B` with the values of `W` by
adding the elements of `W` to the entries in `B` that correspond to the items in
the neighborhood.
1. Repeat for all items in the dataset.
1. `B` is now an un-normalized graph Laplacian matrix. Discard the smallest
eigenvector, and the next `d` eigenvectors are the coordinates of the reduced
dimension. This step is carried out like any other related spectral method, like
Laplacian eigenmaps or diffusion maps.

## Current Status

Very experimental. Will it ever graduate from 'Finicky' to 'Fast'? Its speed
relies on the interaction of:

* Using approximate nearest neighbors to define the local neighborhoods.
* Using a [sparse matrix](https://cran.r-project.org/package=Matrix)
representation of `B`.
* Using truncated SVD via [irlba](https://cran.r-project.org/package=irlba) to
calculate `W`.
* Using [RSpectra](https://cran.r-project.org/package=RSpectra) to efficiently
get the bottom eigenvectors of `B`.

It's not a great idea to use large values of `k` to define the neighborhoods:
you will get something that approaches a "global" SVD/PCA at much greater
effort. If you insist on doing this (e.g. set `n_neighbors = 1000`), the initial
sparse matrix allocation can be surprisingly slow.

The approximate nearest neighbor search (which can be exact if you want) *is*
parallelized, but that's the only thing that is.

The other bottleneck is the trucated SVD for each `W`. This *could* be done in
parallel, but I haven't got round to it yet. The better your [underlying linear
algebra library](https://csantill.github.io/RPerformanceWBLAS/), the better a
time you will have.

The final eigendecomposition of `B` can't be done in parallel, but RSpectra is
fast when it works. Please note that there a variety of failure states:

* RSpectra doesn't converge. You'll see an error if this happens. This can be
due to a bad choice of options, like the `tol` (the tolerance) or `ncv` (the
number of Lanczos basis vectors).
* RSpectra does converge, but provides results that look weird. This could
also be due to insufficiently tight convergence criteria. Some `B` matrices
have several eigenvectors with very very similar eigenvalues, which will then
be returned in the wrong order.
* RSpectra hangs, doing nothing, very slowly increasing its memory usage. This
seems to happen [during sparse
factorization](https://github.com/yixuan/spectra/issues/126). Unfortunately, I
don't have a solution for this.

If you use `method = "irlba"` or `method = "svdr"` then different functions in
irlba will be used instead of RSpectra. These are more likely to finish their
calculations successfully under circumstances where RSpectra stalls, but they
can require a lot more effort to give similarly converged results in other
scenarios. For the same amount of CPU time, the `svdr` setting may do the better
than `irlba` if you suspect there are several eigenvectors with very similar
eigenvalues and the ordering is incorrect (don't ask me how you would diagnose
that in general however).

## See Also

[Rdimtools](https://github.com/kisungyou/Rdimtools) also contains an R
implementation of LTSA.

## Further Reading

A fairly arbitrary selection of papers that contained at least one thing I found
interesting while writing this package.

Zhang, Z., & Zha, H. (2004). 
Principal manifolds and nonlinear dimensionality reduction via tangent space alignment. 
*SIAM journal on scientific computing*, *26*(1), 313-338.
<https://doi.org/10.1137/S1064827502419154>

Xiang, S., Nie, F., Pan, C., & Zhang, C. (2011). 
Regression reformulations of LLE and LTSA with locally linear transformation. 
*IEEE Transactions on Systems, Man, and Cybernetics, Part B (Cybernetics)*, *41*(5), 1250-1262.
<https://doi.org/10.1109/TSMCB.2011.2123886>

Sun, W., Halevy, A., Benedetto, J. J., Czaja, W., Li, W., Liu, C., ... & Wang, R. (2013). 
Nonlinear dimensionality reduction via the ENH-LTSA method for hyperspectral image classification.
*IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing*, *7*(2), 375-388.
<https://doi.org/10.1109/JSTARS.2013.2238890>

Hong, D., Yokoya, N., & Zhu, X. X. (2017). 
Learning a robust local manifold representation for hyperspectral dimensionality reduction. 
*IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing*, *10*(6), 2960-2975.
<https://doi.org/10.1109/JSTARS.2017.2682189>

Zhang, S., Ma, Z., & Tan, H. (2017).
On the Equivalence of HLLE and LTSA.
*IEEE transactions on cybernetics*, *48*(2), 742-753.
<https://doi.org/10.1109/TCYB.2017.2655338>

Ting, D., & Jordan, M. I. (2018).
On Nonlinear Dimensionality Reduction, Linear Smoothing and Autoencoding
*arXiv preprint* *arXiv*:1803.02432.
<https://arxiv.org/abs/1803.02432>

For some (perhaps under-appreciated) theoretical issues around spectral methods
for manifold learning, see:

Gerber, S., Tasdizen, T., & Whitaker, R. (2007, June). 
Robust non-linear dimensionality reduction using successive 1-dimensional Laplacian eigenmaps. 
In *Proceedings of the 24th international conference on Machine learning* (pp. 281-288).
<https://doi.org/10.1145/1273496.1273532>

Goldberg, Y., Zakai, A., Kushnir, D., & Ritov, Y. A. (2008). 
Manifold learning: The price of normalization.
*Journal of Machine Learning Research*, *9*(8).
<https://www.jmlr.org/papers/v9/goldberg08a.html>

Ting, D., & Jordan, M. I. (2020).
Manifold Learning via Manifold Deflation.
*arXiv preprint* *arXiv*:2007.03315.
<https://arxiv.org/abs/2007.03315>

The key paper remains the one by Goldberg and co-workers: if the "aspect ratio"
of the manifold has a variance that is very high in one axis over another (e.g.
a very long and thin manifold) then the normalization constraint inherent to
spectral methods like LTSA will fail to reproduce the manifold.

`flotsam` goes to various efforts to shift matrices to find *largest* eigenvalues
inspired by discussion between [Aaron Lun](https://github.com/LTLA), 
[Kevin Doherty](https://github.com/keevindoherty) and 
[Yixuan Qiu](https://github.com/yixuan) at [this RSpectra issue](https://github.com/yixuan/spectra/issues/126).
