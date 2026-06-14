# Local Tangent Space Alignment

Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha,
2004) for dimensionality reduction.

## Usage

``` r
ltsa(
  X,
  n_neighbors = 15,
  ndim = 2,
  nn_method = "nnd",
  eig_method = "rspectra",
  include_self = TRUE,
  normalize = FALSE,
  ret_B = FALSE,
  n_threads = 0,
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  The input data matrix or dataframe with one observation per row.

- n_neighbors:

  The size of local neighborhood (in terms of number of neighboring
  sample points) used for manifold approximation.

- ndim:

  The dimension of the space to embed into.

- nn_method:

  Method for finding nearest neighbors. Can be one of:

  - `"nnd"` Approximate nearest neighbors by Nearest Neighbor Descent.

  - `"exact"` Exact nearest neighbors by exhaustively comparing all
    items. Slow for large datasets.

- eig_method:

  How to carry out the eigendecomposition. Possible values are:

  - `"rspectra"` Use
    [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html).

  - `"irlba"` Use
    [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html).

  - `"svdr"` Use
    [`irlba::svdr()`](https://rdrr.io/pkg/irlba/man/svdr.html).

  - `"fullsvd"` Use the [`base::svd()`](https://rdrr.io/r/base/svd.html)
    function. This is only feasible for small datasets and should be
    used for diagnostic purposes only.

  - `"eig"` or `"eigen"` Use the
    [`base::eigen()`](https://rdrr.io/r/base/eigen.html) function. This
    is only feasible for small datasets and should be used for
    diagnostic purposes only.

- include_self:

  Should an item be part of its own neighborhood? This has a minor
  effect on most results, but work by Zhang and co-workers (2017)
  suggests that this is in effect the main difference between LTSA and
  the Hessian Locally Linear Embedding (HLLE) method, so setting this to
  `FALSE` may allow emulating the HLLE method.

- normalize:

  If `TRUE` calculate the eigendecomposition on a normalized version of
  the Laplacian. This may be slightly easier to converge while giving
  similar results to the un-normalized case. It may also have suitable
  better properties if clustering is to be carried out on the
  eigenvectors.

- ret_B:

  If `TRUE`, return the matrix instead of the eigenvectors. This is
  mainly useful for diagnostic purposes if eigendecomposition is
  failing.

- n_threads:

  Number of threads to use. Applies only to the nearest neighbor
  calculation.

- verbose:

  If `TRUE` log information about progress to the console.

- ...:

  Extra arguments to be passed to the eigendecomposition method
  specified by `eig_method`. For `"rspectra"`, arguments are passed to
  the `opts` list. Suitable parameters include:

  - `ncv` Number of Lanzcos vectors to use.

  - `tol` Tolerance.

  - `maxitr` Maximum number of iterations.

  For `"irlba"` suitable arguments are:

  - `work` Working subspace dimension size.

  - `tol` Tolerance.

  - `maxit` Maximum number of iterations.

  For `"svdr"` suitable arguments are:

  - `extra` Number of extra vectors to use.

  - `tol` Tolerance.

  - `it` Maximum number of iterations.

  For more details see the documentation for
  [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
  [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html) and
  [`irlba::svdr()`](https://rdrr.io/pkg/irlba/man/svdr.html) functions,
  respectively. Don't pass other arguments unless you know what you are
  doing, as it may cause the `ltsa` to fail.

## References

Zhang, Z., & Zha, H. (2004). Principal manifolds and nonlinear
dimensionality reduction via tangent space alignment. *SIAM journal on
scientific computing*, *26*(1), 313-338.
<https://doi.org/10.1137/S1064827502419154>

Zhang, S., Ma, Z., & Tan, H. (2017). On the Equivalence of HLLE and
LTSA. *IEEE transactions on cybernetics*, *48*(2), 742-753.
<https://doi.org/10.1109/TCYB.2017.2655338>

## Examples

``` r
n <- 1000
max_z <- 10

phi <- stats::runif(n, min = 1.5 * pi, max = 4.5 * pi)
x <- phi * cos(phi)
y <- phi * sin(phi)
z <- stats::runif(n, max = max_z)
swiss_roll <- data.frame(x, y, z)

# unroll it
swiss_ltsa <- ltsa(swiss_roll)
plot(swiss_ltsa, col = phi)


# compare with PCA
swiss_pca <- stats::prcomp(swiss_roll, rank. = 2, scale = FALSE, retx = TRUE)$x
plot(swiss_pca, col = phi)
```
