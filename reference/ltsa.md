# Local Tangent Space Alignment

Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha,
2004) for dimensionality reduction.

## Usage

``` r
ltsa(
  X,
  n_neighbors = NULL,
  ndim = 2,
  nn_method = "nnd",
  eig_method = "rspectra",
  include_self = TRUE,
  normalize = FALSE,
  ret_B = FALSE,
  n_threads = 0,
  n_assembly_threads = 1,
  copy_max_mib = 256,
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  The input data matrix or data frame with one observation per row. If a
  data frame is supplied, non-numeric columns are ignored. At least one
  numeric column is required.

- n_neighbors:

  The size of local neighborhood (in terms of number of neighboring
  sample points) used for manifold approximation. If `NULL`, the default
  is `15` when neighbors are computed, or inferred when a precomputed
  graph is supplied as `nn_method`.

- ndim:

  The dimension of the space to embed into.

- nn_method:

  Method for finding nearest neighbors, or a precomputed
  nearest-neighbor graph. Can be one of:

  - `"nnd"` Approximate nearest neighbors by Nearest Neighbor Descent.

  - `"exact"` Exact nearest neighbors by exhaustively comparing all
    items. Slow for large datasets.

  - A precomputed nearest-neighbor index matrix.

  - A nearest-neighbor result object with an `idx` matrix.

  Precomputed index matrices should match the raw neighbor-search output
  used by `rnndescent`: one row per observation and 1-based indices into
  `X`. When `include_self = TRUE`, `ncol(nn_method)` must equal
  `n_neighbors` and each row must contain its own row index. When
  `include_self = FALSE`, `ncol(nn_method)` must equal `n_neighbors + 1`
  and the first column must contain the row index, because that column
  is dropped before LTSA assembly.

- eig_method:

  How to carry out the final eigendecomposition. Possible values are:

  - `"rspectra"` Use
    [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html)
    to generate candidate vectors, then apply null-aware Rayleigh-Ritz
    postprocessing. This is the default and can report hard backend
    convergence failures through RSpectra's `nconv` metadata.

  - `"irlba"` Use
    [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html) to
    generate candidate vectors, then apply the same null-aware
    Rayleigh-Ritz postprocessing. This backend does not report
    RSpectra-style convergence counts, so convergence issues are
    detected with post-hoc residual diagnostics.

  - `"svdr"` Use
    [`irlba::svdr()`](https://rdrr.io/pkg/irlba/man/svdr.html) to
    generate candidate vectors, then apply the same null-aware
    Rayleigh-Ritz postprocessing. This backend does not report
    RSpectra-style convergence counts, so convergence issues are
    detected with post-hoc residual diagnostics.

  - `"fullsvd"` Use the [`base::svd()`](https://rdrr.io/r/base/svd.html)
    function. This is only feasible for small datasets and should be
    used for diagnostic purposes only.

  - `"eig"` or `"eigen"` Use the
    [`base::eigen()`](https://rdrr.io/r/base/eigen.html) function. This
    is only feasible for small datasets and should be used for
    diagnostic purposes only. Dense `"eig"` is the better diagnostic
    reference when algebraic eigenvalue ordering matters.

- include_self:

  Should an item be part of its own neighborhood? This has a minor
  effect on most results, but work by Zhang and co-workers (2017)
  suggests that this is in effect the main difference between LTSA and
  the Hessian Locally Linear Embedding (HLLE) method, so setting this to
  `FALSE` may allow emulating the HLLE method.

- normalize:

  If `TRUE`, calculate the final decomposition on a normalized version
  of the LTSA alignment matrix. For iterative methods (`"rspectra"`,
  `"irlba"`, and `"svdr"`), the normalized operator uses the same
  null-aware Rayleigh-Ritz selection as the unnormalized path, then
  transforms the selected vectors back to output coordinates. This may
  be slightly easier to converge while giving similar results to the
  unnormalized case. It may also have better properties if clustering is
  to be carried out on the eigenvectors.

- ret_B:

  If `TRUE`, return the matrix instead of the eigenvectors. This is
  mainly useful for diagnostic purposes if eigendecomposition is
  failing.

- n_threads:

  Number of threads to use. Applies only to the nearest neighbor
  calculation.

- n_assembly_threads:

  Number of threads to use when constructing the LTSA alignment matrix
  `B` after nearest neighbors are computed. The default `1` preserves
  the serial assembly path. Values greater than `1` opt into parallel
  construction of `B`, which can be faster but may increase peak memory
  use.

- copy_max_mib:

  Maximum size, in MiB, of the optional row-major dense copy of `X` used
  to speed up high-dimensional local Gram assembly. If this cap is
  exceeded, no row-major copy is made. The default is 256 MiB. Set to
  `0` to disable this copy. This code path is only used when the number
  of numeric columns in `X` is greater than the number of neighbors.
  Outside of synthetic datasets, you are quite likely to hit this code
  path.

- verbose:

  If `TRUE` log information about progress to the console.

- ...:

  Extra arguments to be passed to the eigendecomposition method
  specified by `eig_method`. For `"rspectra"`, arguments are passed to
  the `opts` list used by
  [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html).
  Suitable parameters include:

  - `ncv` Number of Lanczos vectors to use.

  - `tol` Tolerance.

  - `maxitr` Maximum number of iterations.

  For `"irlba"`, arguments are passed to
  [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html). Suitable
  arguments include:

  - `work` Working subspace dimension size.

  - `tol` Tolerance.

  - `maxit` Maximum number of iterations.

  - `reorth` Whether to use full reorthogonalization.

  For `"svdr"`, arguments are passed to
  [`irlba::svdr()`](https://rdrr.io/pkg/irlba/man/svdr.html). Suitable
  arguments include:

  - `extra` Number of extra vectors to use.

  - `tol` Tolerance.

  - `it` Maximum number of iterations.

  The iterative methods all use shared Rayleigh-Ritz postprocessing,
  residual diagnostics, and adaptive candidate expansion. However, only
  RSpectra reports explicit convergence counts; `"irlba"` and `"svdr"`
  rely on post-hoc residual diagnostics instead.

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
#> Warning: LTSA eigenanalysis found an ambiguous low-energy eigenspace; embedding may not be unique up to only rotation/sign. Diagnostics: the boundary gap after ndim is weak: global gap = 8.795e-06, local gap = 0.8579, tolerance = 1e-04. This can happen when the neighborhood graph is disconnected or weakly connected, n_neighbors is too small, or ndim cuts through a low-energy eigenspace.
plot(swiss_ltsa, col = phi)


# compare with PCA
swiss_pca <- stats::prcomp(swiss_roll, rank. = 2, scale = FALSE, retx = TRUE)$x
plot(swiss_pca, col = phi)
```
