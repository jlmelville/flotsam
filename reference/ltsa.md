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
  eig_k = NULL,
  output = c("embedding", "result", "B"),
  include_B = FALSE,
  include_self = TRUE,
  normalize = FALSE,
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

- eig_k:

  Explicit fixed-width number of candidate vectors requested from the
  final eigensolver. If `NULL`, the default is
  `min(n - 1L, max(12L, ndim + 2L))`, where `n` is the number of
  observations. Must satisfy `ndim + 1 <= eig_k < n`. Larger values give
  the Rayleigh-Ritz postprocessing a wider candidate span.

- output:

  What to return:

  - `"embedding"` Return the embedding matrix. This is the default.

  - `"result"` Return a list containing the embedding, compact
    eigenanalysis diagnostics, assembly diagnostics, and optionally `B`.

  - `"B"` Return the assembled unnormalized LTSA matrix and skip final
    eigenanalysis.

- include_B:

  If `TRUE` and `output = "result"`, include the assembled unnormalized
  LTSA matrix `B` in the result object. Ignored for other output modes.

- include_self:

  Should an item be part of its own neighborhood? This has a minor
  effect on most results, but work by Zhang and co-workers (2017)
  suggests that this is in effect the main difference between LTSA and
  the Hessian Locally Linear Embedding (HLLE) method, so setting this to
  `FALSE` may allow emulating the HLLE method.

- normalize:

  If `TRUE`, calculate the final decomposition on a normalized LTSA
  formulation rather than the unnormalized alignment matrix. This is a
  separate spectral objective. For iterative methods (`"rspectra"`,
  `"irlba"`, and `"svdr"`), the normalized operator uses the same
  fixed-width null-aware Rayleigh-Ritz selection as the unnormalized
  path, then transforms the selected vectors back to output coordinates.
  The normalized formulation may have different downstream properties,
  for example when clustering on the eigenvectors.

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
  specified by `eig_method`. These backend settings affect how the
  requested fixed-width candidate block is computed. For `"rspectra"`,
  arguments are passed to the `opts` list used by
  [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html).
  Common tuning parameters are:

  - `tol` Tolerance.

  - `maxitr` Maximum number of iterations.

  - `ncv` Number of Lanczos vectors to use.

  For `"irlba"`, arguments are passed to
  [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html). Common
  tuning parameters are:

  - `tol` Tolerance.

  - `maxit` Maximum number of iterations.

  - `reorth` Whether to use full reorthogonalization.

  For `"svdr"`, arguments are passed to
  [`irlba::svdr()`](https://rdrr.io/pkg/irlba/man/svdr.html). Common
  tuning parameters are:

  - `tol` Tolerance.

  - `it` Maximum number of iterations.

  The iterative methods all use shared fixed-width Rayleigh-Ritz
  postprocessing and residual diagnostics. Only RSpectra reports
  explicit convergence counts; `"irlba"` and `"svdr"` rely on post-hoc
  residual diagnostics instead. `resid_tol` and `gap_tol` set the scaled
  residual and global boundary-gap diagnostic tolerances used to
  classify the requested solve.

  For more details see the documentation for
  [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
  [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html) and
  [`irlba::svdr()`](https://rdrr.io/pkg/irlba/man/svdr.html) functions,
  respectively. Don't pass other arguments unless you know what you are
  doing, as it may cause the `ltsa` to fail.

## Details

LTSA forms a sparse alignment matrix `B` and then returns the low-energy
nonconstant eigenvectors of that matrix. The public iterative final
eigenanalysis path is a fixed-width solve: `eig_k` is the number of
candidate vectors requested from the selected backend before
Rayleigh-Ritz selection and diagnostic classification.

After the backend returns its candidate block, `flotsam` projects out
the known constant null vector, orthonormalizes the remaining candidate
span, performs a small Rayleigh-Ritz extraction on the operator being
solved, sorts by the resulting Ritz values, and returns the first `ndim`
nonconstant Ritz vectors. This postprocessing is intended to make
clustered low-energy spectra less fragile than returning the raw backend
vectors.

Use `output = "result"` to inspect compact diagnostics for the requested
solve. The `eigen` component records the backend method, whether the
normalized formulation was used, `eig_k`, selected values, all post-null
Ritz values available inside the candidate span, scaled residuals,
post-null rank, an estimated largest eigenvalue, solver-neutral status,
messages, and backend metadata. The `assembly` component records
nearest-neighbor source and timing alongside matrix assembly
diagnostics. These diagnostics classify the requested fixed-width solve,
but they are not completeness certificates for the full low-energy
eigenspace. If diagnostics look suspicious, rerun with a larger `eig_k`
or stricter backend settings.

The `"rspectra"` path first estimates the largest eigenvalue and solves
a shifted largest-algebraic problem instead of using shift-invert near
zero. RSpectra's `nconv` metadata is treated as a hard backend
convergence check. The `"irlba"` and `"svdr"` paths use the same
Rayleigh-Ritz selection and diagnostics, but they rely on post-hoc
residual checks because those backends do not expose RSpectra-style
convergence counts.

`normalize = TRUE` selects a separate normalized LTSA formulation. For
iterative methods, the normalized operator uses the same fixed-width,
null-aware Rayleigh-Ritz workflow and then transforms the selected
vectors back to output coordinates.

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


# The default return is the embedding matrix.
small_iris <- iris[1:75, 1:4]
iris_ltsa <- ltsa(
  small_iris,
  n_neighbors = 12,
  nn_method = "exact"
)

# Request compact eigenanalysis and assembly diagnostics.
iris_result <- ltsa(
  small_iris,
  n_neighbors = 12,
  nn_method = "exact",
  eig_k = 8,
  output = "result"
)
iris_result$eigen$status
#> [1] "warning"
iris_result$eigen$messages
#> [1] "Only part of a near-zero selected Ritz block is present."     
#> [2] "Increasing eig_k or using stricter backend settings may help."

# Return the assembled unnormalized LTSA matrix without final eigenanalysis.
iris_B <- ltsa(
  small_iris,
  n_neighbors = 12,
  nn_method = "exact",
  output = "B"
)
dim(iris_B)
#> [1] 75 75

# If diagnostics are suspicious, rerun with a wider fixed candidate request.
iris_wider <- ltsa(
  small_iris,
  n_neighbors = 12,
  nn_method = "exact",
  eig_k = 16,
  output = "result"
)

# Or rerun with stricter RSpectra backend settings.
iris_strict <- ltsa(
  small_iris,
  n_neighbors = 12,
  nn_method = "exact",
  eig_k = 8,
  output = "result",
  tol = 1e-8,
  maxitr = 5000,
  ncv = 30
)
```
