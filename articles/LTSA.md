# Local Tangent Space Alignment

## A Brief Description of LTSA

The mathematics and some of the descriptions of LTSA are not always very
clear, but this is what happens:

1.  You have a dataset `X` with `N` items that you would like to reduce
    to `d` dimensions.
2.  Create an empty `N` x `N` matrix B (filled with zeros).
3.  For each item in the dataset, define a neighborhood, i.e. the
    `k`-nearest neighbors.
4.  Center the `k` x `p` neighborhood matrix and calculate its SVD.
5.  Let `U` contain the retained **left** singular vectors, which give
    tangent coordinates for the `k` observations. Create the local
    residual projector
    $`W = I - \mathbf{1}\mathbf{1}^{\mathsf T}/k - UU^{\mathsf T}`$.
6.  Update the elements of the `N` x `N` matrix `B` with the values of
    `W` by adding the elements of `W` to the entries in `B` that
    correspond to the items in the neighborhood.
7.  Repeat for all items in the dataset.
8.  `B` is now a symmetric positive-semidefinite LTSA alignment operator
    with the constant vector as a known null direction. Remove that
    direction and select the lowest `d` nonconstant directions for the
    embedding. The nullspace need not be one-dimensional, so there is
    not necessarily a unique smallest eigenvector.

### Parameters

There’s only one function in `flotsam`: `ltsa`. Some parameters that may
be worth modifying in roughly descending order of importance:

- `ndim`: The number of dimensions to reduce to. Usually 2 for
  visualization.
- `n_neighbors`: The number of neighbors to use for the local
  neighborhood. When computed neighbors are requested and this is
  omitted, the default is the smaller of `15` and the maximum permitted
  by `N` and `include_self` (`N` with self or `N - 1` without self). It
  must still be at least `ndim + 2`; data too small to meet that minimum
  produce an error. A value that is too low can leave disconnected
  effective-neighborhood components or poorly separated low-energy
  directions.
- `nn_method`: The method to use for finding nearest neighbors:
  - `"nnd"`: Use approximate nearest neighbors by Nearest Neighbor
    Descent. This is the default and appropriate for large datasets.
  - `"exact"`: Use exact nearest neighbors by exhaustively comparing all
    items. Slow for large datasets.
  - A precomputed nearest-neighbor index matrix, or a nearest-neighbor
    result object with an `idx` matrix. With `include_self = TRUE`, the
    shape is `N` x `n_neighbors` and every row contains its own index.
    With `include_self = FALSE`, the shape is `N` x (`n_neighbors + 1`),
    self must be in the first column, and that column is removed before
    assembly.
- `include_self`: Whether to include the item itself as a neighbor.
  [Zhang and co-workers](https://doi.org/10.1109/JSTARS.2017.2682189)
  suggest that this is in effect the main difference between LTSA and
  the Hessian Locally Linear Embedding (HLLE) method, so setting this to
  `FALSE` may allow emulating HLLE. Default is `TRUE`.
- `n_threads`: The number of threads to use for the nearest neighbor
  calculation. Default is `1`. The rnndescent backend treats `0` and `1`
  as serial execution; values greater than `1` request multithreaded
  execution. Multithreaded approximate searches need not reproduce an
  identical neighbor graph across runs.
- `n_assembly_threads`: The number of threads to use when constructing
  the LTSA matrix `B`. Default is `1`, which requests serial assembly
  and generally uses less temporary storage. Values greater than `1`
  request parallel construction of `B`, which trades additional storage
  for speed and may compound threaded-BLAS oversubscription. See the
  [threading
  warning](https://jlmelville.github.io/flotsam/articles/numerical-diagnostics.html#threading).
- `eig_method`: The method to use for finding eigenvalues.
  - `"auto"`: Use dense
    [`base::eigen()`](https://rdrr.io/r/base/eigen.html) when
    `N <= dense_n` or `eig_k >= dense_fraction * N`, and otherwise use
    [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html).
    This is the default.
  - `"rspectra"`: Always use
    [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html).
  - `"irlba"`: Always use
    [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html) and
    post-hoc residual diagnostics.
  - `"svdr"`: Always use
    [`irlba::svdr()`](https://rdrr.io/pkg/irlba/man/svdr.html) and
    post-hoc residual diagnostics.
  - `"eigen"` or `"eig"`: Use
    [`base::eigen()`](https://rdrr.io/r/base/eigen.html). This converts
    `B` to a dense matrix and computes the full eigensystem, then
    retains the lowest `eig_k` candidate vectors for Rayleigh–Ritz
    postprocessing. In terms of both memory usage and CPU time this is
    only feasible for small datasets, but it is useful as a diagnostic
    reference.
- `eig_k`: The number of eigenvectors to return before a post-processing
  step and final `ndim` return. If `NULL`, the default is
  `min(N - 1L, max(12, ndim + 2))`. Must satisfy
  `ndim + 1 <= eig_k < n`. Asking for more eigenvectors than necessary
  seems to improve the reliability of the final result: RSpectra may
  otherwise miss eigenvectors that are close to zero, or fail to
  converge if the tolerance is tightened. The default seems sufficiently
  generous for most cases. Dense eigenanalysis still uses `eig_k` to
  limit the candidate span passed to the postprocessing step. The
  “Numerical Diagnostics” article has more details.
- `normalize`: Whether to solve the generalized problem
  $`Bv = \lambda Dv`$, where $`D = \operatorname{diag}(B)`$. This is a
  different problem from standard LTSA, with `D`-weighted orthogonality
  and centering constraints, rather than merely a numerically
  conditioned route to the same coordinates.
- `output`: The type of output to return:
  - `"embedding"`: Return the embedding matrix. This is the default.
  - `"result"`: Return a list containing the embedding, compact
    eigenanalysis diagnostics, assembly diagnostics, and optionally the
    assembled unnormalized LTSA matrix `B`.
  - `"B"`: Skip final eigenanalysis and return the raw alignment matrix
    when `normalize = FALSE`, or the normalized operator
    $`D^{-1/2} B D^{-1/2}`$ when `normalize = TRUE`.
- `include_B`: Whether to include the assembled unnormalized LTSA matrix
  `B` in the result. Under `eig_method = "auto"`, dense
  [`base::eigen()`](https://rdrr.io/r/base/eigen.html) is selected when
  `N <= 100` or `eig_k >= 0.5 * N`; `dense_n` and `dense_fraction`
  change those thresholds. Otherwise automatic selection uses RSpectra.
  Explicit `"rspectra"`, `"irlba"`, and `"svdr"` requests always use the
  named backend, including for small inputs. With `output = "result"`,
  `eigen$method` records `"auto"` or the explicit canonical method
  (`"eigen"` is stored as `"eig"`), while `eigen$backend$name` records
  the backend that actually ran. The package uses `tol = 1e-6` for
  RSpectra, rather than RSpectra’s own `1e-10` default.

`dense_n` and `dense_fraction` are accepted only with `"auto"`. The
curated backend controls are accepted only with their explicit backend,
and `shift_eps` only with an explicit iterative method. `resid_tol` and
`gap_tol` are non-routing diagnostics accepted in every mode; they are
the only controls accepted by `"eig"`/`"eigen"`. `output = "B"` rejects
all eigenanalysis controls in `...`.

For practical solver and troubleshooting guidance, see the [Numerical
Diagnostics](https://jlmelville.github.io/flotsam/articles/numerical-diagnostics.md)
article.
