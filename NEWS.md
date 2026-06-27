# flotsam 0.0.0.9002

* `ltsa()` now uses explicit fixed-width final eigenanalysis. New arguments
  `eig_k`, `output`, and `include_B` control candidate width and return shape.
  The default return remains the embedding matrix; `output = "B"` replaces the
  old matrix-return switch, and `output = "result"` returns compact
  eigenanalysis and assembly diagnostics. Runtime eigenanalysis no longer
  accepts adaptive rescue-policy controls; if diagnostics look suspicious, use
  a larger `eig_k` or stricter backend settings. Diagnostics are not
  completeness certificates.
* LTSA eigenanalysis documentation now explains the fixed-width `eig_k`
  request, null-vector projection, Rayleigh-Ritz postprocessing, compact
  diagnostics, backend-specific tuning settings, and the fact that
  `normalize = TRUE` is a separate normalized LTSA formulation rather than a
  rescue path for the unnormalized objective.

# flotsam 0.0.0.9001

* New parameter: `n_assembly_threads` to control number of threads to construct
  the `B` matrix.
* You can now pass a precomputed nearest-neighbor index (or the full dense
  neighbor graph list created by `rnndescent`) directly to `nn_method`, for
  example `ltsa(X, nn_method = nn$idx)`. Leave `n_neighbors` unset when using
  this because the correct number of neighbors is inferred from the graph.
* LTSA matrix construction is improved for high-dimensional inputs (`ncol(X)`
  larger than `n_neighbors`) by using a row-major centered Gram and a triangular
  matrix taking advantage of B matrix symmetry.
* New parameter: `copy_max_mib` to control the size cap for the optional
  row-major dense copy. The default is 256 MiB.
* Replaced the old sparse slot-search LTSA assembly path with the serial
  append/finalize C++ builder. Now modestly faster and slightly less peak
  memory usage.
* More robust eigenanalysis convergence.
* Iterative eigenanalysis now extracts more eigenvectors than requested and uses
  null-aware Rayleigh-Ritz polishing on the candidate subspace to attempt a
  better final result. `verbose = TRUE` reports residual, rank, and boundary-gap
  diagnostics. The `eig_method = "irlba"` and `"svdr"` paths now share this
  postprocessing, but rely on post-hoc diagnostics because they do not report
  RSpectra-style convergence counts.
* LTSA now warns when final eigenanalysis diagnostics indicate an ambiguous
  low-energy eigenspace, such as extra near-zero nonconstant modes or a weak
  boundary gap after the requested embedding dimensions. This can indicate a
  disconnected or weakly connected neighborhood graph, too-small `n_neighbors`,
  or `ndim` cutting through a low-energy eigenspace.
* If the `ndim` eigenvalues appear to contain only part of a near-zero
  low-energy cluster, a further refinement step is added with tighter settings
  and more candidates in case eigenvectors have been missed.
* Better argument validation and eigenanalysis error reporting.
* Updated GitHub Actions.

# flotsam 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
