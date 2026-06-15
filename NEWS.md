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
* Better argument validation and eigenanalysis error reporting.
* Updated GitHub Actions.

# flotsam 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
