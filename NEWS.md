# flotsam 0.0.0.9001

* LTSA matrix construction is improved for high-dimensional inputs by using
  row-major centered Gram (if the matrix isn't too large) and a triangular
  matrix taking advantage of B matrix symmetry.
* Replaced the old sparse slot-search LTSA assembly path with the serial
  append/finalize C++ builder. Now modestly faster and slightly less peak
  memory usage.
* Better argument validation and eigenanalysis error reporting.
* Updated GitHub Actions.

# flotsam 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
