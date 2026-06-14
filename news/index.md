# Changelog

## flotsam 0.0.0.9001

- Replaced the old sparse slot-search LTSA assembly path with the serial
  append/finalize C++ builder as the only internal assembly
  implementation. Now modestly faster and slightly less peak memory
  usage.
- Better argument validation and eigenanalysis error reporting.
- Updated GitHub Actions.

## flotsam 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
