#include "ltsa_internal.h"

// calculates M * d where M is a sparse matrix with pointer ps and values
// xs, and d is the vector ds
[[cpp11::register]] cpp11::doubles
spm_times_scalar(const cpp11::integers& ps, const cpp11::doubles& xs,
                 const cpp11::doubles& ds) {
  auto nrow = ds.size();
  if (nrow != ps.size() - 1) {
    cpp11::stop("Inconsistent diagonal and pointer lengths");
  }
  if (xs.size() != ps[nrow]) {
    cpp11::stop("Inconsistent value and pointers");
  }
  cpp11::writable::doubles dxs(xs.size());

  for (R_xlen_t i = 0; i < nrow; i++) {
    R_xlen_t begin = ps[i];
    R_xlen_t end = ps[i + 1];
    for (R_xlen_t j = begin; j < end; j++) {
      dxs[j] = xs[j] * ds[i];
    }
  }

  return dxs;
}
