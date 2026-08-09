#include "ltsa_internal.h"

std::size_t checked_triplet_count(std::size_t n_obs, std::size_t n_nbrs,
                                  const char* name) {
  if (n_nbrs == 0) {
    cpp11::stop("%s must be positive", name);
  }
  const auto max_size = std::numeric_limits<std::size_t>::max();
  if (n_obs > max_size / n_nbrs || n_obs * n_nbrs > max_size / n_nbrs) {
    cpp11::stop("Too many LTSA triplets to assemble");
  }
  return n_obs * n_nbrs * n_nbrs;
}

int checked_zero_based_neighbor_index(int idx, std::size_t n_obs) {
  if (idx < 1 || static_cast<std::size_t>(idx) > n_obs) {
    cpp11::stop("Neighborhood index is outside the sparse matrix dimensions");
  }
  return idx - 1;
}

void checked_append_output(int row, double value, std::vector<int>& out_i,
                           std::vector<double>& out_x, std::size_t max_int) {
  const std::size_t max_slots =
      std::min(max_int, std::min(out_i.max_size(), out_x.max_size()));
  if (out_i.size() >= max_slots) {
    cpp11::stop("Too many non-zero slots for a dgCMatrix");
  }
  out_i.push_back(row);
  out_x.push_back(value);
}

void checked_ndim(int ndim) {
  if (ndim < 1) {
    cpp11::stop("ndim must be positive");
  }
}

int checked_lapack_dim(std::size_t value, const char* name) {
  const auto max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (value > max_int) {
    cpp11::stop("%s is too large for LAPACK", name);
  }
  return static_cast<int>(value);
}
