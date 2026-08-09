#include "ltsa_internal.h"

std::size_t checked_size_add(std::size_t lhs, std::size_t rhs,
                             const char* message) {
  if (lhs > std::numeric_limits<std::size_t>::max() - rhs) {
    cpp11::stop("%s", message);
  }
  return lhs + rhs;
}

std::size_t checked_size_mul(std::size_t lhs, std::size_t rhs,
                             const char* message) {
  if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
    cpp11::stop("%s", message);
  }
  return lhs * rhs;
}

std::size_t triangular_pair_count(std::size_t n_nbrs) {
  if (n_nbrs > (std::numeric_limits<std::size_t>::max() - 1) / n_nbrs) {
    cpp11::stop("Too many triangular LTSA contributions to assemble");
  }
  return n_nbrs * (n_nbrs + 1) / 2;
}

std::size_t triangular_pair_offset(std::size_t local_col,
                                   std::size_t local_row) {
  return local_col * (local_col + 1) / 2 + local_row;
}
