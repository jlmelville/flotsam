#include "ltsa_internal.h"

std::size_t checked_size_add(std::size_t lhs, std::size_t rhs,
                             const char* message) {
  if (lhs > std::numeric_limits<std::size_t>::max() - rhs) {
    stop("%s", message);
  }
  return lhs + rhs;
}

std::size_t checked_size_mul(std::size_t lhs, std::size_t rhs,
                             const char* message) {
  if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
    stop("%s", message);
  }
  return lhs * rhs;
}

std::size_t triangular_pair_count(std::size_t n_nbrs) {
  if (n_nbrs > (std::numeric_limits<std::size_t>::max() - 1) / n_nbrs) {
    stop("Too many triangular LTSA contributions to assemble");
  }
  return n_nbrs * (n_nbrs + 1) / 2;
}

std::size_t triangular_pair_offset(std::size_t local_col,
                                   std::size_t local_row) {
  return local_col * (local_col + 1) / 2 + local_row;
}

std::size_t checked_raw_staging_bytes(std::size_t canonical_count,
                                      std::size_t full_count) {
  const std::size_t total_count = checked_size_add(
      canonical_count, full_count, "Too many raw LTSA contributions to stage");
  const std::size_t row_bytes = checked_size_mul(
      total_count, sizeof(int), "Raw LTSA row staging buffer is too large");
  const std::size_t value_bytes =
      checked_size_mul(total_count, sizeof(double),
                       "Raw LTSA value staging buffer is too large");
  return checked_size_add(row_bytes, value_bytes,
                          "Raw LTSA staging buffers are too large");
}

std::string row_major_fallback_reason(bool use_gram_workspace,
                                      bool row_major_used,
                                      bool row_major_within_limit) {
  if (!use_gram_workspace) {
    return "not_applicable_svd_route";
  }
  if (row_major_used) {
    return "";
  }
  return row_major_within_limit ? "allocation_failed"
                                : "copy_size_exceeds_limit";
}
