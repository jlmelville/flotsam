#include "ltsa_internal.h"

TripletAssemblyBuilder::TripletAssemblyBuilder(
    const cpp11::integers& transposed_neighbor_indices, std::size_t n_neighbors,
    std::size_t n_obs, std::size_t max_int)
    : n_obs_(n_obs), n_neighbors_(n_neighbors), max_int_(max_int),
      canonical_columns_(checked_vector_size<std::vector<CompactEntry>>(
          n_obs, "serial LTSA canonical column containers")),
      full_columns_(checked_vector_size<std::vector<CompactEntry>>(
          n_obs, "serial LTSA full column containers")) {
  checked_triplet_count(n_obs_, n_neighbors_, "value_n_nbrs");
  checked_size_mul(n_obs_, triangular_pair_count(n_neighbors_),
                   "Too many triangular LTSA contributions to stage");

  std::vector<std::size_t> canonical_col_counts(
      checked_vector_size<std::size_t>(n_obs_,
                                       "serial LTSA canonical column counts"),
      0);
  std::vector<int> neighbor_indices(checked_vector_size<int>(
      n_neighbors_, "serial LTSA neighborhood indices"));

  for (std::size_t obs = 0; obs < n_obs_; obs++) {
    std::size_t offset = obs * n_neighbors_;
    for (std::size_t local = 0; local < n_neighbors_; local++) {
      const int idx = checked_zero_based_neighbor_index(
          transposed_neighbor_indices[offset + local], n_obs_);
      neighbor_indices[local] = idx;
    }

    for (std::size_t local_col = 0; local_col < n_neighbors_; local_col++) {
      for (std::size_t local_row = 0; local_row <= local_col; local_row++) {
        const int col =
            std::max(neighbor_indices[local_row], neighbor_indices[local_col]);
        canonical_col_counts[col]++;
      }
    }
  }

  for (std::size_t col = 0; col < n_obs_; col++) {
    canonical_columns_[col].reserve(checked_vector_size<CompactEntry>(
        canonical_col_counts[col], "serial LTSA canonical column"));
  }
}

void TripletAssemblyBuilder::append_prechecked_neighborhood(
    const std::vector<int>& neighbor_indices,
    const std::vector<double>& weights) {
  if (finalized_) {
    cpp11::stop("LTSA triplet builder has already been finalized");
  }
  if (n_appended_ >= n_obs_) {
    cpp11::stop("Too many LTSA neighborhoods appended");
  }
  if (neighbor_indices.size() != n_neighbors_) {
    cpp11::stop("Inconsistent value neighborhood dimensions");
  }

  std::size_t value_k2 = checked_triplet_count(1, n_neighbors_, "value_n_nbrs");
  if (weights.size() != value_k2) {
    cpp11::stop("Inconsistent local weight dimensions");
  }

  append_triangular_prechecked(neighbor_indices, weights);

  n_appended_++;
}

SparseComponents TripletAssemblyBuilder::finalize_sparse_components() {
  if (finalized_) {
    cpp11::stop("LTSA triplet builder has already been finalized");
  }
  if (n_appended_ != n_obs_) {
    cpp11::stop("Not all LTSA neighborhoods were appended");
  }

  SparseComponents out;
  out.p.resize(
      checked_vector_size<int>(
          checked_size_add(n_obs_, 1, "Too many serial LTSA sparse columns"),
          "serial LTSA sparse column pointers"),
      0);

  std::vector<double> row_sums(
      checked_vector_size<double>(n_obs_, "serial LTSA row sums"), 0.0);
  std::vector<int> row_seen(
      checked_vector_size<int>(n_obs_, "serial LTSA row markers"), -1);
  std::vector<int> touched_rows;

  expand_canonical_to_full(row_sums, row_seen, touched_rows);
  std::fill(row_seen.begin(), row_seen.end(), -1);

  for (std::size_t col = 0; col < n_obs_; col++) {
    auto& entries = full_columns_[col];
    int col_marker = static_cast<int>(col);
    touched_rows.clear();
    for (const auto& entry : entries) {
      if (row_seen[entry.row] != col_marker) {
        row_seen[entry.row] = col_marker;
        row_sums[entry.row] = 0.0;
        touched_rows.push_back(entry.row);
      }
      row_sums[entry.row] += entry.value;
    }

    std::sort(touched_rows.begin(), touched_rows.end());
    for (const auto& row : touched_rows) {
      double value = row_sums[row];
      if (value != 0.0) {
        checked_append_output(row, value, out.i, out.x, max_int_);
      }
    }
    out.p[col + 1] = static_cast<int>(out.i.size());
  }

  finalized_ = true;
  canonical_columns_.clear();
  canonical_columns_.shrink_to_fit();
  full_columns_.clear();
  full_columns_.shrink_to_fit();

  return out;
}

void TripletAssemblyBuilder::append_triangular_prechecked(
    const std::vector<int>& neighbor_indices,
    const std::vector<double>& weights) {
  for (std::size_t local_col = 0; local_col < n_neighbors_; local_col++) {
    for (std::size_t local_row = 0; local_row <= local_col; local_row++) {
      const int global_row = neighbor_indices[local_row];
      const int global_col = neighbor_indices[local_col];
      const int row = std::min(global_row, global_col);
      const int col = std::max(global_row, global_col);
      const double value = weights[local_col * n_neighbors_ + local_row];
      canonical_columns_[col].push_back(CompactEntry{row, value});
    }
  }
}

void TripletAssemblyBuilder::expand_canonical_to_full(
    std::vector<double>& row_sums, std::vector<int>& row_seen,
    std::vector<int>& touched_rows) {
  for (std::size_t col = 0; col < n_obs_; col++) {
    auto& entries = canonical_columns_[col];
    const int col_marker = static_cast<int>(col);
    touched_rows.clear();

    for (const auto& entry : entries) {
      if (row_seen[entry.row] != col_marker) {
        row_seen[entry.row] = col_marker;
        row_sums[entry.row] = 0.0;
        touched_rows.push_back(entry.row);
      }
      row_sums[entry.row] += entry.value;
    }

    for (const int row : touched_rows) {
      const double value = row_sums[row];
      if (value == 0.0) {
        continue;
      }
      if (full_columns_[col].size() >= full_columns_[col].max_size()) {
        cpp11::stop("serial LTSA full compact column is too large");
      }
      full_columns_[col].push_back(CompactEntry{row, value});
      if (row != static_cast<int>(col)) {
        if (full_columns_[row].size() >= full_columns_[row].max_size()) {
          cpp11::stop("serial LTSA full compact column is too large");
        }
        full_columns_[row].push_back(
            CompactEntry{static_cast<int>(col), value});
      }
    }
  }

  canonical_columns_.clear();
  canonical_columns_.shrink_to_fit();
}
