#include "ltsa_internal.h"

LtsaTripletAssemblyBuilder::LtsaTripletAssemblyBuilder(
    const cpp11::integers& value_nnt, std::size_t value_n_nbrs,
    std::size_t n_obs, std::size_t max_int)
    : n_obs_(n_obs), value_n_nbrs_(value_n_nbrs), max_int_(max_int),
      canonical_columns_(n_obs), full_columns_(n_obs), append_seen_(n_obs, -1) {
  checked_triplet_count(n_obs_, value_n_nbrs_, "value_n_nbrs");

  std::vector<std::size_t> canonical_col_counts(n_obs_, 0);
  std::vector<std::size_t> full_col_counts(n_obs_, 0);
  std::vector<int> nni(value_n_nbrs_);
  std::vector<int> seen(n_obs_, -1);

  for (std::size_t obs = 0; obs < n_obs_; obs++) {
    std::size_t offset = obs * value_n_nbrs_;
    const int marker = static_cast<int>(obs);
    bool has_duplicate = false;
    for (std::size_t local = 0; local < value_n_nbrs_; local++) {
      const int idx = checked_neighbor_index(value_nnt[offset + local], n_obs_);
      nni[local] = idx;
      if (seen[idx] == marker) {
        has_duplicate = true;
      }
      seen[idx] = marker;
    }

    if (has_duplicate) {
      raw_entries_estimate_ += value_n_nbrs_ * value_n_nbrs_;
      for (std::size_t local_col = 0; local_col < value_n_nbrs_; local_col++) {
        full_col_counts[nni[local_col]] += value_n_nbrs_;
      }
    } else {
      raw_entries_estimate_ += value_n_nbrs_ * (value_n_nbrs_ + 1) / 2;
      for (std::size_t local_col = 0; local_col < value_n_nbrs_; local_col++) {
        for (std::size_t local_row = 0; local_row <= local_col; local_row++) {
          const int col = std::max(nni[local_row], nni[local_col]);
          canonical_col_counts[col]++;
        }
      }
    }
  }

  for (std::size_t col = 0; col < n_obs_; col++) {
    canonical_columns_[col].reserve(canonical_col_counts[col]);
    full_columns_[col].reserve(full_col_counts[col]);
  }
}

void LtsaTripletAssemblyBuilder::append_prechecked(
    const std::vector<int>& nni, const std::vector<double>& weights) {
  if (finalized_) {
    cpp11::stop("LTSA triplet builder has already been finalized");
  }
  if (n_appended_ >= n_obs_) {
    cpp11::stop("Too many LTSA neighborhoods appended");
  }
  if (nni.size() != value_n_nbrs_) {
    cpp11::stop("Inconsistent value neighborhood dimensions");
  }

  std::size_t value_k2 =
      checked_triplet_count(1, value_n_nbrs_, "value_n_nbrs");
  if (weights.size() != value_k2) {
    cpp11::stop("Inconsistent local weight dimensions");
  }

  if (has_duplicate_neighbors(nni)) {
    duplicate_fallback_count_++;
    append_full_prechecked(nni, weights);
  } else {
    append_triangular_prechecked(nni, weights);
  }

  n_appended_++;
}

SparseComponents LtsaTripletAssemblyBuilder::finalize_components() {
  if (finalized_) {
    cpp11::stop("LTSA triplet builder has already been finalized");
  }
  if (n_appended_ != n_obs_) {
    cpp11::stop("Not all LTSA neighborhoods were appended");
  }

  SparseComponents out;
  out.p.resize(n_obs_ + 1, 0);

  std::vector<double> row_sums(n_obs_, 0.0);
  std::vector<int> row_seen(n_obs_, -1);
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

std::size_t LtsaTripletAssemblyBuilder::raw_entries_estimate() const {
  return raw_entries_estimate_;
}

std::size_t LtsaTripletAssemblyBuilder::duplicate_fallback_count() const {
  return duplicate_fallback_count_;
}

bool LtsaTripletAssemblyBuilder::has_duplicate_neighbors(
    const std::vector<int>& nni) {
  const int marker = static_cast<int>(n_appended_ + 1);
  for (const int idx : nni) {
    if (append_seen_[idx] == marker) {
      return true;
    }
    append_seen_[idx] = marker;
  }
  return false;
}

void LtsaTripletAssemblyBuilder::append_full_prechecked(
    const std::vector<int>& nni, const std::vector<double>& weights) {
  for (std::size_t local_col = 0; local_col < value_n_nbrs_; local_col++) {
    const int col = nni[local_col];
    for (std::size_t local_row = 0; local_row < value_n_nbrs_; local_row++) {
      const int row = nni[local_row];
      const double value = weights[local_col * value_n_nbrs_ + local_row];
      full_columns_[col].push_back(CompactEntry{row, value});
    }
  }
}

void LtsaTripletAssemblyBuilder::append_triangular_prechecked(
    const std::vector<int>& nni, const std::vector<double>& weights) {
  for (std::size_t local_col = 0; local_col < value_n_nbrs_; local_col++) {
    for (std::size_t local_row = 0; local_row <= local_col; local_row++) {
      const int global_row = nni[local_row];
      const int global_col = nni[local_col];
      const int row = std::min(global_row, global_col);
      const int col = std::max(global_row, global_col);
      const double value = weights[local_col * value_n_nbrs_ + local_row];
      canonical_columns_[col].push_back(CompactEntry{row, value});
    }
  }
}

void LtsaTripletAssemblyBuilder::expand_canonical_to_full(
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
      full_columns_[col].push_back(CompactEntry{row, value});
      if (row != static_cast<int>(col)) {
        full_columns_[row].push_back(
            CompactEntry{static_cast<int>(col), value});
      }
    }
  }

  canonical_columns_.clear();
  canonical_columns_.shrink_to_fit();
}
