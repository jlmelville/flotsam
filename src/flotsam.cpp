#include <algorithm>
#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <limits>
#include <unordered_set>
#include <vector>
using namespace cpp11;

// convert the row/column indices of neighbor idx in ns to the 1D index in the
// sparse is vector, using ps to find the row range
[[cpp11::register]] integers
sparse_idxs(const integers& is, const integers& ps, const integers& ns) {
  // contents of ns are 1-indexed
  // contents of is and ps are 0-indexed

  std::size_t n_nbrs = ns.size();
  writable::integers sp_idx(n_nbrs * n_nbrs);
  auto mbegin = is.cbegin();
  R_xlen_t n_cols = ps.size() - 1;

  // loop over every pair in n including e.g. (j, i) as well as (i, j)
  for (std::size_t i = 0; i < n_nbrs; i++) {
    auto r = ns[i]; // r is 1-indexed
    if (r < 1 || r > n_cols) {
      stop("Neighborhood index is outside the sparse matrix dimensions");
    }

    for (std::size_t j = 0; j < n_nbrs; j++) {
      auto c = ns[j]; // c is 1-indexed
      if (c < 1 || c > n_cols) {
        stop("Neighborhood index is outside the sparse matrix dimensions");
      }

      // cbegin/cend range of data at row r

      R_xlen_t cbegin = ps[r - 1];
      R_xlen_t cend = ps[r];
      if (cbegin < 0 || cend < cbegin || cend > is.size()) {
        stop("Inconsistent sparse matrix pointers");
      }

      auto range_begin = integers::const_iterator(mbegin + cbegin);
      auto range_end = integers::const_iterator(mbegin + cend);
      auto pos = std::find(range_begin, range_end, c - 1);
      if (pos == range_end) {
        stop("Sparse matrix does not contain requested neighborhood pair");
      }
      sp_idx[i * n_nbrs + j] = pos - mbegin + 1;
    }
  }

  return sp_idx;
}

  // calculates M * d where M is a sparse matrix with pointer ps and values
  // xs, and d is the vector ds
  [[cpp11::register]] doubles
  spm_times_scalar(const integers& ps, const doubles& xs, const doubles& ds)
{
  auto nrow = ds.size();
  if (nrow != ps.size() - 1) {
    stop("Inconsistent diagonal and pointer lengths");
  }
  if (xs.size() != ps[nrow]) {
    stop("Inconsistent value and pointers");
  }
  writable::doubles dxs(xs.size());

  for (R_xlen_t i = 0; i < nrow; i++) {
    R_xlen_t begin = ps[i];
    R_xlen_t end = ps[i + 1];
    for (R_xlen_t j = begin; j < end; j++) {
      dxs[j] = xs[j] * ds[i];
    }
  }

  return dxs;
}

std::vector<std::size_t>
create_idx1d(const integers& nnt, std::size_t n_nbrs, std::size_t n_obs)
{
  std::unordered_set<std::size_t> all_idx;
  std::unordered_set<std::size_t> tmp;
  for (std::size_t i = 0; i < n_obs; i++) {
    auto inbrs = i * n_nbrs;
    for (std::size_t j = 0; j < n_nbrs; j++) {
      auto nbrj = nnt[inbrs + j] - 1;
      for (std::size_t k = j; k < n_nbrs; k++) {
        auto nbrk = nnt[inbrs + k] - 1;
        auto idx = nbrj < nbrk ? nbrj * n_obs + nbrk : nbrk * n_obs + nbrj;
        tmp.insert(idx);
      }
    }
    if (tmp.size() > all_idx.size()) {
      all_idx.insert(tmp.begin(), tmp.end());
      tmp.clear();
    }
  }
  if (tmp.size() > 0) {
    all_idx.insert(tmp.begin(), tmp.end());
    tmp.clear();
  }

  // create sorted vector from set
  std::vector<std::size_t> ridx;
  ridx.reserve(all_idx.size());
  ridx.insert(ridx.end(), all_idx.begin(), all_idx.end());

  return ridx;
}

[[cpp11::register]] list
nbrhood_triplets(const integers& nnt, std::size_t n_nbrs) {
  auto n_obs = nnt.size() / n_nbrs;

  // find unique 1D indices (i < j only)
  std::vector<std::size_t> ridx = create_idx1d(nnt, n_nbrs, n_obs);

  std::sort(ridx.begin(), ridx.end());

  // convert to (i, j)
  std::vector<int> cidx(ridx.size());
  for (std::size_t i = 0; i < ridx.size(); i++) {
    auto i1d = ridx[i];
    ridx[i] = i1d % n_obs;
    cidx[i] = static_cast<int>(i1d / n_obs);
  }

  return writable::list({ "i"_nm = ridx, "j"_nm = cidx });
}

struct WeightedEntry {
  int col;
  int row;
  double value;
  std::size_t seq;
};

struct CompactEntry {
  int row;
  double value;
};

std::size_t checked_triplet_count(std::size_t n_obs,
                                  std::size_t n_nbrs,
                                  const char* name) {
  if (n_nbrs == 0) {
    stop("%s must be positive", name);
  }
  const auto max_size = std::numeric_limits<std::size_t>::max();
  if (n_obs > max_size / n_nbrs || n_obs * n_nbrs > max_size / n_nbrs) {
    stop("Too many LTSA triplets to assemble");
  }
  return n_obs * n_nbrs * n_nbrs;
}

int checked_neighbor_index(int idx, std::size_t n_obs) {
  if (idx < 1 || static_cast<std::size_t>(idx) > n_obs) {
    stop("Neighborhood index is outside the sparse matrix dimensions");
  }
  return idx - 1;
}

void checked_append_output(int row,
                           double value,
                           std::vector<int>& out_i,
                           std::vector<double>& out_x,
                           std::size_t max_int) {
  if (out_i.size() >= max_int) {
    stop("Too many non-zero slots for a dgCMatrix");
  }
  out_i.push_back(row);
  out_x.push_back(value);
}

list compact_ltsa_triplet_assembly_components(const integers& value_nnt,
                                               const doubles& weights,
                                               std::size_t value_n_nbrs,
                                               std::size_t n_obs,
                                               std::size_t max_int) {
  std::vector<std::size_t> col_counts(n_obs, 0);
  for (std::size_t obs = 0; obs < n_obs; obs++) {
    std::size_t offset = obs * value_n_nbrs;
    for (std::size_t local_col = 0; local_col < value_n_nbrs; local_col++) {
      int col = checked_neighbor_index(value_nnt[offset + local_col], n_obs);
      col_counts[col] += value_n_nbrs;
    }
  }

  std::vector<std::vector<CompactEntry>> columns(n_obs);
  for (std::size_t col = 0; col < n_obs; col++) {
    columns[col].reserve(col_counts[col]);
  }

  std::size_t value_k2 = value_n_nbrs * value_n_nbrs;
  for (std::size_t obs = 0; obs < n_obs; obs++) {
    std::size_t nnt_offset = obs * value_n_nbrs;
    std::size_t weight_offset = obs * value_k2;
    for (std::size_t local_col = 0; local_col < value_n_nbrs; local_col++) {
      int col = checked_neighbor_index(value_nnt[nnt_offset + local_col], n_obs);
      for (std::size_t local_row = 0; local_row < value_n_nbrs; local_row++) {
        int row = checked_neighbor_index(value_nnt[nnt_offset + local_row], n_obs);
        double value = weights[weight_offset + local_col * value_n_nbrs + local_row];
        columns[col].push_back(CompactEntry{row, value});
      }
    }
  }

  std::vector<int> out_i;
  std::vector<double> out_x;
  std::vector<int> p(n_obs + 1, 0);

  std::vector<double> row_sums(n_obs, 0.0);
  std::vector<int> row_seen(n_obs, -1);
  std::vector<int> touched_rows;

  for (std::size_t col = 0; col < n_obs; col++) {
    auto& entries = columns[col];
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
        checked_append_output(row, value, out_i, out_x, max_int);
      }
    }
    p[col + 1] = static_cast<int>(out_i.size());
  }

  return writable::list({
    "i"_nm = out_i,
    "p"_nm = p,
    "x"_nm = out_x
  });
}

[[cpp11::register]] list
ltsa_triplet_assembly_components(const integers& pattern_nnt,
                                  std::size_t pattern_n_nbrs,
                                  const integers& value_nnt,
                                  const doubles& weights,
                                  std::size_t value_n_nbrs,
                                  bool preserve_pattern) {
  if (value_nnt.size() == 0 || value_n_nbrs == 0) {
    stop("Value neighborhoods must not be empty");
  }
  if (value_nnt.size() % value_n_nbrs != 0) {
    stop("Inconsistent value neighborhood dimensions");
  }

  std::size_t n_obs = value_nnt.size() / value_n_nbrs;

  if (preserve_pattern) {
    if (pattern_nnt.size() == 0 || pattern_n_nbrs == 0) {
      stop("Pattern neighborhoods must not be empty");
    }
    if (pattern_nnt.size() % pattern_n_nbrs != 0) {
      stop("Inconsistent pattern neighborhood dimensions");
    }
    if (pattern_nnt.size() / pattern_n_nbrs != n_obs) {
      stop("Inconsistent pattern neighborhood dimensions");
    }
  }

  if (value_nnt.size() != n_obs * value_n_nbrs) {
    stop("Inconsistent value neighborhood dimensions");
  }

  std::size_t n_pattern_triplets = preserve_pattern
    ? checked_triplet_count(n_obs, pattern_n_nbrs, "pattern_n_nbrs")
    : 0;
  std::size_t n_value_triplets =
    checked_triplet_count(n_obs, value_n_nbrs, "value_n_nbrs");
  if (weights.size() != n_value_triplets) {
    stop("Inconsistent local weight dimensions");
  }

  const auto max_int = static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (n_obs + 1 > max_int) {
    stop("Too many observations for a dgCMatrix");
  }

  if (!preserve_pattern) {
    return compact_ltsa_triplet_assembly_components(
      value_nnt,
      weights,
      value_n_nbrs,
      n_obs,
      max_int
    );
  }

  std::vector<WeightedEntry> entries;
  entries.reserve(n_pattern_triplets + n_value_triplets);

  std::size_t seq = 0;
  if (preserve_pattern) {
    for (std::size_t obs = 0; obs < n_obs; obs++) {
      std::size_t offset = obs * pattern_n_nbrs;
      for (std::size_t local_col = 0; local_col < pattern_n_nbrs; local_col++) {
        int col = checked_neighbor_index(pattern_nnt[offset + local_col], n_obs);
        for (std::size_t local_row = 0; local_row < pattern_n_nbrs; local_row++) {
          int row = checked_neighbor_index(pattern_nnt[offset + local_row], n_obs);
          entries.push_back(WeightedEntry{col, row, 0.0, seq++});
        }
      }
    }
  }

  std::size_t value_k2 = value_n_nbrs * value_n_nbrs;
  for (std::size_t obs = 0; obs < n_obs; obs++) {
    std::size_t nnt_offset = obs * value_n_nbrs;
    std::size_t weight_offset = obs * value_k2;
    for (std::size_t local_col = 0; local_col < value_n_nbrs; local_col++) {
      int col = checked_neighbor_index(value_nnt[nnt_offset + local_col], n_obs);
      for (std::size_t local_row = 0; local_row < value_n_nbrs; local_row++) {
        int row = checked_neighbor_index(value_nnt[nnt_offset + local_row], n_obs);
        double value = weights[weight_offset + local_col * value_n_nbrs + local_row];
        entries.push_back(WeightedEntry{col, row, value, seq++});
      }
    }
  }

  std::sort(entries.begin(), entries.end(), [](const WeightedEntry& lhs,
                                               const WeightedEntry& rhs) {
    if (lhs.col != rhs.col) {
      return lhs.col < rhs.col;
    }
    if (lhs.row != rhs.row) {
      return lhs.row < rhs.row;
    }
    return lhs.seq < rhs.seq;
  });

  std::vector<int> out_i;
  std::vector<double> out_x;
  out_i.reserve(entries.size());
  out_x.reserve(entries.size());
  std::vector<int> p(n_obs + 1, 0);

  std::size_t pos = 0;
  int current_col = 0;
  while (pos < entries.size()) {
    int col = entries[pos].col;
    int row = entries[pos].row;
    double value = 0.0;

    while (current_col < col) {
      p[++current_col] = static_cast<int>(out_i.size());
    }

    while (pos < entries.size() &&
           entries[pos].col == col &&
           entries[pos].row == row) {
      value += entries[pos].value;
      pos++;
    }

    if (preserve_pattern || value != 0.0) {
      checked_append_output(row, value, out_i, out_x, max_int);
    }
  }

  while (static_cast<std::size_t>(current_col) < n_obs) {
    p[++current_col] = static_cast<int>(out_i.size());
  }

  return writable::list({
    "i"_nm = out_i,
    "p"_nm = p,
    "x"_nm = out_x
  });
}
