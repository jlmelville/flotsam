#include <algorithm>
#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/external_pointer.hpp>
#include <cpp11/integers.hpp>
#include <limits>
#include <vector>
using namespace cpp11;

// calculates M * d where M is a sparse matrix with pointer ps and values
// xs, and d is the vector ds
[[cpp11::register]] doubles
spm_times_scalar(const integers& ps, const doubles& xs, const doubles& ds) {
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

struct CompactEntry {
  int row;
  double value;
};

std::size_t checked_triplet_count(std::size_t n_obs, std::size_t n_nbrs,
                                  const char* name);

int checked_neighbor_index(int idx, std::size_t n_obs);

void checked_append_output(int row, double value, std::vector<int>& out_i,
                           std::vector<double>& out_x, std::size_t max_int);

class LtsaTripletAssemblyBuilder {
public:
  LtsaTripletAssemblyBuilder(const integers& value_nnt,
                             std::size_t value_n_nbrs, std::size_t n_obs,
                             std::size_t max_int)
      : n_obs_(n_obs), value_n_nbrs_(value_n_nbrs), max_int_(max_int),
        columns_(n_obs) {
    std::vector<std::size_t> col_counts(n_obs_, 0);
    for (std::size_t obs = 0; obs < n_obs_; obs++) {
      std::size_t offset = obs * value_n_nbrs_;
      for (std::size_t local_col = 0; local_col < value_n_nbrs_; local_col++) {
        int col = checked_neighbor_index(value_nnt[offset + local_col], n_obs_);
        col_counts[col] += value_n_nbrs_;
      }
    }

    for (std::size_t col = 0; col < n_obs_; col++) {
      columns_[col].reserve(col_counts[col]);
    }
  }

  void append(const integers& nni, const doubles& weights) {
    if (finalized_) {
      stop("LTSA triplet builder has already been finalized");
    }
    if (n_appended_ >= n_obs_) {
      stop("Too many LTSA neighborhoods appended");
    }
    if (static_cast<std::size_t>(nni.size()) != value_n_nbrs_) {
      stop("Inconsistent value neighborhood dimensions");
    }

    std::size_t value_k2 =
        checked_triplet_count(1, value_n_nbrs_, "value_n_nbrs");
    if (static_cast<std::size_t>(weights.size()) != value_k2) {
      stop("Inconsistent local weight dimensions");
    }

    for (std::size_t local_col = 0; local_col < value_n_nbrs_; local_col++) {
      int col = checked_neighbor_index(nni[local_col], n_obs_);
      for (std::size_t local_row = 0; local_row < value_n_nbrs_; local_row++) {
        int row = checked_neighbor_index(nni[local_row], n_obs_);
        double value = weights[local_col * value_n_nbrs_ + local_row];
        columns_[col].push_back(CompactEntry{row, value});
      }
    }

    n_appended_++;
  }

  list finalize() {
    if (finalized_) {
      stop("LTSA triplet builder has already been finalized");
    }
    if (n_appended_ != n_obs_) {
      stop("Not all LTSA neighborhoods were appended");
    }

    std::vector<int> out_i;
    std::vector<double> out_x;
    std::vector<int> p(n_obs_ + 1, 0);

    std::vector<double> row_sums(n_obs_, 0.0);
    std::vector<int> row_seen(n_obs_, -1);
    std::vector<int> touched_rows;

    for (std::size_t col = 0; col < n_obs_; col++) {
      auto& entries = columns_[col];
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
          checked_append_output(row, value, out_i, out_x, max_int_);
        }
      }
      p[col + 1] = static_cast<int>(out_i.size());
    }

    finalized_ = true;
    columns_.clear();
    columns_.shrink_to_fit();

    return writable::list({"i"_nm = out_i, "p"_nm = p, "x"_nm = out_x});
  }

private:
  std::size_t n_obs_;
  std::size_t value_n_nbrs_;
  std::size_t max_int_;
  std::size_t n_appended_ = 0;
  bool finalized_ = false;
  std::vector<std::vector<CompactEntry>> columns_;
};

std::size_t checked_triplet_count(std::size_t n_obs, std::size_t n_nbrs,
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

void checked_append_output(int row, double value, std::vector<int>& out_i,
                           std::vector<double>& out_x, std::size_t max_int) {
  if (out_i.size() >= max_int) {
    stop("Too many non-zero slots for a dgCMatrix");
  }
  out_i.push_back(row);
  out_x.push_back(value);
}

using LtsaTripletAssemblyBuilderPtr =
    external_pointer<LtsaTripletAssemblyBuilder>;

LtsaTripletAssemblyBuilder* checked_ltsa_triplet_builder(SEXP builder_xptr) {
  LtsaTripletAssemblyBuilderPtr builder(builder_xptr);
  LtsaTripletAssemblyBuilder* ptr = builder.get();
  if (ptr == nullptr) {
    stop("Invalid LTSA triplet builder");
  }
  return ptr;
}

[[cpp11::register]] sexp ltsa_triplet_builder_create(const integers& value_nnt,
                                                     std::size_t value_n_nbrs) {
  if (value_nnt.size() == 0 || value_n_nbrs == 0) {
    stop("Value neighborhoods must not be empty");
  }
  if (value_nnt.size() % value_n_nbrs != 0) {
    stop("Inconsistent value neighborhood dimensions");
  }

  std::size_t n_obs = value_nnt.size() / value_n_nbrs;
  if (value_nnt.size() != n_obs * value_n_nbrs) {
    stop("Inconsistent value neighborhood dimensions");
  }

  checked_triplet_count(n_obs, value_n_nbrs, "value_n_nbrs");

  const auto max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (n_obs + 1 > max_int) {
    stop("Too many observations for a dgCMatrix");
  }

  LtsaTripletAssemblyBuilderPtr builder(
      new LtsaTripletAssemblyBuilder(value_nnt, value_n_nbrs, n_obs, max_int));
  return sexp(static_cast<SEXP>(builder));
}

[[cpp11::register]] void ltsa_triplet_builder_append(SEXP builder_xptr,
                                                     const integers& nni,
                                                     const doubles& weights) {
  checked_ltsa_triplet_builder(builder_xptr)->append(nni, weights);
}

[[cpp11::register]] list ltsa_triplet_builder_finalize(SEXP builder_xptr) {
  return checked_ltsa_triplet_builder(builder_xptr)->finalize();
}
