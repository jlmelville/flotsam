#include <algorithm>
#include <cmath>
#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/external_pointer.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/matrix.hpp>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "pforr.h"
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

struct SparseComponents {
  std::vector<int> i;
  std::vector<int> p;
  std::vector<double> x;
};

struct LocalWeights {
  std::vector<double> weights;
  int rank = 0;
};

struct GramLocalWeightsWorkspace {
  GramLocalWeightsWorkspace(std::size_t n_nbrs, std::size_t n_dim, int ndim,
                            bool use_row_major);

  std::size_t n_nbrs_size;
  std::size_t n_dim_size;
  int n_nbrs;
  int n_dim;
  int requested;
  std::vector<int> nni;
  std::vector<double> centered;
  std::vector<double> row_buffer;
  std::vector<double> col_means;
  std::vector<double> gram;
  std::vector<double> values;
  std::vector<double> work;
  std::vector<double> weights;
  std::vector<int> keep;
};

// row-major dense copy is capped at 256 MiB
const std::size_t LTSA_ROW_MAJOR_COPY_MAX_BYTES =
    static_cast<std::size_t>(256) * 1024 * 1024;

std::size_t checked_triplet_count(std::size_t n_obs, std::size_t n_nbrs,
                                  const char* name);

int checked_neighbor_index(int idx, std::size_t n_obs);

void checked_append_output(int row, double value, std::vector<int>& out_i,
                           std::vector<double>& out_x, std::size_t max_int);

void checked_ndim(int ndim);

int checked_lapack_dim(std::size_t value, const char* name);

std::vector<int> checked_neighbors(const integers& nni, std::size_t n_obs);

std::vector<int> flat_neighbors_zero_based(const integers& value_nnt,
                                           std::size_t offset,
                                           std::size_t n_nbrs);

void fill_flat_neighbors_zero_based(const integers& value_nnt,
                                    std::size_t offset, std::size_t n_nbrs,
                                    std::vector<int>& out);

void fill_centered_neighborhood(const doubles_matrix<>& x,
                                const std::vector<int>& nni,
                                std::vector<double>& centered);

void fill_centered_neighborhood_ptr(const double* x_data, std::size_t n_obs,
                                    const std::vector<int>& nni,
                                    std::vector<double>& centered,
                                    std::size_t n_dim);

bool try_make_row_major_copy(const double* x_data, std::size_t n_obs,
                             std::size_t n_dim,
                             std::vector<double>& row_major);

void fill_centered_neighborhood_row_major(
    const std::vector<double>& row_major, const std::vector<int>& nni,
    std::vector<double>& row_buffer, std::vector<double>& col_means,
    std::vector<double>& centered, std::size_t n_dim);

void fill_weights_from_basis(std::size_t n_nbrs, const std::vector<int>& keep,
                             const std::vector<double>& basis,
                             std::vector<double>& weights);

int query_dsyev_workspace(int n, std::vector<double>& gram,
                          std::vector<double>& values);

int compute_local_weights_gram_workspace(
    const double* x_data, std::size_t n_obs,
    GramLocalWeightsWorkspace& workspace,
    const std::vector<double>* row_major);

LocalWeights compute_local_weights_svd(const doubles_matrix<>& x,
                                       const std::vector<int>& nni, int ndim);

LocalWeights compute_local_weights_gram(const doubles_matrix<>& x,
                                        const std::vector<int>& nni, int ndim);

LocalWeights compute_local_weights_shape_routed(const doubles_matrix<>& x,
                                                const std::vector<int>& nni,
                                                int ndim);

class LtsaTripletAssemblyBuilder {
public:
  LtsaTripletAssemblyBuilder(const integers& value_nnt,
                             std::size_t value_n_nbrs, std::size_t n_obs,
                             std::size_t max_int)
      : n_obs_(n_obs), value_n_nbrs_(value_n_nbrs), max_int_(max_int),
        canonical_columns_(n_obs), full_columns_(n_obs),
        append_seen_(n_obs, -1) {
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
        for (std::size_t local_col = 0; local_col < value_n_nbrs_;
             local_col++) {
          full_col_counts[nni[local_col]] += value_n_nbrs_;
        }
      } else {
        raw_entries_estimate_ += value_n_nbrs_ * (value_n_nbrs_ + 1) / 2;
        for (std::size_t local_col = 0; local_col < value_n_nbrs_;
             local_col++) {
          for (std::size_t local_row = 0; local_row <= local_col;
               local_row++) {
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

    std::vector<int> checked_nni(nni.size());
    for (R_xlen_t i = 0; i < nni.size(); i++) {
      checked_nni[i] = checked_neighbor_index(nni[i], n_obs_);
    }

    std::vector<double> checked_weights(weights.size());
    for (R_xlen_t i = 0; i < weights.size(); i++) {
      checked_weights[i] = weights[i];
    }
    append_prechecked(checked_nni, checked_weights);
  }

  void append_prechecked(const std::vector<int>& nni,
                         const std::vector<double>& weights) {
    if (finalized_) {
      stop("LTSA triplet builder has already been finalized");
    }
    if (n_appended_ >= n_obs_) {
      stop("Too many LTSA neighborhoods appended");
    }
    if (nni.size() != value_n_nbrs_) {
      stop("Inconsistent value neighborhood dimensions");
    }

    std::size_t value_k2 =
        checked_triplet_count(1, value_n_nbrs_, "value_n_nbrs");
    if (weights.size() != value_k2) {
      stop("Inconsistent local weight dimensions");
    }

    if (has_duplicate_neighbors(nni)) {
      duplicate_fallback_count_++;
      append_full_prechecked(nni, weights);
    } else {
      append_triangular_prechecked(nni, weights);
    }

    n_appended_++;
  }

  SparseComponents finalize_components() {
    if (finalized_) {
      stop("LTSA triplet builder has already been finalized");
    }
    if (n_appended_ != n_obs_) {
      stop("Not all LTSA neighborhoods were appended");
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

  list finalize() {
    SparseComponents out = finalize_components();
    return writable::list({"i"_nm = out.i, "p"_nm = out.p, "x"_nm = out.x});
  }

  std::size_t raw_entries_estimate() const { return raw_entries_estimate_; }

  std::size_t duplicate_fallback_count() const {
    return duplicate_fallback_count_;
  }

private:
  std::size_t n_obs_;
  std::size_t value_n_nbrs_;
  std::size_t max_int_;
  std::size_t n_appended_ = 0;
  std::size_t raw_entries_estimate_ = 0;
  std::size_t duplicate_fallback_count_ = 0;
  bool finalized_ = false;
  std::vector<std::vector<CompactEntry>> canonical_columns_;
  std::vector<std::vector<CompactEntry>> full_columns_;
  std::vector<int> append_seen_;

  bool has_duplicate_neighbors(const std::vector<int>& nni) {
    const int marker = static_cast<int>(n_appended_ + 1);
    for (const int idx : nni) {
      if (append_seen_[idx] == marker) {
        return true;
      }
      append_seen_[idx] = marker;
    }
    return false;
  }

  void append_full_prechecked(const std::vector<int>& nni,
                              const std::vector<double>& weights) {
    for (std::size_t local_col = 0; local_col < value_n_nbrs_; local_col++) {
      const int col = nni[local_col];
      for (std::size_t local_row = 0; local_row < value_n_nbrs_; local_row++) {
        const int row = nni[local_row];
        const double value = weights[local_col * value_n_nbrs_ + local_row];
        full_columns_[col].push_back(CompactEntry{row, value});
      }
    }
  }

  void append_triangular_prechecked(const std::vector<int>& nni,
                                    const std::vector<double>& weights) {
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

  void expand_canonical_to_full(std::vector<double>& row_sums,
                                std::vector<int>& row_seen,
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
          full_columns_[row].push_back(CompactEntry{static_cast<int>(col),
                                                    value});
        }
      }
    }

    canonical_columns_.clear();
    canonical_columns_.shrink_to_fit();
  }
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

void checked_ndim(int ndim) {
  if (ndim < 1) {
    stop("ndim must be positive");
  }
}

int checked_lapack_dim(std::size_t value, const char* name) {
  const auto max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (value > max_int) {
    stop("%s is too large for LAPACK", name);
  }
  return static_cast<int>(value);
}

int query_dsyev_workspace(int n, std::vector<double>& gram,
                          std::vector<double>& values) {
  char jobz = 'V';
  char uplo = 'U';
  int lwork = -1;
  int info = 0;
  double work_query = 0.0;

  F77_CALL(dsyev)(&jobz, &uplo, &n, gram.data(), &n, values.data(), &work_query,
                  &lwork, &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    stop("LAPACK dsyev workspace is too large");
  }

  return std::max(1, static_cast<int>(work_query));
}

int query_dgesdd_workspace(int n_nbrs, int n_dim, int min_dim,
                           std::vector<double>& a, std::vector<double>& d,
                           std::vector<double>& u, std::vector<double>& vt,
                           std::vector<int>& iwork) {
  char jobz = 'S';
  int m = n_nbrs;
  int n = n_dim;
  int lda = n_nbrs;
  int ldu = n_nbrs;
  int ldvt = min_dim;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;

  F77_CALL(dgesdd)(&jobz, &m, &n, a.data(), &lda, d.data(), u.data(), &ldu,
                   vt.data(), &ldvt, &work_query, &lwork, iwork.data(),
                   &info FCONE);
  if (info != 0) {
    stop("LAPACK dgesdd workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    stop("LAPACK dgesdd workspace is too large");
  }

  return std::max(1, static_cast<int>(work_query));
}

GramLocalWeightsWorkspace::GramLocalWeightsWorkspace(std::size_t n_nbrs,
                                                     std::size_t n_dim,
                                                     int ndim,
                                                     bool use_row_major)
    : n_nbrs_size(n_nbrs), n_dim_size(n_dim),
      n_nbrs(checked_lapack_dim(n_nbrs, "n_neighbors")),
      n_dim(checked_lapack_dim(n_dim, "ncol(X)")),
      requested(std::min(ndim, std::min(this->n_nbrs, this->n_dim))),
      nni(n_nbrs), centered(n_nbrs * n_dim), gram(n_nbrs * n_nbrs),
      values(n_nbrs), weights(n_nbrs * n_nbrs) {
  if (use_row_major) {
    row_buffer.resize(n_nbrs * n_dim);
    col_means.resize(n_dim);
  }
  keep.reserve(requested);
  work.resize(query_dsyev_workspace(this->n_nbrs, gram, values));
}

std::vector<int> checked_neighbors(const integers& nni, std::size_t n_obs) {
  if (nni.size() == 0) {
    stop("Neighborhood must not be empty");
  }

  std::vector<int> out(nni.size());
  for (R_xlen_t i = 0; i < nni.size(); i++) {
    out[i] = checked_neighbor_index(nni[i], n_obs);
  }
  return out;
}

std::vector<int> flat_neighbors_zero_based(const integers& value_nnt,
                                           std::size_t offset,
                                           std::size_t n_nbrs) {
  std::vector<int> out(n_nbrs);
  fill_flat_neighbors_zero_based(value_nnt, offset, n_nbrs, out);
  return out;
}

void fill_flat_neighbors_zero_based(const integers& value_nnt,
                                    std::size_t offset, std::size_t n_nbrs,
                                    std::vector<int>& out) {
  out.resize(n_nbrs);
  for (std::size_t local = 0; local < n_nbrs; local++) {
    out[local] = value_nnt[offset + local] - 1;
  }
}

void fill_centered_neighborhood(const doubles_matrix<>& x,
                                const std::vector<int>& nni,
                                std::vector<double>& centered) {
  const std::size_t n_nbrs = nni.size();
  const std::size_t n_dim = x.ncol();
  const std::size_t n_values = n_nbrs * n_dim;
  if (centered.size() != n_values) {
    centered.resize(n_values);
  }

  for (std::size_t col = 0; col < n_dim; col++) {
    double mean = 0.0;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      mean += x(nni[row], col);
    }
    mean /= static_cast<double>(n_nbrs);

    for (std::size_t row = 0; row < n_nbrs; row++) {
      centered[col * n_nbrs + row] = x(nni[row], col) - mean;
    }
  }
}

void fill_centered_neighborhood_ptr(const double* x_data, std::size_t n_obs,
                                    const std::vector<int>& nni,
                                    std::vector<double>& centered,
                                    std::size_t n_dim) {
  const std::size_t n_nbrs = nni.size();

  for (std::size_t col = 0; col < n_dim; col++) {
    const double* col_ptr = x_data + col * n_obs;
    double mean = 0.0;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      mean += col_ptr[nni[row]];
    }
    mean /= static_cast<double>(n_nbrs);

    double* centered_col = centered.data() + col * n_nbrs;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      centered_col[row] = col_ptr[nni[row]] - mean;
    }
  }
}

bool row_major_copy_within_limit(std::size_t n_obs, std::size_t n_dim) {
  const std::size_t max_values =
      std::numeric_limits<std::size_t>::max() / sizeof(double);
  if (n_obs != 0 && n_dim > max_values / n_obs) {
    return false;
  }

  const std::size_t n_values = n_obs * n_dim;
  return n_values * sizeof(double) <= LTSA_ROW_MAJOR_COPY_MAX_BYTES;
}

void make_row_major_copy(const double* x_data, std::size_t n_obs,
                         std::size_t n_dim, std::vector<double>& row_major) {
  row_major.resize(n_obs * n_dim);
  for (std::size_t col = 0; col < n_dim; col++) {
    const double* col_ptr = x_data + col * n_obs;
    for (std::size_t row = 0; row < n_obs; row++) {
      row_major[row * n_dim + col] = col_ptr[row];
    }
  }
}

bool try_make_row_major_copy(const double* x_data, std::size_t n_obs,
                             std::size_t n_dim,
                             std::vector<double>& row_major) {
  // Keep the row-major fast path internal and conservative. The copy is a
  // clear win for tested 6k image-shaped inputs, but large dense matrices can
  // already be memory-bound during sparse staging.
  if (!row_major_copy_within_limit(n_obs, n_dim)) {
    return false;
  }

  try {
    make_row_major_copy(x_data, n_obs, n_dim, row_major);
  } catch (const std::bad_alloc&) {
    row_major.clear();
    return false;
  }
  return true;
}

void fill_centered_neighborhood_row_major(
    const std::vector<double>& row_major, const std::vector<int>& nni,
    std::vector<double>& row_buffer, std::vector<double>& col_means,
    std::vector<double>& centered, std::size_t n_dim) {
  const std::size_t n_nbrs = nni.size();

  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double* src =
        row_major.data() + static_cast<std::size_t>(nni[row]) * n_dim;
    double* dst = row_buffer.data() + row * n_dim;
    std::copy(src, src + n_dim, dst);
  }

  std::fill(col_means.begin(), col_means.end(), 0.0);
  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double* src = row_buffer.data() + row * n_dim;
    for (std::size_t col = 0; col < n_dim; col++) {
      col_means[col] += src[col];
    }
  }
  for (std::size_t col = 0; col < n_dim; col++) {
    col_means[col] /= static_cast<double>(n_nbrs);
  }

  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double* src = row_buffer.data() + row * n_dim;
    for (std::size_t col = 0; col < n_dim; col++) {
      centered[col * n_nbrs + row] = src[col] - col_means[col];
    }
  }
}

void fill_weights_from_basis(std::size_t n_nbrs, const std::vector<int>& keep,
                             const std::vector<double>& basis,
                             std::vector<double>& weights) {
  const std::size_t n_weights = n_nbrs * n_nbrs;
  if (weights.size() != n_weights) {
    weights.resize(n_weights);
  }
  const double constant = 1.0 / static_cast<double>(n_nbrs);

  for (std::size_t col = 0; col < n_nbrs; col++) {
    for (std::size_t row = 0; row < n_nbrs; row++) {
      double projection = constant;
      for (const int basis_col : keep) {
        projection +=
            basis[row + basis_col * n_nbrs] * basis[col + basis_col * n_nbrs];
      }
      weights[col * n_nbrs + row] = -projection;
    }
  }

  for (std::size_t i = 0; i < n_nbrs; i++) {
    weights[i + i * n_nbrs] += 1.0;
  }
}

LocalWeights compute_local_weights_svd(const doubles_matrix<>& x,
                                       const std::vector<int>& nni, int ndim) {
  const std::size_t n_nbrs_size = nni.size();
  const std::size_t n_dim_size = x.ncol();
  const int n_nbrs = checked_lapack_dim(n_nbrs_size, "n_neighbors");
  const int n_dim = checked_lapack_dim(n_dim_size, "ncol(X)");
  const int min_dim = std::min(n_nbrs, n_dim);
  const int max_rank = min_dim;
  const int requested = std::min(ndim, max_rank);

  std::vector<double> centered;
  fill_centered_neighborhood(x, nni, centered);
  std::vector<double> a = centered;
  std::vector<double> d(min_dim);
  std::vector<double> u(static_cast<std::size_t>(n_nbrs) * min_dim);
  std::vector<double> vt(static_cast<std::size_t>(min_dim) * n_dim);
  std::vector<int> iwork(static_cast<std::size_t>(8) * min_dim);

  char jobz = 'S';
  int m = n_nbrs;
  int n = n_dim;
  int lda = n_nbrs;
  int ldu = n_nbrs;
  int ldvt = min_dim;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;

  F77_CALL(dgesdd)(&jobz, &m, &n, a.data(), &lda, d.data(), u.data(), &ldu,
                   vt.data(), &ldvt, &work_query, &lwork, iwork.data(),
                   &info FCONE);
  if (info != 0) {
    stop("LAPACK dgesdd workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    stop("LAPACK dgesdd workspace is too large");
  }

  lwork = std::max(1, static_cast<int>(work_query));
  std::vector<double> work(lwork);
  F77_CALL(dgesdd)(&jobz, &m, &n, a.data(), &lda, d.data(), u.data(), &ldu,
                   vt.data(), &ldvt, work.data(), &lwork, iwork.data(),
                   &info FCONE);
  if (info != 0) {
    stop("LAPACK dgesdd failed with info = %d", info);
  }

  double max_d = 0.0;
  for (int i = 0; i < min_dim; i++) {
    max_d = std::max(max_d, d[i]);
  }

  const double tol = max_d == 0.0
                         ? 0.0
                         : static_cast<double>(std::max(n_nbrs, n_dim)) *
                               max_d * std::numeric_limits<double>::epsilon();

  LocalWeights out;
  if (max_d > 0.0) {
    for (int i = 0; i < min_dim; i++) {
      out.rank += d[i] > tol;
    }
  }

  std::vector<int> keep;
  keep.reserve(requested);
  for (int col = 0; col < requested; col++) {
    if (d[col] > tol) {
      keep.push_back(col);
    }
  }

  fill_weights_from_basis(n_nbrs_size, keep, u, out.weights);
  return out;
}

LocalWeights compute_local_weights_gram(const doubles_matrix<>& x,
                                        const std::vector<int>& nni, int ndim) {
  const std::size_t n_nbrs_size = nni.size();
  const std::size_t n_dim_size = x.ncol();
  const int n_nbrs = checked_lapack_dim(n_nbrs_size, "n_neighbors");
  const int n_dim = checked_lapack_dim(n_dim_size, "ncol(X)");
  const int max_rank = std::min(n_nbrs, n_dim);
  const int requested = std::min(ndim, max_rank);

  std::vector<double> centered;
  fill_centered_neighborhood(x, nni, centered);
  std::vector<double> gram(n_nbrs_size * n_nbrs_size, 0.0);

  char uplo = 'U';
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int n = n_nbrs;
  int k = n_dim;
  int lda = n_nbrs;
  int ldc = n_nbrs;
  F77_CALL(dsyrk)(&uplo, &trans, &n, &k, &alpha, centered.data(), &lda, &beta,
                  gram.data(), &ldc FCONE FCONE);

  std::vector<double> values(n_nbrs);
  char jobz = 'V';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dsyev)(&jobz, &uplo, &n, gram.data(), &n, values.data(), &work_query,
                  &lwork, &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    stop("LAPACK dsyev workspace is too large");
  }

  lwork = std::max(1, static_cast<int>(work_query));
  std::vector<double> work(lwork);
  F77_CALL(dsyev)(&jobz, &uplo, &n, gram.data(), &n, values.data(), work.data(),
                  &lwork, &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev failed with info = %d", info);
  }

  double max_eval = 0.0;
  for (int i = 0; i < n_nbrs; i++) {
    max_eval = std::max(max_eval, values[i]);
  }

  const double eval_tol =
      max_eval <= 0.0 ? 0.0
                      : static_cast<double>(std::max(n_nbrs, n_dim)) *
                            max_eval * std::numeric_limits<double>::epsilon();

  LocalWeights out;
  if (max_eval > 0.0) {
    for (int i = 0; i < n_nbrs; i++) {
      out.rank += values[i] > eval_tol;
    }
  }

  std::vector<int> keep;
  keep.reserve(requested);
  for (int col = 0; col < requested; col++) {
    const int eig_col = n_nbrs - 1 - col;
    if (values[eig_col] > eval_tol) {
      keep.push_back(eig_col);
    }
  }

  fill_weights_from_basis(n_nbrs_size, keep, gram, out.weights);
  return out;
}

int compute_local_weights_gram_workspace(
    const double* x_data, std::size_t n_obs,
    GramLocalWeightsWorkspace& workspace,
    const std::vector<double>* row_major) {
  if (row_major != nullptr) {
    fill_centered_neighborhood_row_major(
        *row_major, workspace.nni, workspace.row_buffer, workspace.col_means,
        workspace.centered, workspace.n_dim_size);
  } else {
    fill_centered_neighborhood_ptr(x_data, n_obs, workspace.nni,
                                   workspace.centered, workspace.n_dim_size);
  }

  char uplo = 'U';
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int n = workspace.n_nbrs;
  int k = workspace.n_dim;
  int lda = workspace.n_nbrs;
  int ldc = workspace.n_nbrs;
  F77_CALL(dsyrk)(&uplo, &trans, &n, &k, &alpha, workspace.centered.data(),
                  &lda, &beta, workspace.gram.data(), &ldc FCONE FCONE);

  char jobz = 'V';
  int info = 0;
  int lwork = static_cast<int>(workspace.work.size());
  F77_CALL(dsyev)(&jobz, &uplo, &n, workspace.gram.data(), &n,
                  workspace.values.data(), workspace.work.data(), &lwork,
                  &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev failed with info = %d", info);
  }

  double max_eval = 0.0;
  for (int i = 0; i < workspace.n_nbrs; i++) {
    max_eval = std::max(max_eval, workspace.values[i]);
  }

  const double eval_tol =
      max_eval <= 0.0
          ? 0.0
          : static_cast<double>(std::max(workspace.n_nbrs, workspace.n_dim)) *
                max_eval * std::numeric_limits<double>::epsilon();

  int rank = 0;
  if (max_eval > 0.0) {
    for (int i = 0; i < workspace.n_nbrs; i++) {
      rank += workspace.values[i] > eval_tol;
    }
  }

  workspace.keep.clear();
  for (int col = 0; col < workspace.requested; col++) {
    const int eig_col = workspace.n_nbrs - 1 - col;
    if (workspace.values[eig_col] > eval_tol) {
      workspace.keep.push_back(eig_col);
    }
  }

  fill_weights_from_basis(workspace.n_nbrs_size, workspace.keep, workspace.gram,
                          workspace.weights);
  return rank;
}

LocalWeights compute_local_weights_shape_routed(const doubles_matrix<>& x,
                                                const std::vector<int>& nni,
                                                int ndim) {
  if (x.ncol() == 0) {
    stop("X must contain at least one column");
  }
  if (static_cast<std::size_t>(x.ncol()) <= nni.size()) {
    return compute_local_weights_svd(x, nni, ndim);
  }
  return compute_local_weights_gram(x, nni, ndim);
}

struct ParallelLocalWeightsWorkspace {
  ParallelLocalWeightsWorkspace(std::size_t n_nbrs, std::size_t n_dim, int ndim,
                                bool use_svd, bool use_row_major)
      : n_nbrs_size(n_nbrs), n_dim_size(n_dim),
        n_nbrs(checked_lapack_dim(n_nbrs, "n_neighbors")),
        n_dim(checked_lapack_dim(n_dim, "ncol(X)")), route_svd(use_svd),
        min_dim(std::min(this->n_nbrs, this->n_dim)),
        requested(std::min(ndim, min_dim)), nni(n_nbrs),
        centered(n_nbrs * n_dim), weights(n_nbrs * n_nbrs) {
    keep.reserve(requested);

    if (route_svd) {
      svd_a.resize(static_cast<std::size_t>(n_nbrs) * n_dim);
      d.resize(min_dim);
      u.resize(static_cast<std::size_t>(n_nbrs) * min_dim);
      vt.resize(static_cast<std::size_t>(min_dim) * n_dim);
      iwork.resize(static_cast<std::size_t>(8) * min_dim);
      svd_work.resize(query_dgesdd_workspace(
          this->n_nbrs, this->n_dim, min_dim, svd_a, d, u, vt, iwork));
    } else {
      if (use_row_major) {
        row_buffer.resize(n_nbrs * n_dim);
        col_means.resize(n_dim);
      }
      gram.resize(n_nbrs * n_nbrs);
      values.resize(n_nbrs);
      gram_work.resize(query_dsyev_workspace(this->n_nbrs, gram, values));
    }
  }

  std::size_t n_nbrs_size;
  std::size_t n_dim_size;
  int n_nbrs;
  int n_dim;
  bool route_svd;
  int min_dim;
  int requested;
  std::vector<int> nni;
  std::vector<int> keep;
  std::vector<double> centered;
  std::vector<double> weights;
  std::vector<double> row_buffer;
  std::vector<double> col_means;
  std::vector<double> gram;
  std::vector<double> values;
  std::vector<double> gram_work;
  std::vector<double> svd_a;
  std::vector<double> d;
  std::vector<double> u;
  std::vector<double> vt;
  std::vector<double> svd_work;
  std::vector<int> iwork;
};

struct ParallelWorkerDiagnostics {
  int rank_deficient_count = 0;
  int min_local_rank = std::numeric_limits<int>::max();
  int duplicate_fallback_count = 0;
  int failed_step = 0;
  int failed_info = 0;
  int failed_obs = -1;
};

struct TriangularSlotPlan {
  std::vector<std::size_t> canonical_column_starts;
  std::vector<std::size_t> canonical_column_counts;
  std::vector<std::size_t> canonical_slot_offsets;
  std::vector<std::size_t> full_column_starts;
  std::vector<std::size_t> full_column_counts;
  std::vector<std::size_t> full_slot_offsets;
  std::vector<unsigned char> duplicate_flags;
  std::size_t raw_entries = 0;
  std::size_t raw_bytes = 0;
  std::size_t duplicate_fallback_count = 0;
};

struct ReduceWorkspace {
  explicit ReduceWorkspace(std::size_t n_obs)
      : row_sums(n_obs, 0.0), row_seen(n_obs, -1) {
    touched_rows.reserve(1024);
  }

  std::vector<double> row_sums;
  std::vector<int> row_seen;
  std::vector<int> touched_rows;
};

void fill_flat_neighbors_zero_based_ptr(const int* value_ptr,
                                        std::size_t offset,
                                        std::size_t n_nbrs,
                                        std::vector<int>& out) {
  for (std::size_t local = 0; local < n_nbrs; local++) {
    out[local] = value_ptr[offset + local] - 1;
  }
}

int compute_svd_weights_workspace(ParallelLocalWeightsWorkspace& workspace,
                                  int& rank) {
  std::copy(workspace.centered.begin(), workspace.centered.end(),
            workspace.svd_a.begin());

  char jobz = 'S';
  int m = workspace.n_nbrs;
  int n = workspace.n_dim;
  int lda = workspace.n_nbrs;
  int ldu = workspace.n_nbrs;
  int ldvt = workspace.min_dim;
  int lwork = static_cast<int>(workspace.svd_work.size());
  int info = 0;

  F77_CALL(dgesdd)(&jobz, &m, &n, workspace.svd_a.data(), &lda,
                   workspace.d.data(), workspace.u.data(), &ldu,
                   workspace.vt.data(), &ldvt, workspace.svd_work.data(),
                   &lwork, workspace.iwork.data(), &info FCONE);
  if (info != 0) {
    return info;
  }

  double max_d = 0.0;
  for (int i = 0; i < workspace.min_dim; i++) {
    max_d = std::max(max_d, workspace.d[i]);
  }

  const double tol =
      max_d == 0.0
          ? 0.0
          : static_cast<double>(std::max(workspace.n_nbrs, workspace.n_dim)) *
                max_d * std::numeric_limits<double>::epsilon();

  rank = 0;
  if (max_d > 0.0) {
    for (int i = 0; i < workspace.min_dim; i++) {
      rank += workspace.d[i] > tol;
    }
  }

  workspace.keep.clear();
  for (int col = 0; col < workspace.requested; col++) {
    if (workspace.d[col] > tol) {
      workspace.keep.push_back(col);
    }
  }

  fill_weights_from_basis(workspace.n_nbrs_size, workspace.keep, workspace.u,
                          workspace.weights);
  return 0;
}

int compute_gram_weights_workspace_info(
    const double* x_data, std::size_t n_obs,
    ParallelLocalWeightsWorkspace& workspace,
    const std::vector<double>* row_major, int& rank) {
  if (row_major != nullptr) {
    fill_centered_neighborhood_row_major(
        *row_major, workspace.nni, workspace.row_buffer, workspace.col_means,
        workspace.centered, workspace.n_dim_size);
  } else {
    fill_centered_neighborhood_ptr(x_data, n_obs, workspace.nni,
                                   workspace.centered, workspace.n_dim_size);
  }

  char uplo = 'U';
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int n = workspace.n_nbrs;
  int k = workspace.n_dim;
  int lda = workspace.n_nbrs;
  int ldc = workspace.n_nbrs;
  F77_CALL(dsyrk)(&uplo, &trans, &n, &k, &alpha, workspace.centered.data(),
                  &lda, &beta, workspace.gram.data(), &ldc FCONE FCONE);

  char jobz = 'V';
  int info = 0;
  int lwork = static_cast<int>(workspace.gram_work.size());
  F77_CALL(dsyev)(&jobz, &uplo, &n, workspace.gram.data(), &n,
                  workspace.values.data(), workspace.gram_work.data(), &lwork,
                  &info FCONE FCONE);
  if (info != 0) {
    return info;
  }

  double max_eval = 0.0;
  for (int i = 0; i < workspace.n_nbrs; i++) {
    max_eval = std::max(max_eval, workspace.values[i]);
  }

  const double eval_tol =
      max_eval <= 0.0
          ? 0.0
          : static_cast<double>(std::max(workspace.n_nbrs, workspace.n_dim)) *
                max_eval * std::numeric_limits<double>::epsilon();

  rank = 0;
  if (max_eval > 0.0) {
    for (int i = 0; i < workspace.n_nbrs; i++) {
      rank += workspace.values[i] > eval_tol;
    }
  }

  workspace.keep.clear();
  for (int col = 0; col < workspace.requested; col++) {
    const int eig_col = workspace.n_nbrs - 1 - col;
    if (workspace.values[eig_col] > eval_tol) {
      workspace.keep.push_back(eig_col);
    }
  }

  fill_weights_from_basis(workspace.n_nbrs_size, workspace.keep,
                          workspace.gram, workspace.weights);
  return 0;
}

int compute_parallel_local_weights(const double* x_data, std::size_t n_obs,
                                   ParallelLocalWeightsWorkspace& workspace,
                                   const std::vector<double>* row_major,
                                   int ndim,
                                   ParallelWorkerDiagnostics& diagnostics,
                                   std::size_t obs) {
  int rank = 0;
  int info = 0;
  if (workspace.route_svd) {
    fill_centered_neighborhood_ptr(x_data, n_obs, workspace.nni,
                                   workspace.centered, workspace.n_dim_size);
    info = compute_svd_weights_workspace(workspace, rank);
  } else {
    info = compute_gram_weights_workspace_info(x_data, n_obs, workspace,
                                               row_major, rank);
  }

  if (info != 0) {
    diagnostics.failed_step = workspace.route_svd ? 1 : 2;
    diagnostics.failed_info = info;
    diagnostics.failed_obs = static_cast<int>(obs + 1);
    return info;
  }

  if (rank < ndim) {
    diagnostics.rank_deficient_count++;
    diagnostics.min_local_rank = std::min(diagnostics.min_local_rank, rank);
  }
  return 0;
}

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
  if (n_nbrs >
      (std::numeric_limits<std::size_t>::max() - 1) / n_nbrs) {
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
  const std::size_t value_bytes = checked_size_mul(
      total_count, sizeof(double), "Raw LTSA value staging buffer is too large");
  return checked_size_add(row_bytes, value_bytes,
                          "Raw LTSA staging buffers are too large");
}

template <typename T>
void checked_resize_vector(std::vector<T>& out, std::size_t n,
                           const char* name) {
  if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
    stop("%s is too large", name);
  }
  try {
    out.resize(n);
  } catch (const std::bad_alloc&) {
    stop("Unable to allocate %s", name);
  } catch (const std::length_error&) {
    stop("%s is too large", name);
  }
}

TriangularSlotPlan assign_triangular_two_pass_slots_flat(
    const int* value_ptr, std::size_t n_obs, std::size_t n_nbrs) {
  checked_triplet_count(n_obs, n_nbrs, "value_n_nbrs");

  TriangularSlotPlan plan;
  const std::size_t tri_count = triangular_pair_count(n_nbrs);
  plan.canonical_column_counts.assign(n_obs, 0);
  plan.full_column_counts.assign(n_obs, 0);
  checked_resize_vector(
      plan.canonical_slot_offsets,
      checked_size_mul(n_obs, tri_count,
                       "Too many triangular LTSA slot offsets"),
      "triangular LTSA slot offsets");
  checked_resize_vector(
      plan.full_slot_offsets,
      checked_size_mul(n_obs, n_nbrs, "Too many full LTSA slot offsets"),
      "full LTSA slot offsets");
  plan.duplicate_flags.assign(n_obs, 0);

  std::vector<int> nni(n_nbrs);
  std::vector<int> seen(n_obs, -1);
  for (std::size_t obs = 0; obs < n_obs; obs++) {
    const std::size_t offset = obs * n_nbrs;
    const int marker = static_cast<int>(obs + 1);
    bool has_duplicate = false;

    for (std::size_t local = 0; local < n_nbrs; local++) {
      const int idx = checked_neighbor_index(value_ptr[offset + local], n_obs);
      nni[local] = idx;
      if (seen[idx] == marker) {
        has_duplicate = true;
      }
      seen[idx] = marker;
    }

    if (has_duplicate) {
      plan.duplicate_flags[obs] = 1;
      plan.duplicate_fallback_count++;
      for (std::size_t local_col = 0; local_col < n_nbrs; local_col++) {
        const int col = nni[local_col];
        plan.full_slot_offsets[offset + local_col] =
            plan.full_column_counts[col];
        plan.full_column_counts[col] = checked_size_add(
            plan.full_column_counts[col], n_nbrs,
            "Too many duplicate LTSA contributions to stage");
      }
      continue;
    }

    const std::size_t obs_tri_offset = obs * tri_count;
    for (std::size_t local_col = 0; local_col < n_nbrs; local_col++) {
      for (std::size_t local_row = 0; local_row <= local_col; local_row++) {
        const int global_row = nni[local_row];
        const int global_col = nni[local_col];
        const int col = std::max(global_row, global_col);
        const std::size_t pair_offset =
            obs_tri_offset + triangular_pair_offset(local_col, local_row);
        plan.canonical_slot_offsets[pair_offset] =
            plan.canonical_column_counts[col];
        plan.canonical_column_counts[col] = checked_size_add(
            plan.canonical_column_counts[col], 1,
            "Too many triangular LTSA contributions to stage");
      }
    }
  }

  plan.canonical_column_starts.assign(n_obs + 1, 0);
  plan.full_column_starts.assign(n_obs + 1, 0);
  for (std::size_t col = 0; col < n_obs; col++) {
    plan.canonical_column_starts[col + 1] = checked_size_add(
        plan.canonical_column_starts[col], plan.canonical_column_counts[col],
        "Too many triangular LTSA contributions to stage");
    plan.full_column_starts[col + 1] =
        checked_size_add(plan.full_column_starts[col],
                         plan.full_column_counts[col],
                         "Too many duplicate LTSA contributions to stage");
  }

  plan.raw_entries = checked_size_add(
      plan.canonical_column_starts[n_obs], plan.full_column_starts[n_obs],
      "Too many raw LTSA contributions to stage");
  plan.raw_bytes = checked_raw_staging_bytes(
      plan.canonical_column_starts[n_obs], plan.full_column_starts[n_obs]);
  return plan;
}

struct ParallelTriangularFillWorker {
  const double* x_data;
  const std::vector<double>* row_major;
  const int* value_ptr;
  std::size_t n_obs;
  std::size_t n_nbrs;
  std::size_t tri_count;
  int ndim;
  const std::vector<std::size_t>* canonical_column_starts;
  const std::vector<std::size_t>* canonical_slot_offsets;
  std::vector<int>* canonical_raw_rows;
  std::vector<double>* canonical_raw_values;
  const std::vector<std::size_t>* full_column_starts;
  const std::vector<std::size_t>* full_slot_offsets;
  std::vector<int>* full_raw_rows;
  std::vector<double>* full_raw_values;
  const std::vector<unsigned char>* duplicate_flags;
  std::vector<ParallelLocalWeightsWorkspace>* workspaces;
  std::vector<ParallelWorkerDiagnostics>* diagnostics;

  void operator()(std::size_t begin, std::size_t end, std::size_t chunk_id) {
    ParallelLocalWeightsWorkspace& workspace = (*workspaces)[chunk_id];
    ParallelWorkerDiagnostics& worker_diagnostics = (*diagnostics)[chunk_id];

    for (std::size_t obs = begin; obs < end; obs++) {
      const std::size_t offset = obs * n_nbrs;
      fill_flat_neighbors_zero_based_ptr(value_ptr, offset, n_nbrs,
                                         workspace.nni);
      if (compute_parallel_local_weights(x_data, n_obs, workspace, row_major,
                                         ndim, worker_diagnostics, obs) != 0) {
        break;
      }

      if ((*duplicate_flags)[obs]) {
        worker_diagnostics.duplicate_fallback_count++;
        for (std::size_t local_col = 0; local_col < n_nbrs; local_col++) {
          const int col = workspace.nni[local_col];
          const std::size_t slot =
              (*full_column_starts)[col] + (*full_slot_offsets)[offset + local_col];
          for (std::size_t local_row = 0; local_row < n_nbrs; local_row++) {
            const std::size_t pos = slot + local_row;
            (*full_raw_rows)[pos] = workspace.nni[local_row];
            (*full_raw_values)[pos] =
                workspace.weights[local_col * n_nbrs + local_row];
          }
        }
      } else {
        const std::size_t obs_tri_offset = obs * tri_count;
        for (std::size_t local_col = 0; local_col < n_nbrs; local_col++) {
          for (std::size_t local_row = 0; local_row <= local_col; local_row++) {
            const int global_row = workspace.nni[local_row];
            const int global_col = workspace.nni[local_col];
            const int row = std::min(global_row, global_col);
            const int col = std::max(global_row, global_col);
            const std::size_t pair_offset =
                obs_tri_offset + triangular_pair_offset(local_col, local_row);
            const std::size_t pos =
                (*canonical_column_starts)[col] +
                (*canonical_slot_offsets)[pair_offset];
            (*canonical_raw_rows)[pos] = row;
            (*canonical_raw_values)[pos] =
                workspace.weights[local_col * n_nbrs + local_row];
          }
        }
      }
    }
  }
};

struct ColumnReduceWorker {
  const std::vector<std::size_t>* column_starts;
  const std::vector<std::size_t>* column_counts;
  const std::vector<int>* raw_rows;
  const std::vector<double>* raw_values;
  std::vector<std::vector<CompactEntry>>* reduced_columns;
  std::vector<ReduceWorkspace>* workspaces;

  void operator()(std::size_t begin, std::size_t end, std::size_t chunk_id) {
    ReduceWorkspace& workspace = (*workspaces)[chunk_id];

    for (std::size_t col = begin; col < end; col++) {
      const int marker = static_cast<int>(col);
      workspace.touched_rows.clear();
      const std::size_t start = (*column_starts)[col];
      const std::size_t count = (*column_counts)[col];
      for (std::size_t pos = start; pos < start + count; pos++) {
        const int row = (*raw_rows)[pos];
        if (workspace.row_seen[row] != marker) {
          workspace.row_seen[row] = marker;
          workspace.row_sums[row] = 0.0;
          workspace.touched_rows.push_back(row);
        }
        workspace.row_sums[row] += (*raw_values)[pos];
      }

      std::sort(workspace.touched_rows.begin(), workspace.touched_rows.end());
      std::vector<CompactEntry>& out = (*reduced_columns)[col];
      out.reserve(workspace.touched_rows.size());
      for (const int row : workspace.touched_rows) {
        const double value = workspace.row_sums[row];
        if (value != 0.0) {
          out.push_back(CompactEntry{row, value});
        }
      }
    }
  }
};

std::vector<std::vector<CompactEntry>> reduce_raw_columns_parallel(
    const std::vector<std::size_t>& column_starts,
    const std::vector<std::size_t>& column_counts,
    const std::vector<int>& raw_rows, const std::vector<double>& raw_values,
    std::size_t n_obs, std::size_t n_threads) {
  std::vector<pforr::IndexRange> ranges =
      pforr::split_input_range(pforr::IndexRange(0, n_obs), n_threads, 1);
  std::vector<std::vector<CompactEntry>> reduced_columns(n_obs);
  std::vector<ReduceWorkspace> workspaces;
  workspaces.reserve(ranges.size());
  for (std::size_t chunk = 0; chunk < ranges.size(); chunk++) {
    workspaces.emplace_back(n_obs);
  }

  ColumnReduceWorker worker{&column_starts, &column_counts, &raw_rows,
                            &raw_values, &reduced_columns, &workspaces};
  pforr::parallel_for_indexed(0, n_obs, worker, n_threads, 1);
  return reduced_columns;
}

void expand_canonical_columns_to_full(
    const std::vector<std::vector<CompactEntry>>& canonical_columns,
    std::vector<std::vector<CompactEntry>>& full_columns) {
  for (std::size_t col = 0; col < canonical_columns.size(); col++) {
    for (const CompactEntry& entry : canonical_columns[col]) {
      full_columns[col].push_back(entry);
      if (entry.row != static_cast<int>(col)) {
        full_columns[entry.row].push_back(
            CompactEntry{static_cast<int>(col), entry.value});
      }
    }
  }
}

void append_reduced_columns(
    const std::vector<std::vector<CompactEntry>>& source,
    std::vector<std::vector<CompactEntry>>& target) {
  for (std::size_t col = 0; col < source.size(); col++) {
    target[col].insert(target[col].end(), source[col].begin(),
                       source[col].end());
  }
}

SparseComponents finalize_compact_columns(
    const std::vector<std::vector<CompactEntry>>& columns, std::size_t n_obs,
    std::size_t max_int) {
  SparseComponents out;
  out.p.resize(n_obs + 1, 0);

  std::vector<double> row_sums(n_obs, 0.0);
  std::vector<int> row_seen(n_obs, -1);
  std::vector<int> touched_rows;
  touched_rows.reserve(1024);

  for (std::size_t col = 0; col < n_obs; col++) {
    const int marker = static_cast<int>(col);
    touched_rows.clear();
    for (const CompactEntry& entry : columns[col]) {
      if (row_seen[entry.row] != marker) {
        row_seen[entry.row] = marker;
        row_sums[entry.row] = 0.0;
        touched_rows.push_back(entry.row);
      }
      row_sums[entry.row] += entry.value;
    }

    std::sort(touched_rows.begin(), touched_rows.end());
    for (const int row : touched_rows) {
      const double value = row_sums[row];
      if (value != 0.0) {
        checked_append_output(row, value, out.i, out.x, max_int);
      }
    }
    out.p[col + 1] = static_cast<int>(out.i.size());
  }

  return out;
}

SparseComponents finalize_triangular_two_pass_raw(
    const TriangularSlotPlan& plan,
    const std::vector<int>& canonical_raw_rows,
    const std::vector<double>& canonical_raw_values,
    const std::vector<int>& full_raw_rows,
    const std::vector<double>& full_raw_values, std::size_t n_obs,
    std::size_t n_threads, std::size_t max_int) {
  std::vector<std::vector<CompactEntry>> canonical_columns =
      reduce_raw_columns_parallel(plan.canonical_column_starts,
                                  plan.canonical_column_counts,
                                  canonical_raw_rows, canonical_raw_values,
                                  n_obs, n_threads);

  std::vector<std::vector<CompactEntry>> full_columns(n_obs);
  expand_canonical_columns_to_full(canonical_columns, full_columns);

  if (!full_raw_rows.empty()) {
    std::vector<std::vector<CompactEntry>> duplicate_columns =
        reduce_raw_columns_parallel(plan.full_column_starts,
                                    plan.full_column_counts, full_raw_rows,
                                    full_raw_values, n_obs, n_threads);
    append_reduced_columns(duplicate_columns, full_columns);
  }

  return finalize_compact_columns(full_columns, n_obs, max_int);
}

ParallelWorkerDiagnostics combine_worker_diagnostics(
    const std::vector<ParallelWorkerDiagnostics>& diagnostics, int ndim) {
  ParallelWorkerDiagnostics out;
  for (const ParallelWorkerDiagnostics& worker : diagnostics) {
    out.rank_deficient_count += worker.rank_deficient_count;
    out.min_local_rank = std::min(out.min_local_rank, worker.min_local_rank);
    out.duplicate_fallback_count += worker.duplicate_fallback_count;
  }
  if (out.min_local_rank == std::numeric_limits<int>::max()) {
    out.min_local_rank = ndim;
  }
  return out;
}

void stop_on_parallel_worker_failure(
    const std::vector<ParallelWorkerDiagnostics>& diagnostics) {
  for (std::size_t worker = 0; worker < diagnostics.size(); worker++) {
    const ParallelWorkerDiagnostics& current = diagnostics[worker];
    if (current.failed_step != 0) {
      const char* routine = current.failed_step == 1 ? "dgesdd" : "dsyev";
      stop("LAPACK %s failed in LTSA assembly worker %d at neighborhood %d "
           "with info = %d",
           routine, static_cast<int>(worker + 1), current.failed_obs,
           current.failed_info);
    }
  }
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
  if (static_cast<std::size_t>(value_nnt.size()) != n_obs * value_n_nbrs) {
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

[[cpp11::register]] list ltsa_local_weights_cpp(const doubles_matrix<>& x,
                                                const integers& nni, int ndim) {
  checked_ndim(ndim);
  std::vector<int> checked_nni = checked_neighbors(nni, x.nrow());
  LocalWeights local = compute_local_weights_shape_routed(x, checked_nni, ndim);

  const std::size_t n_nbrs = checked_nni.size();
  writable::doubles_matrix<> weights(n_nbrs, n_nbrs);
  for (std::size_t col = 0; col < n_nbrs; col++) {
    for (std::size_t row = 0; row < n_nbrs; row++) {
      weights(row, col) = local.weights[col * n_nbrs + row];
    }
  }

  return writable::list({"Wi"_nm = weights, "rank"_nm = local.rank});
}

[[cpp11::register]] list ltsa_assemble_local_weights(const doubles_matrix<>& x,
                                                     const integers& value_nnt,
                                                     std::size_t value_n_nbrs,
                                                     int ndim) {
  checked_ndim(ndim);
  if (value_nnt.size() == 0 || value_n_nbrs == 0) {
    stop("Value neighborhoods must not be empty");
  }
  if (value_nnt.size() % value_n_nbrs != 0) {
    stop("Inconsistent value neighborhood dimensions");
  }

  std::size_t n_obs = value_nnt.size() / value_n_nbrs;
  if (static_cast<std::size_t>(x.nrow()) != n_obs) {
    stop("Inconsistent input and neighborhood dimensions");
  }

  const auto max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (n_obs + 1 > max_int) {
    stop("Too many observations for a dgCMatrix");
  }

  LtsaTripletAssemblyBuilder builder(value_nnt, value_n_nbrs, n_obs, max_int);

  int rank_deficient_count = 0;
  int min_local_rank = ndim;
  const bool use_gram_workspace =
      x.ncol() != 0 && static_cast<std::size_t>(x.ncol()) > value_n_nbrs;
  const double* x_data = use_gram_workspace ? REAL(x.data()) : nullptr;
  std::vector<double> row_major_x;
  bool use_row_major_gram = false;
  bool row_major_within_limit = false;
  std::unique_ptr<GramLocalWeightsWorkspace> gram_workspace;
  if (use_gram_workspace) {
    row_major_within_limit =
        row_major_copy_within_limit(n_obs, static_cast<std::size_t>(x.ncol()));
    if (row_major_within_limit) {
      try {
        make_row_major_copy(x_data, n_obs, static_cast<std::size_t>(x.ncol()),
                            row_major_x);
        use_row_major_gram = true;
      } catch (const std::bad_alloc&) {
        row_major_x.clear();
      } catch (const std::length_error&) {
        row_major_x.clear();
      }
    }
    gram_workspace.reset(new GramLocalWeightsWorkspace(
        value_n_nbrs, static_cast<std::size_t>(x.ncol()), ndim,
        use_row_major_gram));
  }

  for (std::size_t obs = 0; obs < n_obs; obs++) {
    const std::size_t offset = obs * value_n_nbrs;

    if (use_gram_workspace) {
      fill_flat_neighbors_zero_based(value_nnt, offset, value_n_nbrs,
                                     gram_workspace->nni);
      int rank = compute_local_weights_gram_workspace(
          x_data, n_obs, *gram_workspace,
          use_row_major_gram ? &row_major_x : nullptr);
      if (rank < ndim) {
        rank_deficient_count++;
        min_local_rank = std::min(min_local_rank, rank);
      }
      builder.append_prechecked(gram_workspace->nni, gram_workspace->weights);
    } else {
      std::vector<int> local_nni =
          flat_neighbors_zero_based(value_nnt, offset, value_n_nbrs);
      LocalWeights local =
          compute_local_weights_shape_routed(x, local_nni, ndim);
      if (local.rank < ndim) {
        rank_deficient_count++;
        min_local_rank = std::min(min_local_rank, local.rank);
      }
      builder.append_prechecked(local_nni, local.weights);
    }
  }

  SparseComponents components = builder.finalize_components();
  const std::size_t raw_entries = builder.raw_entries_estimate();
  const std::size_t raw_bytes = checked_raw_staging_bytes(raw_entries, 0);
  const std::string fallback_reason = row_major_fallback_reason(
      use_gram_workspace, use_row_major_gram, row_major_within_limit);

  return writable::list({"i"_nm = components.i, "p"_nm = components.p,
                         "x"_nm = components.x,
                         "rank_deficient_count"_nm = rank_deficient_count,
                         "min_local_rank"_nm = min_local_rank,
                         "assembly_route"_nm = "serial_triangular",
                         "requested_assembly_threads"_nm = 1,
                         "effective_assembly_threads"_nm = 1,
                         "raw_entries_estimate"_nm =
                             static_cast<double>(raw_entries),
                         "raw_bytes_estimate"_nm =
                             static_cast<double>(raw_bytes),
                         "duplicate_fallback_count"_nm =
                             static_cast<int>(
                                 builder.duplicate_fallback_count()),
                         "row_major_used"_nm = use_row_major_gram,
                         "row_major_fallback_reason"_nm = fallback_reason,
                         "parallel_fallback_reason"_nm = "not_requested"});
}

[[cpp11::register]] list
ltsa_assemble_local_weights_parallel(const doubles_matrix<>& x,
                                     const integers& value_nnt,
                                     std::size_t value_n_nbrs, int ndim,
                                     int requested_threads) {
  checked_ndim(ndim);
  if (requested_threads < 1) {
    stop("n_assembly_threads must be positive");
  }
  if (value_nnt.size() == 0 || value_n_nbrs == 0) {
    stop("Value neighborhoods must not be empty");
  }
  if (value_nnt.size() % value_n_nbrs != 0) {
    stop("Inconsistent value neighborhood dimensions");
  }

  const std::size_t n_obs = value_nnt.size() / value_n_nbrs;
  if (static_cast<std::size_t>(x.nrow()) != n_obs) {
    stop("Inconsistent input and neighborhood dimensions");
  }
  if (x.ncol() == 0) {
    stop("X must contain at least one column");
  }

  const auto max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (n_obs + 1 > max_int) {
    stop("Too many observations for a dgCMatrix");
  }

  const int* value_ptr = INTEGER(value_nnt.data());
  TriangularSlotPlan slot_plan =
      assign_triangular_two_pass_slots_flat(value_ptr, n_obs, value_n_nbrs);

  const std::size_t requested_thread_count =
      static_cast<std::size_t>(requested_threads);
  const std::vector<pforr::IndexRange> obs_ranges =
      pforr::split_input_range(pforr::IndexRange(0, n_obs),
                               requested_thread_count, 1);
  const std::size_t effective_threads = obs_ranges.size();

  const bool use_svd_route =
      static_cast<std::size_t>(x.ncol()) <= value_n_nbrs;
  const double* x_data = REAL(x.data());
  std::vector<double> row_major_x;
  const std::vector<double>* row_major_ptr = nullptr;
  bool row_major_within_limit = false;
  if (!use_svd_route) {
    row_major_within_limit =
        row_major_copy_within_limit(n_obs, static_cast<std::size_t>(x.ncol()));
    if (row_major_within_limit) {
      try {
        make_row_major_copy(x_data, n_obs, static_cast<std::size_t>(x.ncol()),
                            row_major_x);
        row_major_ptr = &row_major_x;
      } catch (const std::bad_alloc&) {
        row_major_x.clear();
      } catch (const std::length_error&) {
        row_major_x.clear();
      }
    }
  }

  std::vector<ParallelWorkerDiagnostics> worker_diagnostics(effective_threads);
  std::vector<ParallelLocalWeightsWorkspace> workspaces;
  workspaces.reserve(effective_threads);
  for (std::size_t chunk = 0; chunk < effective_threads; chunk++) {
    workspaces.emplace_back(value_n_nbrs, static_cast<std::size_t>(x.ncol()),
                            ndim, use_svd_route, row_major_ptr != nullptr);
  }

  const std::size_t canonical_count =
      slot_plan.canonical_column_starts[n_obs];
  const std::size_t full_count = slot_plan.full_column_starts[n_obs];
  std::vector<int> canonical_raw_rows;
  std::vector<double> canonical_raw_values;
  std::vector<int> full_raw_rows;
  std::vector<double> full_raw_values;
  checked_resize_vector(canonical_raw_rows, canonical_count,
                        "canonical raw LTSA row buffer");
  checked_resize_vector(canonical_raw_values, canonical_count,
                        "canonical raw LTSA value buffer");
  checked_resize_vector(full_raw_rows, full_count,
                        "duplicate raw LTSA row buffer");
  checked_resize_vector(full_raw_values, full_count,
                        "duplicate raw LTSA value buffer");

  ParallelTriangularFillWorker worker{
      x_data,
      row_major_ptr,
      value_ptr,
      n_obs,
      value_n_nbrs,
      triangular_pair_count(value_n_nbrs),
      ndim,
      &slot_plan.canonical_column_starts,
      &slot_plan.canonical_slot_offsets,
      &canonical_raw_rows,
      &canonical_raw_values,
      &slot_plan.full_column_starts,
      &slot_plan.full_slot_offsets,
      &full_raw_rows,
      &full_raw_values,
      &slot_plan.duplicate_flags,
      &workspaces,
      &worker_diagnostics};

  pforr::parallel_for_indexed(0, n_obs, worker, requested_thread_count, 1);
  stop_on_parallel_worker_failure(worker_diagnostics);

  SparseComponents components = finalize_triangular_two_pass_raw(
      slot_plan, canonical_raw_rows, canonical_raw_values, full_raw_rows,
      full_raw_values, n_obs, requested_thread_count, max_int);
  ParallelWorkerDiagnostics diagnostics =
      combine_worker_diagnostics(worker_diagnostics, ndim);
  const std::string fallback_reason =
      row_major_fallback_reason(!use_svd_route, row_major_ptr != nullptr,
                                row_major_within_limit);

  return writable::list({"i"_nm = components.i, "p"_nm = components.p,
                         "x"_nm = components.x,
                         "rank_deficient_count"_nm =
                             diagnostics.rank_deficient_count,
                         "min_local_rank"_nm = diagnostics.min_local_rank,
                         "assembly_route"_nm =
                             "parallel_triangular_two_pass",
                         "requested_assembly_threads"_nm = requested_threads,
                         "effective_assembly_threads"_nm =
                             static_cast<int>(effective_threads),
                         "raw_entries_estimate"_nm =
                             static_cast<double>(slot_plan.raw_entries),
                         "raw_bytes_estimate"_nm =
                             static_cast<double>(slot_plan.raw_bytes),
                         "duplicate_fallback_count"_nm =
                             diagnostics.duplicate_fallback_count,
                         "row_major_used"_nm = row_major_ptr != nullptr,
                         "row_major_fallback_reason"_nm = fallback_reason,
                         "parallel_fallback_reason"_nm = ""});
}
