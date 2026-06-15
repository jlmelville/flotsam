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
#include <vector>

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
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
        for (std::size_t local_col = 0; local_col < value_n_nbrs_;
             local_col++) {
          full_col_counts[nni[local_col]] += value_n_nbrs_;
        }
      } else {
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

private:
  std::size_t n_obs_;
  std::size_t value_n_nbrs_;
  std::size_t max_int_;
  std::size_t n_appended_ = 0;
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
  std::unique_ptr<GramLocalWeightsWorkspace> gram_workspace;
  if (use_gram_workspace) {
    use_row_major_gram = try_make_row_major_copy(
        x_data, n_obs, static_cast<std::size_t>(x.ncol()), row_major_x);
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
  return writable::list({"i"_nm = components.i, "p"_nm = components.p,
                         "x"_nm = components.x,
                         "rank_deficient_count"_nm = rank_deficient_count,
                         "min_local_rank"_nm = min_local_rank});
}
