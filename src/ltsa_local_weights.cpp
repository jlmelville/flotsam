#include "ltsa_internal.h"

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
    cpp11::stop("LAPACK dsyev workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    cpp11::stop("LAPACK dsyev workspace is too large");
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
    cpp11::stop("LAPACK dgesdd workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    cpp11::stop("LAPACK dgesdd workspace is too large");
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

std::vector<int> flat_neighbors_zero_based(const cpp11::integers& value_nnt,
                                           std::size_t offset,
                                           std::size_t n_nbrs) {
  std::vector<int> out(n_nbrs);
  fill_flat_neighbors_zero_based(value_nnt, offset, n_nbrs, out);
  return out;
}

void fill_flat_neighbors_zero_based(const cpp11::integers& value_nnt,
                                    std::size_t offset, std::size_t n_nbrs,
                                    std::vector<int>& out) {
  out.resize(n_nbrs);
  for (std::size_t local = 0; local < n_nbrs; local++) {
    out[local] = value_nnt[offset + local] - 1;
  }
}

namespace {

void fill_centered_neighborhood(const cpp11::doubles_matrix<>& x,
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

} // namespace

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

bool row_major_copy_within_limit(std::size_t n_obs, std::size_t n_dim,
                                 std::size_t max_bytes) {
  const std::size_t max_values =
      std::numeric_limits<std::size_t>::max() / sizeof(double);
  if (n_obs != 0 && n_dim > max_values / n_obs) {
    return false;
  }

  const std::size_t n_values = n_obs * n_dim;
  return n_values * sizeof(double) <= max_bytes;
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

void fill_centered_neighborhood_row_major(const std::vector<double>& row_major,
                                          const std::vector<int>& nni,
                                          std::vector<double>& row_buffer,
                                          std::vector<double>& col_means,
                                          std::vector<double>& centered,
                                          std::size_t n_dim) {
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

namespace {

LocalWeights compute_local_weights_svd(const cpp11::doubles_matrix<>& x,
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
    cpp11::stop("LAPACK dgesdd workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    cpp11::stop("LAPACK dgesdd workspace is too large");
  }

  lwork = std::max(1, static_cast<int>(work_query));
  std::vector<double> work(lwork);
  F77_CALL(dgesdd)(&jobz, &m, &n, a.data(), &lda, d.data(), u.data(), &ldu,
                   vt.data(), &ldvt, work.data(), &lwork, iwork.data(),
                   &info FCONE);
  if (info != 0) {
    cpp11::stop("LAPACK dgesdd failed with info = %d", info);
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

LocalWeights compute_local_weights_gram(const cpp11::doubles_matrix<>& x,
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
    cpp11::stop("LAPACK dsyev workspace query failed with info = %d", info);
  }
  if (work_query > std::numeric_limits<int>::max()) {
    cpp11::stop("LAPACK dsyev workspace is too large");
  }

  lwork = std::max(1, static_cast<int>(work_query));
  std::vector<double> work(lwork);
  F77_CALL(dsyev)(&jobz, &uplo, &n, gram.data(), &n, values.data(), work.data(),
                  &lwork, &info FCONE FCONE);
  if (info != 0) {
    cpp11::stop("LAPACK dsyev failed with info = %d", info);
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

} // namespace

int compute_local_weights_gram_workspace(const double* x_data,
                                         std::size_t n_obs,
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
    cpp11::stop("LAPACK dsyev failed with info = %d", info);
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

LocalWeights compute_local_weights_shape_routed(const cpp11::doubles_matrix<>& x,
                                                const std::vector<int>& nni,
                                                int ndim) {
  if (x.ncol() == 0) {
    cpp11::stop("X must contain at least one column");
  }
  if (static_cast<std::size_t>(x.ncol()) <= nni.size()) {
    return compute_local_weights_svd(x, nni, ndim);
  }
  return compute_local_weights_gram(x, nni, ndim);
}
