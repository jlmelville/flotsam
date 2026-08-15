#include "ltsa_internal.h"

int query_dsyev_workspace(int n, std::vector<double> &gram,
                          std::vector<double> &values) {
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

int query_dgesdd_workspace(int n_nbrs, int n_features, int min_dim,
                           std::vector<double> &a, std::vector<double> &d,
                           std::vector<double> &u, std::vector<double> &vt,
                           std::vector<int> &iwork) {
  char jobz = 'S';
  int m = n_nbrs;
  int n = n_features;
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
                                                     std::size_t n_features,
                                                     int ndim,
                                                     bool use_row_major)
    : n_nbrs_size(n_nbrs), n_features_size(n_features),
      n_nbrs(checked_lapack_dim(n_nbrs, "n_neighbors")),
      n_features(checked_lapack_dim(n_features, "ncol(X)")),
      requested_basis_size(
          std::min(ndim, std::min(this->n_nbrs, this->n_features))),
      neighbor_indices(
          checked_vector_size<int>(n_nbrs, "LTSA neighborhood indices")),
      centered(checked_vector_size_mul<double>(
          n_nbrs, n_features, "LTSA centered neighborhood workspace")),
      gram(checked_vector_size_mul<double>(n_nbrs, n_nbrs,
                                           "LTSA Gram workspace")),
      values(checked_vector_size<double>(n_nbrs, "LTSA Gram eigenvalues")),
      weights(checked_vector_size_mul<double>(n_nbrs, n_nbrs,
                                              "LTSA local weights")) {
  if (use_row_major) {
    row_buffer.resize(checked_vector_size_mul<double>(
        n_nbrs, n_features, "LTSA row-major neighborhood workspace"));
    col_means.resize(
        checked_vector_size<double>(n_features, "LTSA column means"));
  }
  basis_columns.reserve(
      checked_vector_size<int>(static_cast<std::size_t>(requested_basis_size),
                               "LTSA retained local basis"));
  work.resize(checked_vector_size<double>(
      static_cast<std::size_t>(
          query_dsyev_workspace(this->n_nbrs, gram, values)),
      "LTSA dsyev workspace"));
}

std::vector<int>
flat_neighbors_zero_based(const cpp11::integers &transposed_neighbor_indices,
                          std::size_t offset, std::size_t n_nbrs) {
  std::vector<int> out(
      checked_vector_size<int>(n_nbrs, "LTSA neighborhood indices"));
  fill_flat_neighbors_zero_based(transposed_neighbor_indices, offset, n_nbrs,
                                 out);
  return out;
}

void fill_flat_neighbors_zero_based(
    const cpp11::integers &transposed_neighbor_indices, std::size_t offset,
    std::size_t n_nbrs, std::vector<int> &out) {
  out.resize(checked_vector_size<int>(n_nbrs, "LTSA neighborhood indices"));
  for (std::size_t local = 0; local < n_nbrs; local++) {
    out[local] = transposed_neighbor_indices[offset + local] - 1;
  }
}

namespace {

void fill_centered_neighborhood(const cpp11::doubles_matrix<> &x,
                                const std::vector<int> &neighbor_indices,
                                std::vector<double> &centered) {
  const std::size_t n_nbrs = neighbor_indices.size();
  const std::size_t n_features = x.ncol();
  const std::size_t n_values = checked_vector_size_mul<double>(
      n_nbrs, n_features, "LTSA centered neighborhood");
  if (centered.size() != n_values) {
    centered.resize(n_values);
  }

  for (std::size_t col = 0; col < n_features; col++) {
    double mean = 0.0;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      mean += x(neighbor_indices[row], col);
    }
    mean /= static_cast<double>(n_nbrs);

    for (std::size_t row = 0; row < n_nbrs; row++) {
      centered[col * n_nbrs + row] = x(neighbor_indices[row], col) - mean;
    }
  }
}

} // namespace

void fill_centered_neighborhood_column_major(
    const double *x_data, std::size_t n_obs,
    const std::vector<int> &neighbor_indices, std::vector<double> &centered,
    std::size_t n_features) {
  const std::size_t n_nbrs = neighbor_indices.size();

  for (std::size_t col = 0; col < n_features; col++) {
    const double *col_ptr = x_data + col * n_obs;
    double mean = 0.0;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      mean += col_ptr[neighbor_indices[row]];
    }
    mean /= static_cast<double>(n_nbrs);

    double *centered_col = centered.data() + col * n_nbrs;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      centered_col[row] = col_ptr[neighbor_indices[row]] - mean;
    }
  }
}

bool row_major_copy_within_limit(std::size_t n_obs, std::size_t n_features,
                                 std::size_t max_bytes) {
  const std::size_t max_values =
      std::numeric_limits<std::size_t>::max() / sizeof(double);
  if (n_obs != 0 && n_features > max_values / n_obs) {
    return false;
  }

  const std::size_t n_values = n_obs * n_features;
  return n_values * sizeof(double) <= max_bytes;
}

void make_row_major_copy(const double *x_data, std::size_t n_obs,
                         std::size_t n_features,
                         std::vector<double> &row_major) {
  row_major.resize(checked_vector_size_mul<double>(
      n_obs, n_features, "LTSA row-major input copy"));
  for (std::size_t col = 0; col < n_features; col++) {
    const double *col_ptr = x_data + col * n_obs;
    for (std::size_t row = 0; row < n_obs; row++) {
      row_major[row * n_features + col] = col_ptr[row];
    }
  }
}

void fill_centered_neighborhood_row_major(
    const std::vector<double> &row_major,
    const std::vector<int> &neighbor_indices, std::vector<double> &row_buffer,
    std::vector<double> &col_means, std::vector<double> &centered,
    std::size_t n_features) {
  const std::size_t n_nbrs = neighbor_indices.size();

  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double *src =
        row_major.data() +
        static_cast<std::size_t>(neighbor_indices[row]) * n_features;
    double *dst = row_buffer.data() + row * n_features;
    std::copy(src, src + n_features, dst);
  }

  std::fill(col_means.begin(), col_means.end(), 0.0);
  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double *src = row_buffer.data() + row * n_features;
    for (std::size_t col = 0; col < n_features; col++) {
      col_means[col] += src[col];
    }
  }
  for (std::size_t col = 0; col < n_features; col++) {
    col_means[col] /= static_cast<double>(n_nbrs);
  }

  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double *src = row_buffer.data() + row * n_features;
    for (std::size_t col = 0; col < n_features; col++) {
      centered[col * n_nbrs + row] = src[col] - col_means[col];
    }
  }
}

void fill_weights_from_basis(std::size_t n_nbrs,
                             const std::vector<int> &basis_columns,
                             const std::vector<double> &basis,
                             std::vector<double> &weights) {
  const std::size_t n_weights =
      checked_vector_size_mul<double>(n_nbrs, n_nbrs, "LTSA local weights");
  if (weights.size() != n_weights) {
    weights.resize(n_weights);
  }
  const double constant = 1.0 / static_cast<double>(n_nbrs);

  for (std::size_t col = 0; col < n_nbrs; col++) {
    for (std::size_t row = 0; row < n_nbrs; row++) {
      double projection = constant;
      for (const int basis_col : basis_columns) {
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

int select_local_basis_columns(const std::vector<double> &values, int n_values,
                               int n_nbrs, int n_features,
                               int requested_basis_size, bool values_ascending,
                               std::vector<int> &basis_columns) {
  double max_value = 0.0;
  for (int i = 0; i < n_values; i++) {
    max_value = std::max(max_value, values[i]);
  }

  const double tol =
      max_value <= 0.0 ? 0.0
                       : static_cast<double>(std::max(n_nbrs, n_features)) *
                             max_value * std::numeric_limits<double>::epsilon();

  int rank = 0;
  if (max_value > 0.0) {
    for (int i = 0; i < n_values; i++) {
      rank += values[i] > tol;
    }
  }

  basis_columns.clear();
  const int n_keep_candidates = std::min(requested_basis_size, n_values);
  for (int col = 0; col < n_keep_candidates; col++) {
    const int basis_col = values_ascending ? n_values - 1 - col : col;
    if (values[basis_col] > tol) {
      basis_columns.push_back(basis_col);
    }
  }

  return rank;
}

namespace {

LocalWeights compute_local_weights_svd(const cpp11::doubles_matrix<> &x,
                                       const std::vector<int> &neighbor_indices,
                                       int ndim) {
  const std::size_t n_nbrs_size = neighbor_indices.size();
  const std::size_t n_features_size = x.ncol();
  const int n_nbrs = checked_lapack_dim(n_nbrs_size, "n_neighbors");
  const int n_features = checked_lapack_dim(n_features_size, "ncol(X)");
  const int min_dim = std::min(n_nbrs, n_features);
  const int max_rank = min_dim;
  const int requested_basis_size = std::min(ndim, max_rank);

  std::vector<double> centered;
  fill_centered_neighborhood(x, neighbor_indices, centered);
  std::vector<double> a = centered;
  std::vector<double> d(checked_vector_size<double>(
      static_cast<std::size_t>(min_dim), "LTSA singular values"));
  std::vector<double> u(checked_vector_size_mul<double>(
      static_cast<std::size_t>(n_nbrs), static_cast<std::size_t>(min_dim),
      "LTSA left singular vectors"));
  std::vector<double> vt(checked_vector_size_mul<double>(
      static_cast<std::size_t>(min_dim), static_cast<std::size_t>(n_features),
      "LTSA right singular vectors"));
  std::vector<int> iwork(checked_vector_size_mul<int>(
      8, static_cast<std::size_t>(min_dim), "LTSA dgesdd integer workspace"));

  char jobz = 'S';
  int m = n_nbrs;
  int n = n_features;
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
  std::vector<double> work(checked_vector_size<double>(
      static_cast<std::size_t>(lwork), "LTSA dgesdd workspace"));
  F77_CALL(dgesdd)(&jobz, &m, &n, a.data(), &lda, d.data(), u.data(), &ldu,
                   vt.data(), &ldvt, work.data(), &lwork, iwork.data(),
                   &info FCONE);
  if (info != 0) {
    cpp11::stop("LAPACK dgesdd failed with info = %d", info);
  }

  LocalWeights out;
  std::vector<int> basis_columns;
  basis_columns.reserve(
      checked_vector_size<int>(static_cast<std::size_t>(requested_basis_size),
                               "LTSA retained local basis"));
  out.rank =
      select_local_basis_columns(d, min_dim, n_nbrs, n_features,
                                 requested_basis_size, false, basis_columns);

  fill_weights_from_basis(n_nbrs_size, basis_columns, u, out.weights);
  return out;
}

LocalWeights
compute_local_weights_gram(const cpp11::doubles_matrix<> &x,
                           const std::vector<int> &neighbor_indices, int ndim) {
  const std::size_t n_nbrs_size = neighbor_indices.size();
  const std::size_t n_features_size = x.ncol();
  const int n_nbrs = checked_lapack_dim(n_nbrs_size, "n_neighbors");
  const int n_features = checked_lapack_dim(n_features_size, "ncol(X)");
  const int max_rank = std::min(n_nbrs, n_features);
  const int requested_basis_size = std::min(ndim, max_rank);

  std::vector<double> centered;
  fill_centered_neighborhood(x, neighbor_indices, centered);
  std::vector<double> gram(checked_vector_size_mul<double>(
                               n_nbrs_size, n_nbrs_size, "LTSA Gram workspace"),
                           0.0);

  char uplo = 'U';
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int n = n_nbrs;
  int k = n_features;
  int lda = n_nbrs;
  int ldc = n_nbrs;
  F77_CALL(dsyrk)(&uplo, &trans, &n, &k, &alpha, centered.data(), &lda, &beta,
                  gram.data(), &ldc FCONE FCONE);

  std::vector<double> values(checked_vector_size<double>(
      static_cast<std::size_t>(n_nbrs), "LTSA Gram eigenvalues"));
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
  std::vector<double> work(checked_vector_size<double>(
      static_cast<std::size_t>(lwork), "LTSA dsyev workspace"));
  F77_CALL(dsyev)(&jobz, &uplo, &n, gram.data(), &n, values.data(), work.data(),
                  &lwork, &info FCONE FCONE);
  if (info != 0) {
    cpp11::stop("LAPACK dsyev failed with info = %d", info);
  }

  LocalWeights out;
  std::vector<int> basis_columns;
  basis_columns.reserve(
      checked_vector_size<int>(static_cast<std::size_t>(requested_basis_size),
                               "LTSA retained local basis"));
  out.rank =
      select_local_basis_columns(values, n_nbrs, n_nbrs, n_features,
                                 requested_basis_size, true, basis_columns);

  fill_weights_from_basis(n_nbrs_size, basis_columns, gram, out.weights);
  return out;
}

} // namespace

int compute_local_weights_gram_workspace(const double *x_data,
                                         std::size_t n_obs,
                                         GramLocalWeightsWorkspace &workspace,
                                         const std::vector<double> *row_major) {
  if (row_major != nullptr) {
    fill_centered_neighborhood_row_major(
        *row_major, workspace.neighbor_indices, workspace.row_buffer,
        workspace.col_means, workspace.centered, workspace.n_features_size);
  } else {
    fill_centered_neighborhood_column_major(
        x_data, n_obs, workspace.neighbor_indices, workspace.centered,
        workspace.n_features_size);
  }

  char uplo = 'U';
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int n = workspace.n_nbrs;
  int k = workspace.n_features;
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

  const int rank = select_local_basis_columns(
      workspace.values, workspace.n_nbrs, workspace.n_nbrs,
      workspace.n_features, workspace.requested_basis_size, true,
      workspace.basis_columns);

  fill_weights_from_basis(workspace.n_nbrs_size, workspace.basis_columns,
                          workspace.gram, workspace.weights);
  return rank;
}

LocalWeights
compute_local_weights_by_shape(const cpp11::doubles_matrix<> &x,
                               const std::vector<int> &neighbor_indices,
                               int ndim) {
  if (x.ncol() == 0) {
    cpp11::stop("X must contain at least one column");
  }
  if (static_cast<std::size_t>(x.ncol()) <= neighbor_indices.size()) {
    return compute_local_weights_svd(x, neighbor_indices, ndim);
  }
  return compute_local_weights_gram(x, neighbor_indices, ndim);
}
