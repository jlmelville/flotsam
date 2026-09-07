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

SvdLocalWeightsWorkspace::SvdLocalWeightsWorkspace(std::size_t n_nbrs,
                                                   std::size_t n_features,
                                                   int ndim)
    : n_nbrs_size(n_nbrs), n_features_size(n_features),
      n_nbrs(checked_lapack_dim(n_nbrs, "n_neighbors")),
      n_features(checked_lapack_dim(n_features, "ncol(X)")),
      min_dim(std::min(this->n_nbrs, this->n_features)),
      requested_basis_size(std::min(ndim, min_dim)),
      neighbor_indices(
          checked_vector_size<int>(n_nbrs, "LTSA neighborhood indices")),
      centered(checked_vector_size_mul<double>(
          n_nbrs, n_features, "LTSA centered neighborhood workspace")),
      a(checked_vector_size_mul<double>(n_nbrs, n_features,
                                        "LTSA dgesdd matrix workspace")),
      d(checked_vector_size<double>(static_cast<std::size_t>(min_dim),
                                    "LTSA singular values")),
      u(checked_vector_size_mul<double>(n_nbrs,
                                        static_cast<std::size_t>(min_dim),
                                        "LTSA left singular vectors")),
      vt(checked_vector_size_mul<double>(static_cast<std::size_t>(min_dim),
                                         n_features,
                                         "LTSA right singular vectors")),
      iwork(checked_vector_size_mul<int>(8, static_cast<std::size_t>(min_dim),
                                         "LTSA dgesdd integer workspace")),
      weights(checked_vector_size_mul<double>(n_nbrs, n_nbrs,
                                              "LTSA local weights")) {
  basis_columns.reserve(
      checked_vector_size<int>(static_cast<std::size_t>(requested_basis_size),
                               "LTSA retained local basis"));
  work.resize(checked_vector_size<double>(
      static_cast<std::size_t>(query_dgesdd_workspace(
          this->n_nbrs, this->n_features, min_dim, a, d, u, vt, iwork)),
      "LTSA dgesdd workspace"));
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
  }
  basis_columns.reserve(
      checked_vector_size<int>(static_cast<std::size_t>(requested_basis_size),
                               "LTSA retained local basis"));
  work.resize(checked_vector_size<double>(
      static_cast<std::size_t>(
          query_dsyev_workspace(this->n_nbrs, gram, values)),
      "LTSA dsyev workspace"));
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

template <typename ValueAt>
bool fill_centered_column(std::size_t n_nbrs, ValueAt value_at,
                          double *centered_col) {
  const double anchor = value_at(0);
  if (!std::isfinite(anchor)) {
    return false;
  }

  long double sum = 0.0;
  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double value = value_at(row);
    const double delta = value - anchor;
    if (!std::isfinite(value) || !std::isfinite(delta)) {
      return false;
    }
    centered_col[row] = delta;
    sum += static_cast<long double>(delta);
    if (!std::isfinite(sum)) {
      return false;
    }
  }

  const double mean_delta =
      static_cast<double>(sum / static_cast<long double>(n_nbrs));
  if (!std::isfinite(mean_delta)) {
    return false;
  }
  for (std::size_t row = 0; row < n_nbrs; row++) {
    centered_col[row] -= mean_delta;
    if (!std::isfinite(centered_col[row])) {
      return false;
    }
  }
  return true;
}

bool scale_centered_neighborhood(std::vector<double> &centered) {
  double scale = 0.0;
  for (const double value : centered) {
    if (!std::isfinite(value)) {
      return false;
    }
    scale = std::max(scale, std::abs(value));
  }

  if (scale == 0.0) {
    return true;
  }
  for (double &value : centered) {
    value /= scale;
  }
  return true;
}

} // namespace

bool fill_centered_neighborhood_column_major(
    const double *x_data, std::size_t n_obs,
    const std::vector<int> &neighbor_indices, std::vector<double> &centered,
    std::size_t n_features) {
  const std::size_t n_nbrs = neighbor_indices.size();

  for (std::size_t col = 0; col < n_features; col++) {
    const double *col_ptr = x_data + col * n_obs;
    double *centered_col = centered.data() + col * n_nbrs;
    const auto value_at = [&](std::size_t row) {
      return col_ptr[neighbor_indices[row]];
    };
    if (!fill_centered_column(n_nbrs, value_at, centered_col)) {
      return false;
    }
  }
  return scale_centered_neighborhood(centered);
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
  // This copy is completed on the main thread before workers are launched.
  for (std::size_t col = 0; col < n_features; col++) {
    if (col % 64 == 0) {
      cpp11::check_user_interrupt();
    }
    const double *col_ptr = x_data + col * n_obs;
    for (std::size_t row = 0; row < n_obs; row++) {
      row_major[row * n_features + col] = col_ptr[row];
    }
  }
}

bool fill_centered_neighborhood_row_major(
    const std::vector<double> &row_major,
    const std::vector<int> &neighbor_indices, std::vector<double> &row_buffer,
    std::vector<double> &centered, std::size_t n_features) {
  const std::size_t n_nbrs = neighbor_indices.size();

  for (std::size_t row = 0; row < n_nbrs; row++) {
    const double *src =
        row_major.data() +
        static_cast<std::size_t>(neighbor_indices[row]) * n_features;
    double *dst = row_buffer.data() + row * n_features;
    std::copy(src, src + n_features, dst);
  }

  for (std::size_t col = 0; col < n_features; col++) {
    double *centered_col = centered.data() + col * n_nbrs;
    const auto value_at = [&](std::size_t row) {
      return row_buffer[row * n_features + col];
    };
    if (!fill_centered_column(n_nbrs, value_at, centered_col)) {
      return false;
    }
  }
  return scale_centered_neighborhood(centered);
}

int clean_local_basis(std::size_t n_nbrs, std::vector<int> &basis_columns,
                      std::vector<double> &basis) {
  std::size_t n_retained = 0;
  const std::size_t n_selected = basis_columns.size();
  const double drop_tolerance =
      std::sqrt(std::numeric_limits<double>::epsilon());

  for (const int basis_col : basis_columns) {
    double *candidate =
        basis.data() + static_cast<std::size_t>(basis_col) * n_nbrs;
    double original_norm = 0.0;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      if (!std::isfinite(candidate[row])) {
        return -1;
      }
      original_norm = std::hypot(original_norm, candidate[row]);
    }

    // LAPACK supplies unit vectors. The cutoff below detects a direction that
    // numerical mean projection annihilates; it is unrelated to data scale.
    for (int pass = 0; pass < 2; pass++) {
      long double sum = 0.0;
      for (std::size_t row = 0; row < n_nbrs; row++) {
        sum += static_cast<long double>(candidate[row]);
      }
      const double mean =
          static_cast<double>(sum / static_cast<long double>(n_nbrs));
      for (std::size_t row = 0; row < n_nbrs; row++) {
        candidate[row] -= mean;
      }

      for (std::size_t retained = 0; retained < n_retained; retained++) {
        const int retained_col = basis_columns[retained];
        const double *previous =
            basis.data() + static_cast<std::size_t>(retained_col) * n_nbrs;
        long double dot = 0.0;
        for (std::size_t row = 0; row < n_nbrs; row++) {
          dot += static_cast<long double>(previous[row]) * candidate[row];
        }
        const double projection = static_cast<double>(dot);
        for (std::size_t row = 0; row < n_nbrs; row++) {
          candidate[row] -= projection * previous[row];
        }
      }
    }

    double norm = 0.0;
    for (std::size_t row = 0; row < n_nbrs; row++) {
      if (!std::isfinite(candidate[row])) {
        return -1;
      }
      norm = std::hypot(norm, candidate[row]);
    }
    if (norm <= drop_tolerance * std::max(1.0, original_norm)) {
      continue;
    }
    for (std::size_t row = 0; row < n_nbrs; row++) {
      candidate[row] /= norm;
    }
    basis_columns[n_retained] = basis_col;
    n_retained++;
  }

  const int dropped = static_cast<int>(n_selected - n_retained);
  basis_columns.resize(n_retained);
  return dropped;
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
[[noreturn]] void stop_local_weights_computation(int status, int neighborhood) {
  if (status == LOCAL_WEIGHTS_COMPUTATION_NONFINITE) {
    cpp11::stop("LTSA encountered non-finite neighborhood arithmetic during "
                "centering, scaling, or basis cleanup at neighborhood %d",
                neighborhood);
  }
  cpp11::stop("LTSA local weight computation failed at neighborhood %d",
              neighborhood);
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

int compute_local_weights_svd_workspace(const double *x_data, std::size_t n_obs,
                                        SvdLocalWeightsWorkspace &workspace,
                                        int &rank, int &computation_status) {
  rank = 0;
  computation_status = LOCAL_WEIGHTS_COMPUTATION_OK;

  if (!fill_centered_neighborhood_column_major(
          x_data, n_obs, workspace.neighbor_indices, workspace.centered,
          workspace.n_features_size)) {
    computation_status = LOCAL_WEIGHTS_COMPUTATION_NONFINITE;
    return 0;
  }
  std::copy(workspace.centered.begin(), workspace.centered.end(),
            workspace.a.begin());

  char jobz = 'S';
  int m = workspace.n_nbrs;
  int n = workspace.n_features;
  int lda = workspace.n_nbrs;
  int ldu = workspace.n_nbrs;
  int ldvt = workspace.min_dim;
  int lwork = static_cast<int>(workspace.work.size());
  int info = 0;
  F77_CALL(dgesdd)(&jobz, &m, &n, workspace.a.data(), &lda, workspace.d.data(),
                   workspace.u.data(), &ldu, workspace.vt.data(), &ldvt,
                   workspace.work.data(), &lwork, workspace.iwork.data(),
                   &info FCONE);
  if (info != 0) {
    return info;
  }

  rank = select_local_basis_columns(
      workspace.d, workspace.min_dim, workspace.n_nbrs, workspace.n_features,
      workspace.requested_basis_size, false, workspace.basis_columns);
  const int dropped = clean_local_basis(workspace.n_nbrs_size,
                                        workspace.basis_columns, workspace.u);
  if (dropped < 0) {
    computation_status = LOCAL_WEIGHTS_COMPUTATION_NONFINITE;
    return 0;
  }
  rank = std::max(0, rank - dropped);

  fill_weights_from_basis(workspace.n_nbrs_size, workspace.basis_columns,
                          workspace.u, workspace.weights);
  return 0;
}

int compute_local_weights_gram_workspace(const double *x_data,
                                         std::size_t n_obs,
                                         GramLocalWeightsWorkspace &workspace,
                                         const std::vector<double> *row_major,
                                         int &rank, int &computation_status) {
  rank = 0;
  computation_status = LOCAL_WEIGHTS_COMPUTATION_OK;
  bool centered_ok = false;
  if (row_major != nullptr) {
    centered_ok = fill_centered_neighborhood_row_major(
        *row_major, workspace.neighbor_indices, workspace.row_buffer,
        workspace.centered, workspace.n_features_size);
  } else {
    centered_ok = fill_centered_neighborhood_column_major(
        x_data, n_obs, workspace.neighbor_indices, workspace.centered,
        workspace.n_features_size);
  }
  if (!centered_ok) {
    computation_status = LOCAL_WEIGHTS_COMPUTATION_NONFINITE;
    return 0;
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
    return info;
  }

  rank = select_local_basis_columns(workspace.values, workspace.n_nbrs,
                                    workspace.n_nbrs, workspace.n_features,
                                    workspace.requested_basis_size, true,
                                    workspace.basis_columns);
  const int dropped = clean_local_basis(
      workspace.n_nbrs_size, workspace.basis_columns, workspace.gram);
  if (dropped < 0) {
    computation_status = LOCAL_WEIGHTS_COMPUTATION_NONFINITE;
    return 0;
  }
  rank = std::max(0, rank - dropped);

  fill_weights_from_basis(workspace.n_nbrs_size, workspace.basis_columns,
                          workspace.gram, workspace.weights);
  return 0;
}
