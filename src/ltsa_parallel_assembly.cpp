#include "ltsa_internal.h"

namespace {

struct ParallelLocalWeightsWorkspace {
  ParallelLocalWeightsWorkspace(std::size_t n_nbrs, std::size_t n_dim, int ndim,
                                bool use_svd, bool use_row_major)
      : n_nbrs_size(n_nbrs), n_dim_size(n_dim),
        n_nbrs(checked_lapack_dim(n_nbrs, "n_neighbors")),
        n_dim(checked_lapack_dim(n_dim, "ncol(X)")), route_svd(use_svd),
        min_dim(std::min(this->n_nbrs, this->n_dim)),
        requested(std::min(ndim, min_dim)),
        nni(checked_vector_size<int>(n_nbrs,
                                     "parallel LTSA neighborhood indices")),
        centered(checked_vector_size_mul<double>(
            n_nbrs, n_dim, "parallel LTSA centered neighborhood workspace")),
        weights(checked_vector_size_mul<double>(
            n_nbrs, n_nbrs, "parallel LTSA local weights")) {
    keep.reserve(
        checked_vector_size<int>(static_cast<std::size_t>(requested),
                                 "parallel LTSA retained local basis"));

    if (route_svd) {
      svd_a.resize(checked_vector_size_mul<double>(
          n_nbrs, n_dim, "parallel LTSA dgesdd matrix workspace"));
      d.resize(checked_vector_size<double>(static_cast<std::size_t>(min_dim),
                                           "parallel LTSA singular values"));
      u.resize(checked_vector_size_mul<double>(
          n_nbrs, static_cast<std::size_t>(min_dim),
          "parallel LTSA left singular vectors"));
      vt.resize(checked_vector_size_mul<double>(
          static_cast<std::size_t>(min_dim), n_dim,
          "parallel LTSA right singular vectors"));
      iwork.resize(checked_vector_size_mul<int>(
          8, static_cast<std::size_t>(min_dim),
          "parallel LTSA dgesdd integer workspace"));
      svd_work.resize(checked_vector_size<double>(
          static_cast<std::size_t>(query_dgesdd_workspace(
              this->n_nbrs, this->n_dim, min_dim, svd_a, d, u, vt, iwork)),
          "parallel LTSA dgesdd workspace"));
    } else {
      if (use_row_major) {
        row_buffer.resize(checked_vector_size_mul<double>(
            n_nbrs, n_dim, "parallel LTSA row-major neighborhood workspace"));
        col_means.resize(
            checked_vector_size<double>(n_dim, "parallel LTSA column means"));
      }
      gram.resize(checked_vector_size_mul<double>(
          n_nbrs, n_nbrs, "parallel LTSA Gram workspace"));
      values.resize(checked_vector_size<double>(
          n_nbrs, "parallel LTSA Gram eigenvalues"));
      gram_work.resize(checked_vector_size<double>(
          static_cast<std::size_t>(
              query_dsyev_workspace(this->n_nbrs, gram, values)),
          "parallel LTSA dsyev workspace"));
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
  int failed_step = 0;
  int failed_info = 0;
  int failed_obs = -1;
};

struct TriangularSlotPlan {
  std::vector<std::size_t> column_starts;
  std::vector<std::size_t> column_counts;
  std::vector<std::size_t> slot_offsets;
};

struct ReduceWorkspace {
  explicit ReduceWorkspace(std::size_t n_obs)
      : row_sums(checked_vector_size<double>(
                     n_obs, "parallel LTSA reduction row sums"),
                 0.0),
        row_seen(checked_vector_size<int>(
                     n_obs, "parallel LTSA reduction row markers"),
                 -1) {
    touched_rows.reserve(checked_vector_size<int>(
        std::min(n_obs, static_cast<std::size_t>(1024)),
        "parallel LTSA reduction touched rows"));
  }

  std::vector<double> row_sums;
  std::vector<int> row_seen;
  std::vector<int> touched_rows;
};

void fill_flat_neighbors_zero_based_ptr(const int* value_ptr,
                                        std::size_t offset, std::size_t n_nbrs,
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

  rank = select_local_basis_columns(workspace.d, workspace.min_dim,
                                    workspace.n_nbrs, workspace.n_dim,
                                    workspace.requested, false, workspace.keep);

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

  rank = select_local_basis_columns(workspace.values, workspace.n_nbrs,
                                    workspace.n_nbrs, workspace.n_dim,
                                    workspace.requested, true, workspace.keep);

  fill_weights_from_basis(workspace.n_nbrs_size, workspace.keep, workspace.gram,
                          workspace.weights);
  return 0;
}

int compute_parallel_local_weights(const double* x_data, std::size_t n_obs,
                                   ParallelLocalWeightsWorkspace& workspace,
                                   const std::vector<double>* row_major,
                                   int ndim,
                                   ParallelWorkerDiagnostics& diagnostics,
                                   std::size_t obs, int& rank) {
  rank = 0;
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

template <typename T>
void checked_resize_vector(std::vector<T>& out, std::size_t n,
                           const char* name) {
  checked_vector_size<T>(n, name);
  try {
    out.resize(n);
  } catch (const std::bad_alloc&) {
    cpp11::stop("Unable to allocate %s", name);
  } catch (const std::length_error&) {
    cpp11::stop("%s is too large", name);
  }
}

TriangularSlotPlan assign_triangular_two_pass_slots_flat(const int* value_ptr,
                                                         std::size_t n_obs,
                                                         std::size_t n_nbrs) {
  checked_triplet_count(n_obs, n_nbrs, "value_n_nbrs");

  TriangularSlotPlan plan;
  const std::size_t tri_count = triangular_pair_count(n_nbrs);
  plan.column_counts.assign(
      checked_vector_size<std::size_t>(n_obs, "triangular LTSA column counts"),
      0);
  checked_resize_vector(
      plan.slot_offsets,
      checked_size_mul(n_obs, tri_count,
                       "Too many triangular LTSA slot offsets"),
      "triangular LTSA slot offsets");

  std::vector<int> nni(
      checked_vector_size<int>(n_nbrs, "parallel LTSA neighborhood indices"));
  for (std::size_t obs = 0; obs < n_obs; obs++) {
    const std::size_t offset = obs * n_nbrs;

    for (std::size_t local = 0; local < n_nbrs; local++) {
      const int idx = checked_neighbor_index(value_ptr[offset + local], n_obs);
      nni[local] = idx;
    }

    const std::size_t obs_tri_offset = obs * tri_count;
    for (std::size_t local_col = 0; local_col < n_nbrs; local_col++) {
      for (std::size_t local_row = 0; local_row <= local_col; local_row++) {
        const int global_row = nni[local_row];
        const int global_col = nni[local_col];
        const int col = std::max(global_row, global_col);
        const std::size_t pair_offset =
            obs_tri_offset + triangular_pair_offset(local_col, local_row);
        plan.slot_offsets[pair_offset] = plan.column_counts[col];
        plan.column_counts[col] =
            checked_size_add(plan.column_counts[col], 1,
                             "Too many triangular LTSA contributions to stage");
      }
    }
  }

  plan.column_starts.assign(
      checked_vector_size<std::size_t>(
          checked_size_add(n_obs, 1, "Too many triangular LTSA column starts"),
          "triangular LTSA column starts"),
      0);
  for (std::size_t col = 0; col < n_obs; col++) {
    plan.column_starts[col + 1] =
        checked_size_add(plan.column_starts[col], plan.column_counts[col],
                         "Too many triangular LTSA contributions to stage");
  }
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
  const std::vector<std::size_t>* column_starts;
  const std::vector<std::size_t>* slot_offsets;
  std::vector<int>* raw_rows;
  std::vector<double>* raw_values;
  std::vector<int>* local_ranks;
  std::vector<ParallelLocalWeightsWorkspace>* workspaces;
  std::vector<ParallelWorkerDiagnostics>* diagnostics;

  void operator()(std::size_t begin, std::size_t end, std::size_t chunk_id) {
    ParallelLocalWeightsWorkspace& workspace = (*workspaces)[chunk_id];
    ParallelWorkerDiagnostics& worker_diagnostics = (*diagnostics)[chunk_id];

    for (std::size_t obs = begin; obs < end; obs++) {
      const std::size_t offset = obs * n_nbrs;
      fill_flat_neighbors_zero_based_ptr(value_ptr, offset, n_nbrs,
                                         workspace.nni);
      int local_rank = 0;
      if (compute_parallel_local_weights(x_data, n_obs, workspace, row_major,
                                         ndim, worker_diagnostics, obs,
                                         local_rank) != 0) {
        break;
      }
      (*local_ranks)[obs] = local_rank;

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
              (*column_starts)[col] + (*slot_offsets)[pair_offset];
          (*raw_rows)[pos] = row;
          (*raw_values)[pos] =
              workspace.weights[local_col * n_nbrs + local_row];
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

std::vector<std::vector<CompactEntry>>
reduce_raw_columns_parallel(const std::vector<std::size_t>& column_starts,
                            const std::vector<std::size_t>& column_counts,
                            const std::vector<int>& raw_rows,
                            const std::vector<double>& raw_values,
                            std::size_t n_obs, std::size_t n_threads) {
  checked_vector_size<CompactEntry>(n_obs, "parallel LTSA reduced column");
  std::vector<pforr::IndexRange> ranges =
      pforr::split_input_range(pforr::IndexRange(0, n_obs), n_threads, 1);
  std::vector<std::vector<CompactEntry>> reduced_columns(
      checked_vector_size<std::vector<CompactEntry>>(
          n_obs, "parallel LTSA reduced column containers"));
  std::vector<ReduceWorkspace> workspaces;
  workspaces.reserve(checked_vector_size<ReduceWorkspace>(
      ranges.size(), "parallel LTSA reduction workspaces"));
  for (std::size_t chunk = 0; chunk < ranges.size(); chunk++) {
    workspaces.emplace_back(n_obs);
  }

  ColumnReduceWorker worker{&column_starts, &column_counts,   &raw_rows,
                            &raw_values,    &reduced_columns, &workspaces};
  pforr::parallel_for_indexed(0, n_obs, worker, n_threads, 1);
  return reduced_columns;
}

void expand_canonical_columns_to_full(
    const std::vector<std::vector<CompactEntry>>& canonical_columns,
    std::vector<std::vector<CompactEntry>>& full_columns) {
  std::vector<std::size_t> full_column_counts(
      checked_vector_size<std::size_t>(canonical_columns.size(),
                                       "parallel LTSA full column counts"),
      0);
  for (std::size_t col = 0; col < canonical_columns.size(); col++) {
    for (const CompactEntry& entry : canonical_columns[col]) {
      full_column_counts[col] =
          checked_size_add(full_column_counts[col], 1,
                           "Too many parallel LTSA compact column entries");
      if (entry.row != static_cast<int>(col)) {
        full_column_counts[entry.row] =
            checked_size_add(full_column_counts[entry.row], 1,
                             "Too many parallel LTSA compact column entries");
      }
    }
  }
  for (std::size_t col = 0; col < full_columns.size(); col++) {
    full_columns[col].reserve(checked_vector_size<CompactEntry>(
        full_column_counts[col], "parallel LTSA full compact column"));
  }

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

SparseComponents
finalize_compact_columns(const std::vector<std::vector<CompactEntry>>& columns,
                         std::size_t n_obs, std::size_t max_int) {
  SparseComponents out;
  out.p.resize(
      checked_vector_size<int>(
          checked_size_add(n_obs, 1, "Too many parallel LTSA sparse columns"),
          "parallel LTSA sparse column pointers"),
      0);

  std::vector<double> row_sums(
      checked_vector_size<double>(n_obs, "parallel LTSA row sums"), 0.0);
  std::vector<int> row_seen(
      checked_vector_size<int>(n_obs, "parallel LTSA row markers"), -1);
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
    const TriangularSlotPlan& plan, const std::vector<int>& raw_rows,
    const std::vector<double>& raw_values, std::size_t n_obs,
    std::size_t n_threads, std::size_t max_int) {
  std::vector<std::vector<CompactEntry>> triangular_columns =
      reduce_raw_columns_parallel(plan.column_starts, plan.column_counts,
                                  raw_rows, raw_values, n_obs, n_threads);

  std::vector<std::vector<CompactEntry>> full_columns(
      checked_vector_size<std::vector<CompactEntry>>(
          n_obs, "parallel LTSA full column containers"));
  expand_canonical_columns_to_full(triangular_columns, full_columns);

  return finalize_compact_columns(full_columns, n_obs, max_int);
}

ParallelWorkerDiagnostics combine_worker_diagnostics(
    const std::vector<ParallelWorkerDiagnostics>& diagnostics, int ndim) {
  ParallelWorkerDiagnostics out;
  for (const ParallelWorkerDiagnostics& worker : diagnostics) {
    out.rank_deficient_count += worker.rank_deficient_count;
    out.min_local_rank = std::min(out.min_local_rank, worker.min_local_rank);
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
      cpp11::stop("LAPACK %s failed in LTSA assembly worker %d at neighborhood "
                  "%d with info = %d",
                  routine, static_cast<int>(worker + 1), current.failed_obs,
                  current.failed_info);
    }
  }
}

} // namespace

[[cpp11::register]] cpp11::list ltsa_assemble_local_weights_parallel(
    const cpp11::doubles_matrix<>& x, const cpp11::integers& value_nnt,
    std::size_t value_n_nbrs, int ndim, int requested_threads) {
  checked_ndim(ndim);
  if (requested_threads < 1) {
    cpp11::stop("n_assembly_threads must be positive");
  }
  if (value_nnt.size() == 0 || value_n_nbrs == 0) {
    cpp11::stop("Value neighborhoods must not be empty");
  }
  if (value_nnt.size() % value_n_nbrs != 0) {
    cpp11::stop("Inconsistent value neighborhood dimensions");
  }

  const std::size_t n_obs = value_nnt.size() / value_n_nbrs;
  if (static_cast<std::size_t>(x.nrow()) != n_obs) {
    cpp11::stop("Inconsistent input and neighborhood dimensions");
  }
  if (x.ncol() == 0) {
    cpp11::stop("X must contain at least one column");
  }

  const auto max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (n_obs >= max_int) {
    cpp11::stop("Too many observations for a dgCMatrix");
  }

  const int* value_ptr = INTEGER(value_nnt.data());
  TriangularSlotPlan slot_plan =
      assign_triangular_two_pass_slots_flat(value_ptr, n_obs, value_n_nbrs);

  const std::size_t requested_thread_count =
      static_cast<std::size_t>(requested_threads);
  const std::vector<pforr::IndexRange> obs_ranges = pforr::split_input_range(
      pforr::IndexRange(0, n_obs), requested_thread_count, 1);
  const std::size_t effective_threads = obs_ranges.size();

  const bool use_svd_route = static_cast<std::size_t>(x.ncol()) <= value_n_nbrs;
  const double* x_data = REAL(x.data());
  std::vector<double> row_major_x;
  const std::vector<double>* row_major_ptr = nullptr;
  if (!use_svd_route) {
    if (row_major_copy_within_limit(n_obs, static_cast<std::size_t>(x.ncol()),
                                    LTSA_ROW_MAJOR_COPY_MAX_BYTES)) {
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

  std::vector<ParallelWorkerDiagnostics> worker_diagnostics(
      checked_vector_size<ParallelWorkerDiagnostics>(
          effective_threads, "parallel LTSA worker diagnostics"));
  std::vector<ParallelLocalWeightsWorkspace> workspaces;
  workspaces.reserve(checked_vector_size<ParallelLocalWeightsWorkspace>(
      effective_threads, "parallel LTSA worker workspaces"));
  for (std::size_t chunk = 0; chunk < effective_threads; chunk++) {
    workspaces.emplace_back(value_n_nbrs, static_cast<std::size_t>(x.ncol()),
                            ndim, use_svd_route, row_major_ptr != nullptr);
  }

  const std::size_t raw_count = slot_plan.column_starts[n_obs];
  std::vector<int> raw_rows;
  std::vector<double> raw_values;
  std::vector<int> local_ranks(
      checked_vector_size<int>(n_obs, "parallel LTSA local ranks"), 0);
  checked_resize_vector(raw_rows, raw_count, "raw LTSA row buffer");
  checked_resize_vector(raw_values, raw_count, "raw LTSA value buffer");

  ParallelTriangularFillWorker worker{x_data,
                                      row_major_ptr,
                                      value_ptr,
                                      n_obs,
                                      value_n_nbrs,
                                      triangular_pair_count(value_n_nbrs),
                                      ndim,
                                      &slot_plan.column_starts,
                                      &slot_plan.slot_offsets,
                                      &raw_rows,
                                      &raw_values,
                                      &local_ranks,
                                      &workspaces,
                                      &worker_diagnostics};

  pforr::parallel_for_indexed(0, n_obs, worker, requested_thread_count, 1);
  stop_on_parallel_worker_failure(worker_diagnostics);

  SparseComponents components = finalize_triangular_two_pass_raw(
      slot_plan, raw_rows, raw_values, n_obs, requested_thread_count, max_int);
  ParallelWorkerDiagnostics diagnostics =
      combine_worker_diagnostics(worker_diagnostics, ndim);

  // Convert through const references before named_arg's by-value assignment so
  // the R payloads do not overlap temporary copies of the C++ vectors.
  cpp11::sexp r_i(cpp11::as_sexp(components.i));
  cpp11::sexp r_p(cpp11::as_sexp(components.p));
  cpp11::sexp r_x(cpp11::as_sexp(components.x));
  cpp11::sexp r_local_ranks(cpp11::as_sexp(local_ranks));

  return cpp11::writable::list(
      {cpp11::named_arg("i") = r_i, cpp11::named_arg("p") = r_p,
       cpp11::named_arg("x") = r_x,
       cpp11::named_arg("rank_deficient_count") =
           diagnostics.rank_deficient_count,
       cpp11::named_arg("min_local_rank") = diagnostics.min_local_rank,
       cpp11::named_arg("local_ranks") = r_local_ranks,
       cpp11::named_arg("local_solver_route") = use_svd_route ? "svd" : "gram",
       cpp11::named_arg("assembly_route") = "parallel_triangular_two_pass",
       cpp11::named_arg("requested_assembly_threads") = requested_threads,
       cpp11::named_arg("effective_assembly_threads") =
           static_cast<int>(effective_threads)});
}
