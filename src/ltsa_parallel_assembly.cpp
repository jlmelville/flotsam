#include "ltsa_internal.h"

#include <initializer_list>

namespace {

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
      svd_work.resize(query_dgesdd_workspace(this->n_nbrs, this->n_dim, min_dim,
                                             svd_a, d, u, vt, iwork));
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
  int failed_step = 0;
  int failed_info = 0;
  int failed_obs = -1;
};

struct TriangularSlotPlan {
  std::vector<std::size_t> column_starts;
  std::vector<std::size_t> column_counts;
  std::vector<std::size_t> slot_offsets;
  std::size_t raw_entries = 0;
  std::size_t raw_bytes = 0;
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

struct AssemblyRouteMemoryEstimate {
  std::size_t effective_workers = 1;
  std::size_t raw_entries = 0;
  std::size_t full_compact_slots_bound = 0;
  std::size_t final_sparse_slots_bound = 0;
  bool row_major_copy_included = false;
  std::size_t neighborhood_index_staging_bytes = 0;
  std::size_t sparse_slot_offsets_bytes = 0;
  std::size_t raw_row_staging_bytes = 0;
  std::size_t raw_value_staging_bytes = 0;
  std::size_t column_counters_bytes = 0;
  std::size_t column_starts_bytes = 0;
  std::size_t local_rank_staging_bytes = 0;
  std::size_t worker_workspace_bytes_each = 0;
  std::size_t worker_workspaces_bytes = 0;
  std::size_t reduction_workspaces_bytes_bound = 0;
  std::size_t column_containers_bytes = 0;
  std::size_t canonical_compact_staging_bytes_bound = 0;
  std::size_t full_compact_staging_bytes_bound = 0;
  std::size_t row_major_copy_bytes = 0;
  std::size_t final_sparse_output_bytes_bound = 0;
  std::size_t cpp_to_r_output_copy_bytes_bound = 0;
  std::size_t cpp_to_r_local_rank_copy_bytes = 0;
  std::size_t accepted_rank_diagnostics_bytes_bound = 0;
  std::size_t rank_diagnostics_workspaces_bytes_bound = 0;
  std::size_t accepted_component_diagnostics_bytes_bound = 0;
  std::size_t component_discovery_workspaces_bytes_bound = 0;
  std::size_t sparse_validation_staging_bytes_bound = 0;
  std::size_t reverse_occurrence_bytes = 0;
  std::size_t control_objects_bytes = 0;
  std::size_t fill_phase_bytes_bound = 0;
  std::size_t reduction_phase_bytes_bound = 0;
  std::size_t expansion_phase_bytes_bound = 0;
  std::size_t finalization_phase_bytes_bound = 0;
  std::size_t return_phase_bytes_bound = 0;
  std::size_t r_rank_diagnostics_phase_bytes_bound = 0;
  std::size_t r_component_diagnostics_phase_bytes_bound = 0;
  std::size_t r_sparse_validation_phase_bytes_bound = 0;
  std::size_t r_return_phase_bytes_bound = 0;
  std::size_t modeled_peak_bytes_bound = 0;
};

std::size_t memory_bytes(std::size_t count, std::size_t element_size) {
  return checked_size_mul(count, element_size,
                          "LTSA assembly memory estimate overflowed");
}

std::size_t memory_sum(std::initializer_list<std::size_t> values) {
  std::size_t total = 0;
  for (const std::size_t value : values) {
    total = checked_size_add(total, value,
                             "LTSA assembly memory estimate overflowed");
  }
  return total;
}

std::size_t serial_workspace_bytes(std::size_t n_nbrs, std::size_t n_dim,
                                   int ndim, bool use_svd, bool use_row_major) {
  const std::size_t min_dim = std::min(n_nbrs, n_dim);
  const std::size_t requested =
      std::min(static_cast<std::size_t>(ndim), min_dim);
  const std::size_t nbr_dim = checked_size_mul(
      n_nbrs, n_dim, "LTSA assembly memory estimate overflowed");
  const std::size_t nbr_squared = checked_size_mul(
      n_nbrs, n_nbrs, "LTSA assembly memory estimate overflowed");

  if (!use_svd) {
    const std::size_t fixed_bytes = memory_sum(
        {sizeof(GramLocalWeightsWorkspace), memory_bytes(n_nbrs, sizeof(int)),
         memory_bytes(requested, sizeof(int)),
         memory_bytes(nbr_dim, sizeof(double)),
         memory_bytes(nbr_squared, sizeof(double)),
         memory_bytes(n_nbrs, sizeof(double)),
         memory_bytes(nbr_squared, sizeof(double)),
         use_row_major ? memory_bytes(nbr_dim, sizeof(double)) : 0,
         use_row_major ? memory_bytes(n_dim, sizeof(double)) : 0});
    const int n_nbrs_int = checked_lapack_dim(n_nbrs, "n_neighbors");
    const std::size_t work =
        static_cast<std::size_t>(query_dsyev_workspace_size(n_nbrs_int));
    return memory_sum({fixed_bytes, memory_bytes(work, sizeof(double))});
  }

  const std::size_t vector_objects = memory_sum(
      {sizeof(LocalWeights), memory_bytes(6, sizeof(std::vector<double>)),
       memory_bytes(3, sizeof(std::vector<int>))});
  const std::size_t fixed_bytes = memory_sum(
      {vector_objects, memory_bytes(n_nbrs, sizeof(int)),
       memory_bytes(nbr_dim, sizeof(double)),
       memory_bytes(nbr_dim, sizeof(double)),
       memory_bytes(min_dim, sizeof(double)),
       memory_bytes(
           checked_size_mul(n_nbrs, min_dim,
                            "LTSA assembly memory estimate overflowed"),
           sizeof(double)),
       memory_bytes(
           checked_size_mul(min_dim, n_dim,
                            "LTSA assembly memory estimate overflowed"),
           sizeof(double)),
       memory_bytes(checked_size_mul(
                        8, min_dim, "LTSA assembly memory estimate overflowed"),
                    sizeof(int)),
       memory_bytes(requested, sizeof(int)),
       memory_bytes(nbr_squared, sizeof(double))});
  const int n_nbrs_int = checked_lapack_dim(n_nbrs, "n_neighbors");
  const int n_dim_int = checked_lapack_dim(n_dim, "ncol(X)");
  const int min_dim_int = checked_lapack_dim(min_dim, "local matrix rank");
  const std::size_t work = static_cast<std::size_t>(
      query_dgesdd_workspace_size(n_nbrs_int, n_dim_int, min_dim_int));
  return memory_sum({fixed_bytes, memory_bytes(work, sizeof(double))});
}

std::size_t parallel_workspace_bytes(std::size_t n_nbrs, std::size_t n_dim,
                                     int ndim, bool use_svd,
                                     bool use_row_major) {
  const std::size_t min_dim = std::min(n_nbrs, n_dim);
  const std::size_t requested =
      std::min(static_cast<std::size_t>(ndim), min_dim);
  const std::size_t nbr_dim = checked_size_mul(
      n_nbrs, n_dim, "LTSA assembly memory estimate overflowed");
  const std::size_t nbr_squared = checked_size_mul(
      n_nbrs, n_nbrs, "LTSA assembly memory estimate overflowed");
  std::size_t bytes = memory_sum({sizeof(ParallelLocalWeightsWorkspace),
                                  memory_bytes(n_nbrs, sizeof(int)),
                                  memory_bytes(requested, sizeof(int)),
                                  memory_bytes(nbr_dim, sizeof(double)),
                                  memory_bytes(nbr_squared, sizeof(double))});

  if (use_svd) {
    bytes = memory_sum(
        {bytes, memory_bytes(nbr_dim, sizeof(double)),
         memory_bytes(min_dim, sizeof(double)),
         memory_bytes(
             checked_size_mul(n_nbrs, min_dim,
                              "LTSA assembly memory estimate overflowed"),
             sizeof(double)),
         memory_bytes(
             checked_size_mul(min_dim, n_dim,
                              "LTSA assembly memory estimate overflowed"),
             sizeof(double)),
         memory_bytes(
             checked_size_mul(8, min_dim,
                              "LTSA assembly memory estimate overflowed"),
             sizeof(int))});
    const int n_nbrs_int = checked_lapack_dim(n_nbrs, "n_neighbors");
    const int n_dim_int = checked_lapack_dim(n_dim, "ncol(X)");
    const int min_dim_int = checked_lapack_dim(min_dim, "local matrix rank");
    const std::size_t work = static_cast<std::size_t>(
        query_dgesdd_workspace_size(n_nbrs_int, n_dim_int, min_dim_int));
    return memory_sum({bytes, memory_bytes(work, sizeof(double))});
  }

  bytes = memory_sum({bytes, memory_bytes(nbr_squared, sizeof(double)),
                      memory_bytes(n_nbrs, sizeof(double)),
                      use_row_major ? memory_bytes(nbr_dim, sizeof(double)) : 0,
                      use_row_major ? memory_bytes(n_dim, sizeof(double)) : 0});
  const int n_nbrs_int = checked_lapack_dim(n_nbrs, "n_neighbors");
  const std::size_t work =
      static_cast<std::size_t>(query_dsyev_workspace_size(n_nbrs_int));
  return memory_sum({bytes, memory_bytes(work, sizeof(double))});
}

std::size_t reduction_workspace_bytes(std::size_t n_obs) {
  const std::size_t touched_capacity =
      std::max(n_obs, static_cast<std::size_t>(1024));
  return memory_sum({sizeof(ReduceWorkspace),
                     memory_bytes(n_obs, sizeof(double)),
                     memory_bytes(n_obs, sizeof(int)),
                     memory_bytes(touched_capacity, sizeof(int))});
}

cpp11::writable::list
route_memory_estimate_list(const AssemblyRouteMemoryEstimate& estimate) {
  return cpp11::writable::list(
      {cpp11::named_arg("estimate_kind") = "modeled_storage_bound",
       cpp11::named_arg("effective_worker_count") =
           static_cast<double>(estimate.effective_workers),
       cpp11::named_arg("raw_entries") =
           static_cast<double>(estimate.raw_entries),
       cpp11::named_arg("full_compact_slots_bound") =
           static_cast<double>(estimate.full_compact_slots_bound),
       cpp11::named_arg("final_sparse_slots_bound") =
           static_cast<double>(estimate.final_sparse_slots_bound),
       cpp11::named_arg("row_major_copy_included") =
           estimate.row_major_copy_included,
       cpp11::named_arg("components_bytes") = cpp11::writable::list(
           {cpp11::named_arg("neighborhood_index_staging") =
                static_cast<double>(estimate.neighborhood_index_staging_bytes),
            cpp11::named_arg("sparse_slot_offsets") =
                static_cast<double>(estimate.sparse_slot_offsets_bytes),
            cpp11::named_arg("raw_row_staging") =
                static_cast<double>(estimate.raw_row_staging_bytes),
            cpp11::named_arg("raw_value_staging") =
                static_cast<double>(estimate.raw_value_staging_bytes),
            cpp11::named_arg("column_counters") =
                static_cast<double>(estimate.column_counters_bytes),
            cpp11::named_arg("column_starts") =
                static_cast<double>(estimate.column_starts_bytes),
            cpp11::named_arg("local_rank_staging") =
                static_cast<double>(estimate.local_rank_staging_bytes),
            cpp11::named_arg("worker_workspace_each_bound") =
                static_cast<double>(estimate.worker_workspace_bytes_each),
            cpp11::named_arg("worker_workspaces") =
                static_cast<double>(estimate.worker_workspaces_bytes),
            cpp11::named_arg("reduction_workspaces_bound") =
                static_cast<double>(estimate.reduction_workspaces_bytes_bound),
            cpp11::named_arg("compact_column_containers") =
                static_cast<double>(estimate.column_containers_bytes),
            cpp11::named_arg("canonical_compact_staging_bound") =
                static_cast<double>(
                    estimate.canonical_compact_staging_bytes_bound),
            cpp11::named_arg("full_compact_staging_bound") =
                static_cast<double>(estimate.full_compact_staging_bytes_bound),
            cpp11::named_arg("optional_row_major_copy") =
                static_cast<double>(estimate.row_major_copy_bytes),
            cpp11::named_arg("final_sparse_output_bound") =
                static_cast<double>(estimate.final_sparse_output_bytes_bound),
            cpp11::named_arg("cpp_to_r_output_copy_bound") =
                static_cast<double>(estimate.cpp_to_r_output_copy_bytes_bound),
            cpp11::named_arg("cpp_to_r_local_rank_copy") =
                static_cast<double>(estimate.cpp_to_r_local_rank_copy_bytes),
            cpp11::named_arg("accepted_rank_diagnostics_bound") =
                static_cast<double>(
                    estimate.accepted_rank_diagnostics_bytes_bound),
            cpp11::named_arg("rank_diagnostics_workspaces_bound") =
                static_cast<double>(
                    estimate.rank_diagnostics_workspaces_bytes_bound),
            cpp11::named_arg("accepted_component_diagnostics_bound") =
                static_cast<double>(
                    estimate.accepted_component_diagnostics_bytes_bound),
            cpp11::named_arg("component_discovery_workspaces_bound") =
                static_cast<double>(
                    estimate.component_discovery_workspaces_bytes_bound),
            cpp11::named_arg("sparse_validation_staging_bound") =
                static_cast<double>(
                    estimate.sparse_validation_staging_bytes_bound),
            cpp11::named_arg("reverse_occurrence") =
                static_cast<double>(estimate.reverse_occurrence_bytes),
            cpp11::named_arg("control_objects") =
                static_cast<double>(estimate.control_objects_bytes)}),
       cpp11::named_arg("phase_bytes_bound") = cpp11::writable::list(
           {cpp11::named_arg("fill") =
                static_cast<double>(estimate.fill_phase_bytes_bound),
            cpp11::named_arg("reduction") =
                static_cast<double>(estimate.reduction_phase_bytes_bound),
            cpp11::named_arg("expansion") =
                static_cast<double>(estimate.expansion_phase_bytes_bound),
            cpp11::named_arg("finalization") =
                static_cast<double>(estimate.finalization_phase_bytes_bound),
            cpp11::named_arg("return_copy") =
                static_cast<double>(estimate.return_phase_bytes_bound),
            cpp11::named_arg("r_rank_diagnostics") = static_cast<double>(
                estimate.r_rank_diagnostics_phase_bytes_bound),
            cpp11::named_arg("r_component_diagnostics") = static_cast<double>(
                estimate.r_component_diagnostics_phase_bytes_bound),
            cpp11::named_arg("r_sparse_validation") = static_cast<double>(
                estimate.r_sparse_validation_phase_bytes_bound),
            cpp11::named_arg("r_return") =
                static_cast<double>(estimate.r_return_phase_bytes_bound)}),
       cpp11::named_arg("modeled_peak_bytes_bound") =
           static_cast<double>(estimate.modeled_peak_bytes_bound)});
}

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
  if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
    cpp11::stop("%s is too large", name);
  }
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
  plan.column_counts.assign(n_obs, 0);
  checked_resize_vector(
      plan.slot_offsets,
      checked_size_mul(n_obs, tri_count,
                       "Too many triangular LTSA slot offsets"),
      "triangular LTSA slot offsets");

  std::vector<int> nni(n_nbrs);
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

  plan.column_starts.assign(n_obs + 1, 0);
  for (std::size_t col = 0; col < n_obs; col++) {
    plan.column_starts[col + 1] =
        checked_size_add(plan.column_starts[col], plan.column_counts[col],
                         "Too many triangular LTSA contributions to stage");
  }

  plan.raw_entries = plan.column_starts[n_obs];
  plan.raw_bytes = checked_raw_staging_bytes(plan.raw_entries);
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
  std::vector<pforr::IndexRange> ranges =
      pforr::split_input_range(pforr::IndexRange(0, n_obs), n_threads, 1);
  std::vector<std::vector<CompactEntry>> reduced_columns(n_obs);
  std::vector<ReduceWorkspace> workspaces;
  workspaces.reserve(ranges.size());
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
    const TriangularSlotPlan& plan, const std::vector<int>& raw_rows,
    const std::vector<double>& raw_values, std::size_t n_obs,
    std::size_t n_threads, std::size_t max_int) {
  std::vector<std::vector<CompactEntry>> triangular_columns =
      reduce_raw_columns_parallel(plan.column_starts, plan.column_counts,
                                  raw_rows, raw_values, n_obs, n_threads);

  std::vector<std::vector<CompactEntry>> full_columns(n_obs);
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

[[cpp11::register]] cpp11::list ltsa_assembly_memory_estimates_cpp(
    std::size_t n_obs, std::size_t n_nbrs, std::size_t n_dim, int ndim,
    int requested_threads, bool include_self, double row_major_copy_max_bytes) {
  checked_ndim(ndim);
  if (n_obs == 0 || n_nbrs == 0 || n_dim == 0) {
    cpp11::stop("LTSA assembly dimensions must be positive");
  }
  if (requested_threads < 1) {
    cpp11::stop("n_assembly_threads must be positive");
  }
  const std::size_t max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (n_obs >= max_int) {
    cpp11::stop("Too many observations for a dgCMatrix");
  }
  checked_lapack_dim(n_nbrs, "n_neighbors");
  checked_lapack_dim(n_dim, "ncol(X)");
  checked_triplet_count(n_obs, n_nbrs, "value_n_nbrs");

  const std::size_t row_major_copy_max =
      checked_row_major_copy_max_bytes(row_major_copy_max_bytes);
  const bool use_svd = n_dim <= n_nbrs;
  const bool row_major_copy_included =
      !use_svd && row_major_copy_within_limit(n_obs, n_dim, row_major_copy_max);
  const std::size_t row_major_copy_bytes =
      row_major_copy_included
          ? memory_bytes(
                checked_size_mul(n_obs, n_dim,
                                 "LTSA assembly memory estimate overflowed"),
                sizeof(double))
          : 0;

  const std::size_t neighbor_entries = checked_size_mul(
      n_obs, n_nbrs, "LTSA assembly memory estimate overflowed");
  const std::size_t neighborhood_index_staging_bytes =
      memory_bytes(checked_size_mul(neighbor_entries, include_self ? 1 : 2,
                                    "LTSA assembly memory estimate overflowed"),
                   sizeof(int));
  const std::size_t triangular_count = triangular_pair_count(n_nbrs);
  const std::size_t raw_entries = checked_size_mul(
      n_obs, triangular_count, "LTSA assembly memory estimate overflowed");
  const std::size_t doubled_raw_entries = checked_size_mul(
      raw_entries, 2, "LTSA assembly memory estimate overflowed");
  const std::size_t dense_slots = checked_size_mul(
      n_obs, n_obs, "LTSA assembly memory estimate overflowed");
  const std::size_t full_compact_slots_bound =
      std::min(doubled_raw_entries, dense_slots);
  const std::size_t final_sparse_slots_bound =
      std::min(full_compact_slots_bound, max_int);
  const std::size_t final_sparse_output_bytes_bound = memory_sum(
      {memory_bytes(checked_size_add(
                        n_obs, 1, "LTSA assembly memory estimate overflowed"),
                    sizeof(int)),
       memory_bytes(final_sparse_slots_bound, sizeof(int)),
       memory_bytes(final_sparse_slots_bound, sizeof(double))});
  const std::size_t local_rank_staging_bytes = memory_bytes(n_obs, sizeof(int));
  const std::size_t column_containers_bytes = memory_bytes(
      checked_size_mul(n_obs, 2, "LTSA assembly memory estimate overflowed"),
      sizeof(std::vector<CompactEntry>));
  const std::size_t canonical_compact_staging_bytes_bound =
      memory_bytes(raw_entries, sizeof(CompactEntry));
  const std::size_t full_compact_staging_bytes_bound =
      memory_bytes(full_compact_slots_bound, sizeof(CompactEntry));
  const std::size_t one_reduction_workspace = reduction_workspace_bytes(n_obs);
  const std::size_t max_local_rank = std::min(n_nbrs, n_dim);
  const std::size_t rank_histogram_entries = checked_size_add(
      max_local_rank, 1, "LTSA assembly memory estimate overflowed");
  const std::size_t accepted_rank_diagnostics_bytes_bound =
      memory_sum({memory_bytes(n_obs, sizeof(int)),
                  memory_bytes(rank_histogram_entries, sizeof(int)),
                  memory_bytes(rank_histogram_entries, sizeof(SEXP))});
  const std::size_t rank_diagnostics_workspaces_bytes_bound = memory_bytes(
      checked_size_mul(n_obs, 3, "LTSA assembly memory estimate overflowed"),
      sizeof(int));
  const std::size_t accepted_component_diagnostics_bytes_bound = memory_sum(
      {memory_bytes(checked_size_mul(
                        n_obs, 2, "LTSA assembly memory estimate overflowed"),
                    sizeof(int)),
       sizeof(int)});
  const std::size_t component_discovery_workspaces_bytes_bound = memory_bytes(
      checked_size_mul(n_obs, 4, "LTSA assembly memory estimate overflowed"),
      sizeof(int));
  const std::size_t sparse_validation_staging_bytes_bound = std::max(
      memory_bytes(checked_size_mul(n_obs, 2,
                                    "LTSA assembly memory estimate overflowed"),
                   sizeof(int)),
      memory_bytes(final_sparse_slots_bound, sizeof(int)));
  const std::size_t reverse_occurrence_bytes = memory_bytes(n_obs, sizeof(int));

  AssemblyRouteMemoryEstimate serial;
  serial.raw_entries = raw_entries;
  serial.full_compact_slots_bound = full_compact_slots_bound;
  serial.final_sparse_slots_bound = final_sparse_slots_bound;
  serial.row_major_copy_included = row_major_copy_included;
  serial.neighborhood_index_staging_bytes = neighborhood_index_staging_bytes;
  serial.column_counters_bytes = memory_bytes(n_obs, sizeof(std::size_t));
  serial.local_rank_staging_bytes = local_rank_staging_bytes;
  serial.worker_workspace_bytes_each = serial_workspace_bytes(
      n_nbrs, n_dim, ndim, use_svd, row_major_copy_included);
  serial.worker_workspaces_bytes = serial.worker_workspace_bytes_each;
  serial.reduction_workspaces_bytes_bound = one_reduction_workspace;
  serial.column_containers_bytes = column_containers_bytes;
  serial.canonical_compact_staging_bytes_bound =
      canonical_compact_staging_bytes_bound;
  serial.full_compact_staging_bytes_bound = full_compact_staging_bytes_bound;
  serial.row_major_copy_bytes = row_major_copy_bytes;
  serial.final_sparse_output_bytes_bound = final_sparse_output_bytes_bound;
  serial.cpp_to_r_output_copy_bytes_bound = final_sparse_output_bytes_bound;
  serial.cpp_to_r_local_rank_copy_bytes = local_rank_staging_bytes;
  serial.accepted_rank_diagnostics_bytes_bound =
      accepted_rank_diagnostics_bytes_bound;
  serial.rank_diagnostics_workspaces_bytes_bound =
      rank_diagnostics_workspaces_bytes_bound;
  serial.accepted_component_diagnostics_bytes_bound =
      accepted_component_diagnostics_bytes_bound;
  serial.component_discovery_workspaces_bytes_bound =
      component_discovery_workspaces_bytes_bound;
  serial.sparse_validation_staging_bytes_bound =
      sparse_validation_staging_bytes_bound;
  serial.reverse_occurrence_bytes = reverse_occurrence_bytes;
  serial.control_objects_bytes = memory_sum(
      {sizeof(LtsaTripletAssemblyBuilder), sizeof(SparseComponents)});
  const std::size_t serial_base = memory_sum(
      {serial.neighborhood_index_staging_bytes, serial.local_rank_staging_bytes,
       serial.worker_workspaces_bytes, serial.column_containers_bytes,
       serial.row_major_copy_bytes, serial.control_objects_bytes});
  serial.fill_phase_bytes_bound =
      memory_sum({serial_base, serial.column_counters_bytes,
                  serial.canonical_compact_staging_bytes_bound});
  serial.reduction_phase_bytes_bound =
      memory_sum({serial_base, serial.canonical_compact_staging_bytes_bound,
                  serial.full_compact_staging_bytes_bound,
                  serial.reduction_workspaces_bytes_bound});
  serial.expansion_phase_bytes_bound = serial.reduction_phase_bytes_bound;
  serial.finalization_phase_bytes_bound =
      memory_sum({serial_base, serial.full_compact_staging_bytes_bound,
                  serial.reduction_workspaces_bytes_bound,
                  serial.final_sparse_output_bytes_bound});
  serial.return_phase_bytes_bound =
      memory_sum({serial_base, serial.final_sparse_output_bytes_bound,
                  serial.cpp_to_r_output_copy_bytes_bound,
                  serial.cpp_to_r_local_rank_copy_bytes});
  const std::size_t serial_r_base = memory_sum(
      {serial.neighborhood_index_staging_bytes, serial.local_rank_staging_bytes,
       serial.final_sparse_output_bytes_bound});
  serial.r_rank_diagnostics_phase_bytes_bound =
      memory_sum({serial_r_base, serial.accepted_rank_diagnostics_bytes_bound,
                  serial.rank_diagnostics_workspaces_bytes_bound});
  serial.r_component_diagnostics_phase_bytes_bound =
      memory_sum({serial_r_base, serial.accepted_rank_diagnostics_bytes_bound,
                  serial.accepted_component_diagnostics_bytes_bound,
                  serial.component_discovery_workspaces_bytes_bound});
  serial.r_sparse_validation_phase_bytes_bound =
      memory_sum({serial_r_base, serial.accepted_rank_diagnostics_bytes_bound,
                  serial.accepted_component_diagnostics_bytes_bound,
                  serial.sparse_validation_staging_bytes_bound});
  serial.r_return_phase_bytes_bound =
      memory_sum({serial_r_base, serial.accepted_rank_diagnostics_bytes_bound,
                  serial.accepted_component_diagnostics_bytes_bound,
                  serial.reverse_occurrence_bytes});
  serial.modeled_peak_bytes_bound = std::max(
      {serial.fill_phase_bytes_bound, serial.reduction_phase_bytes_bound,
       serial.expansion_phase_bytes_bound,
       serial.finalization_phase_bytes_bound, serial.return_phase_bytes_bound,
       serial.r_rank_diagnostics_phase_bytes_bound,
       serial.r_component_diagnostics_phase_bytes_bound,
       serial.r_sparse_validation_phase_bytes_bound,
       serial.r_return_phase_bytes_bound});

  AssemblyRouteMemoryEstimate parallel;
  parallel.effective_workers =
      std::min(n_obs, static_cast<std::size_t>(requested_threads));
  parallel.raw_entries = raw_entries;
  parallel.full_compact_slots_bound = full_compact_slots_bound;
  parallel.final_sparse_slots_bound = final_sparse_slots_bound;
  parallel.row_major_copy_included = row_major_copy_included;
  parallel.neighborhood_index_staging_bytes = neighborhood_index_staging_bytes;
  parallel.sparse_slot_offsets_bytes =
      memory_bytes(raw_entries, sizeof(std::size_t));
  parallel.raw_row_staging_bytes = memory_bytes(raw_entries, sizeof(int));
  parallel.raw_value_staging_bytes = memory_bytes(raw_entries, sizeof(double));
  parallel.column_counters_bytes = memory_bytes(n_obs, sizeof(std::size_t));
  parallel.column_starts_bytes = memory_bytes(
      checked_size_add(n_obs, 1, "LTSA assembly memory estimate overflowed"),
      sizeof(std::size_t));
  parallel.local_rank_staging_bytes = local_rank_staging_bytes;
  parallel.worker_workspace_bytes_each = parallel_workspace_bytes(
      n_nbrs, n_dim, ndim, use_svd, row_major_copy_included);
  parallel.worker_workspaces_bytes = checked_size_mul(
      parallel.worker_workspace_bytes_each, parallel.effective_workers,
      "LTSA assembly memory estimate overflowed");
  parallel.reduction_workspaces_bytes_bound =
      checked_size_mul(one_reduction_workspace, parallel.effective_workers,
                       "LTSA assembly memory estimate overflowed");
  parallel.column_containers_bytes = column_containers_bytes;
  parallel.canonical_compact_staging_bytes_bound =
      canonical_compact_staging_bytes_bound;
  parallel.full_compact_staging_bytes_bound = full_compact_staging_bytes_bound;
  parallel.row_major_copy_bytes = row_major_copy_bytes;
  parallel.final_sparse_output_bytes_bound = final_sparse_output_bytes_bound;
  parallel.cpp_to_r_output_copy_bytes_bound = final_sparse_output_bytes_bound;
  parallel.cpp_to_r_local_rank_copy_bytes = local_rank_staging_bytes;
  parallel.accepted_rank_diagnostics_bytes_bound =
      accepted_rank_diagnostics_bytes_bound;
  parallel.rank_diagnostics_workspaces_bytes_bound =
      rank_diagnostics_workspaces_bytes_bound;
  parallel.accepted_component_diagnostics_bytes_bound =
      accepted_component_diagnostics_bytes_bound;
  parallel.component_discovery_workspaces_bytes_bound =
      component_discovery_workspaces_bytes_bound;
  parallel.sparse_validation_staging_bytes_bound =
      sparse_validation_staging_bytes_bound;
  parallel.reverse_occurrence_bytes = reverse_occurrence_bytes;
  parallel.control_objects_bytes = memory_sum(
      {sizeof(TriangularSlotPlan), sizeof(SparseComponents),
       memory_bytes(parallel.effective_workers,
                    sizeof(ParallelWorkerDiagnostics)),
       memory_bytes(
           checked_size_mul(parallel.effective_workers, 3,
                            "LTSA assembly memory estimate overflowed"),
           sizeof(pforr::IndexRange)),
       memory_bytes(parallel.effective_workers, sizeof(std::thread))});
  const std::size_t parallel_base = memory_sum(
      {parallel.neighborhood_index_staging_bytes,
       parallel.sparse_slot_offsets_bytes, parallel.column_counters_bytes,
       parallel.column_starts_bytes, parallel.local_rank_staging_bytes,
       parallel.worker_workspaces_bytes, parallel.row_major_copy_bytes,
       parallel.control_objects_bytes});
  const std::size_t parallel_raw_staging = memory_sum(
      {parallel.raw_row_staging_bytes, parallel.raw_value_staging_bytes});
  parallel.fill_phase_bytes_bound =
      memory_sum({parallel_base, parallel_raw_staging});
  parallel.reduction_phase_bytes_bound =
      memory_sum({parallel_base, parallel_raw_staging,
                  parallel.reduction_workspaces_bytes_bound,
                  parallel.column_containers_bytes,
                  parallel.canonical_compact_staging_bytes_bound});
  parallel.expansion_phase_bytes_bound = memory_sum(
      {parallel_base, parallel_raw_staging, parallel.column_containers_bytes,
       parallel.canonical_compact_staging_bytes_bound,
       parallel.full_compact_staging_bytes_bound});
  parallel.finalization_phase_bytes_bound = memory_sum(
      {parallel_base, parallel_raw_staging, parallel.column_containers_bytes,
       parallel.canonical_compact_staging_bytes_bound,
       parallel.full_compact_staging_bytes_bound, one_reduction_workspace,
       parallel.final_sparse_output_bytes_bound});
  parallel.return_phase_bytes_bound =
      memory_sum({parallel_base, parallel_raw_staging,
                  parallel.final_sparse_output_bytes_bound,
                  parallel.cpp_to_r_output_copy_bytes_bound,
                  parallel.cpp_to_r_local_rank_copy_bytes});
  const std::size_t parallel_r_base =
      memory_sum({parallel.neighborhood_index_staging_bytes,
                  parallel.local_rank_staging_bytes,
                  parallel.final_sparse_output_bytes_bound});
  parallel.r_rank_diagnostics_phase_bytes_bound = memory_sum(
      {parallel_r_base, parallel.accepted_rank_diagnostics_bytes_bound,
       parallel.rank_diagnostics_workspaces_bytes_bound});
  parallel.r_component_diagnostics_phase_bytes_bound = memory_sum(
      {parallel_r_base, parallel.accepted_rank_diagnostics_bytes_bound,
       parallel.accepted_component_diagnostics_bytes_bound,
       parallel.component_discovery_workspaces_bytes_bound});
  parallel.r_sparse_validation_phase_bytes_bound = memory_sum(
      {parallel_r_base, parallel.accepted_rank_diagnostics_bytes_bound,
       parallel.accepted_component_diagnostics_bytes_bound,
       parallel.sparse_validation_staging_bytes_bound});
  parallel.r_return_phase_bytes_bound = memory_sum(
      {parallel_r_base, parallel.accepted_rank_diagnostics_bytes_bound,
       parallel.accepted_component_diagnostics_bytes_bound,
       parallel.reverse_occurrence_bytes});
  parallel.modeled_peak_bytes_bound = std::max(
      {parallel.fill_phase_bytes_bound, parallel.reduction_phase_bytes_bound,
       parallel.expansion_phase_bytes_bound,
       parallel.finalization_phase_bytes_bound,
       parallel.return_phase_bytes_bound,
       parallel.r_rank_diagnostics_phase_bytes_bound,
       parallel.r_component_diagnostics_phase_bytes_bound,
       parallel.r_sparse_validation_phase_bytes_bound,
       parallel.r_return_phase_bytes_bound});

  return cpp11::writable::list(
      {cpp11::named_arg("serial") = route_memory_estimate_list(serial),
       cpp11::named_arg("parallel") = route_memory_estimate_list(parallel),
       cpp11::named_arg("sizeof_bytes") = cpp11::writable::list(
           {cpp11::named_arg("int") = static_cast<double>(sizeof(int)),
            cpp11::named_arg("double") = static_cast<double>(sizeof(double)),
            cpp11::named_arg("size_t") =
                static_cast<double>(sizeof(std::size_t)),
            cpp11::named_arg("r_sexp_pointer") =
                static_cast<double>(sizeof(SEXP)),
            cpp11::named_arg("compact_entry") =
                static_cast<double>(sizeof(CompactEntry)),
            cpp11::named_arg("compact_column_container") =
                static_cast<double>(sizeof(std::vector<CompactEntry>)),
            cpp11::named_arg("serial_gram_workspace_object") =
                static_cast<double>(sizeof(GramLocalWeightsWorkspace)),
            cpp11::named_arg("parallel_workspace_object") =
                static_cast<double>(sizeof(ParallelLocalWeightsWorkspace)),
            cpp11::named_arg("reduction_workspace_object") =
                static_cast<double>(sizeof(ReduceWorkspace)),
            cpp11::named_arg("thread_object") =
                static_cast<double>(sizeof(std::thread))})});
}

[[cpp11::register]] cpp11::list ltsa_assemble_local_weights_parallel(
    const cpp11::doubles_matrix<>& x, const cpp11::integers& value_nnt,
    std::size_t value_n_nbrs, int ndim, int requested_threads,
    double row_major_copy_max_bytes) {
  checked_ndim(ndim);
  const std::size_t row_major_copy_max =
      checked_row_major_copy_max_bytes(row_major_copy_max_bytes);
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
  if (n_obs + 1 > max_int) {
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
  bool row_major_within_limit = false;
  if (!use_svd_route) {
    row_major_within_limit = row_major_copy_within_limit(
        n_obs, static_cast<std::size_t>(x.ncol()), row_major_copy_max);
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

  const std::size_t raw_count = slot_plan.column_starts[n_obs];
  std::vector<int> raw_rows;
  std::vector<double> raw_values;
  std::vector<int> local_ranks(n_obs, 0);
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
  const std::string fallback_reason = row_major_fallback_reason(
      !use_svd_route, row_major_ptr != nullptr, row_major_within_limit);

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
           static_cast<int>(effective_threads),
       cpp11::named_arg("raw_entries_estimate") =
           static_cast<double>(slot_plan.raw_entries),
       cpp11::named_arg("raw_bytes_estimate") =
           static_cast<double>(slot_plan.raw_bytes),
       cpp11::named_arg("row_major_used") = row_major_ptr != nullptr,
       cpp11::named_arg("row_major_fallback_reason") = fallback_reason,
       cpp11::named_arg("parallel_fallback_reason") = ""});
}
