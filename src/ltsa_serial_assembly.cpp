#include "ltsa_internal.h"

[[cpp11::register]] list
ltsa_assemble_local_weights(const doubles_matrix<>& x,
                            const integers& value_nnt, std::size_t value_n_nbrs,
                            int ndim, double row_major_copy_max_bytes) {
  checked_ndim(ndim);
  const std::size_t row_major_copy_max =
      checked_row_major_copy_max_bytes(row_major_copy_max_bytes);
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
    row_major_within_limit = row_major_copy_within_limit(
        n_obs, static_cast<std::size_t>(x.ncol()), row_major_copy_max);
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

  return writable::list(
      {"i"_nm = components.i, "p"_nm = components.p, "x"_nm = components.x,
       "rank_deficient_count"_nm = rank_deficient_count,
       "min_local_rank"_nm = min_local_rank,
       "assembly_route"_nm = "serial_triangular",
       "requested_assembly_threads"_nm = 1, "effective_assembly_threads"_nm = 1,
       "raw_entries_estimate"_nm = static_cast<double>(raw_entries),
       "raw_bytes_estimate"_nm = static_cast<double>(raw_bytes),
       "duplicate_fallback_count"_nm =
           static_cast<int>(builder.duplicate_fallback_count()),
       "row_major_used"_nm = use_row_major_gram,
       "row_major_fallback_reason"_nm = fallback_reason,
       "parallel_fallback_reason"_nm = "not_requested"});
}
