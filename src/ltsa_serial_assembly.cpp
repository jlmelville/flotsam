#include "ltsa_internal.h"

[[cpp11::register]] cpp11::list
assemble_local_weights(const cpp11::doubles_matrix<> &x,
                       const cpp11::integers &transposed_neighbor_indices,
                       std::size_t n_neighbors, int ndim) {
  checked_ndim(ndim);
  if (transposed_neighbor_indices.size() == 0 || n_neighbors == 0) {
    cpp11::stop("Value neighborhoods must not be empty");
  }
  if (transposed_neighbor_indices.size() % n_neighbors != 0) {
    cpp11::stop("Inconsistent value neighborhood dimensions");
  }

  std::size_t n_obs = transposed_neighbor_indices.size() / n_neighbors;
  if (static_cast<std::size_t>(x.nrow()) != n_obs) {
    cpp11::stop("Inconsistent input and neighborhood dimensions");
  }

  const auto max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (n_obs >= max_int) {
    cpp11::stop("Too many observations for a dgCMatrix");
  }

  TripletAssemblyBuilder builder(transposed_neighbor_indices, n_neighbors,
                                 n_obs, max_int);

  int rank_deficient_count = 0;
  int min_local_rank = ndim;
  const bool use_gram_workspace =
      x.ncol() != 0 && static_cast<std::size_t>(x.ncol()) > n_neighbors;
  const double *x_data = use_gram_workspace ? REAL(x.data()) : nullptr;
  std::vector<double> row_major_x;
  bool use_row_major_gram = false;
  std::unique_ptr<GramLocalWeightsWorkspace> gram_workspace;
  if (use_gram_workspace) {
    if (row_major_copy_within_limit(n_obs, static_cast<std::size_t>(x.ncol()),
                                    ROW_MAJOR_COPY_LIMIT_BYTES)) {
      try {
        make_row_major_copy(x_data, n_obs, static_cast<std::size_t>(x.ncol()),
                            row_major_x);
        use_row_major_gram = true;
      } catch (const std::bad_alloc &) {
        row_major_x.clear();
      } catch (const std::length_error &) {
        row_major_x.clear();
      }
    }
    gram_workspace.reset(new GramLocalWeightsWorkspace(
        n_neighbors, static_cast<std::size_t>(x.ncol()), ndim,
        use_row_major_gram));
  }

  for (std::size_t obs = 0; obs < n_obs; obs++) {
    const std::size_t offset = obs * n_neighbors;

    if (use_gram_workspace) {
      fill_flat_neighbors_zero_based(transposed_neighbor_indices, offset,
                                     n_neighbors,
                                     gram_workspace->neighbor_indices);
      int rank = compute_local_weights_gram_workspace(
          x_data, n_obs, *gram_workspace,
          use_row_major_gram ? &row_major_x : nullptr);
      if (rank < ndim) {
        rank_deficient_count++;
        min_local_rank = std::min(min_local_rank, rank);
      }
      builder.append_prechecked_neighborhood(gram_workspace->neighbor_indices,
                                             gram_workspace->weights);
    } else {
      std::vector<int> local_neighbor_indices = flat_neighbors_zero_based(
          transposed_neighbor_indices, offset, n_neighbors);
      LocalWeights local =
          compute_local_weights_by_shape(x, local_neighbor_indices, ndim);
      if (local.rank < ndim) {
        rank_deficient_count++;
        min_local_rank = std::min(min_local_rank, local.rank);
      }
      builder.append_prechecked_neighborhood(local_neighbor_indices,
                                             local.weights);
    }
  }

  SparseComponents components = builder.finalize_sparse_components();

  // Convert through const references before named_arg's by-value assignment so
  // the R payloads do not overlap temporary copies of the C++ vectors.
  cpp11::sexp r_i(cpp11::as_sexp(components.i));
  cpp11::sexp r_p(cpp11::as_sexp(components.p));
  cpp11::sexp r_x(cpp11::as_sexp(components.x));

  return cpp11::writable::list(
      {cpp11::named_arg("i") = r_i, cpp11::named_arg("p") = r_p,
       cpp11::named_arg("x") = r_x,
       cpp11::named_arg("rank_deficient_count") = rank_deficient_count,
       cpp11::named_arg("min_local_rank") = min_local_rank,
       cpp11::named_arg("assembly_route") = "serial_triangular",
       cpp11::named_arg("requested_assembly_threads") = 1,
       cpp11::named_arg("effective_assembly_threads") = 1});
}
