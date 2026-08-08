expect_sparse_equivalent <- function(candidate, reference, tolerance = 1e-12) {
  candidate <- Matrix::drop0(candidate)
  reference <- Matrix::drop0(reference)
  expect_s4_class(candidate, "dgCMatrix")
  expect_s4_class(reference, "dgCMatrix")
  expect_identical(candidate@Dim, reference@Dim)
  expect_equal(
    as.matrix(candidate),
    as.matrix(reference),
    tolerance = tolerance
  )
}

expect_ltsa_assembly_parallel_matches <- function(
  X,
  nn_idx,
  ndim,
  include_self,
  n_assembly_threads = 3L,
  tolerance = 1e-11
) {
  serial <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = ndim,
    include_self = include_self,
    n_assembly_threads = 1L
  )
  parallel <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = ndim,
    include_self = include_self,
    n_assembly_threads = n_assembly_threads
  )

  expect_sparse_equivalent(parallel$B, serial$B, tolerance = tolerance)
  expect_identical(
    parallel$rank_deficient_count,
    serial$rank_deficient_count
  )
  expect_identical(parallel$min_local_rank, serial$min_local_rank)
  expect_identical(
    parallel$diagnostics$local_solver_route,
    serial$diagnostics$local_solver_route
  )
  expect_identical(
    parallel$diagnostics$local_rank_histogram,
    serial$diagnostics$local_rank_histogram
  )
  expect_identical(
    parallel$diagnostics$rank_deficient_neighborhood_indices,
    serial$diagnostics$rank_deficient_neighborhood_indices
  )
  expect_identical(
    parallel$diagnostics$component_count,
    serial$diagnostics$component_count
  )
  expect_identical(
    parallel$diagnostics$component_sizes,
    serial$diagnostics$component_sizes
  )
  expect_identical(
    parallel$diagnostics$component_membership,
    serial$diagnostics$component_membership
  )
  expect_identical(serial$diagnostics$assembly_route, "serial_triangular")
  expect_identical(
    parallel$diagnostics$assembly_route,
    "parallel_triangular_two_pass"
  )
  expect_identical(
    parallel$diagnostics$requested_assembly_threads,
    as.integer(n_assembly_threads)
  )
  expect_false("raw_bytes_estimate" %in% names(serial$diagnostics))
  expect_true("raw_bytes_estimate" %in% names(parallel$diagnostics))
  expect_identical(parallel$diagnostics$memory$requested_route, "parallel")
  expect_identical(parallel$diagnostics$memory$selected_route, "parallel")
  expect_true(Matrix::isSymmetric(parallel$B))
  expect_equal(sum(parallel$B@x == 0), 0)

  invisible(parallel)
}

expect_embedding_subspace_equivalent <- function(
  candidate,
  reference,
  tolerance = 1e-8
) {
  expect_same_subspace(candidate, reference, tolerance = tolerance)
}

exact_nn_idx <- function(X, n_neighbors, include_self) {
  nn <- rnndescent::brute_force_knn(
    data = X,
    k = ifelse(include_self, n_neighbors, n_neighbors + 1L),
    n_threads = 0L,
    verbose = FALSE
  )
  mode(nn$idx) <- "integer"
  nn$idx
}

ltsa_local_weights_r_reference <- function(X, nni, ndim) {
  Xi <- scale(X[nni, , drop = FALSE], center = TRUE, scale = FALSE)
  max_rank <- min(dim(Xi))
  ndim <- min(ndim, max_rank)
  res <- svd(Xi, nu = ndim, nv = 0)

  if (length(res$d) == 0 || max(res$d) == 0) {
    rank <- 0L
    tol <- 0
  } else {
    tol <- max(dim(Xi)) * max(res$d) * .Machine$double.eps
    rank <- sum(res$d > tol)
  }

  keep <- seq_len(min(ndim, length(res$d), ncol(res$u)))
  keep <- keep[res$d[keep] > tol]

  k <- length(nni)
  Gi <- cbind(1 / sqrt(k), res$u[, keep, drop = FALSE])
  Wi <- -tcrossprod(Gi)
  diag(Wi) <- diag(Wi) + 1.0

  list(Wi = Wi, rank = rank)
}

assemble_ltsa_B_r_triplet_reference <- function(X, nn_idx, ndim, include_self) {
  n <- nrow(X)
  weight_idx <- if (include_self) {
    nn_idx
  } else {
    nn_idx[, -1L, drop = FALSE]
  }
  k <- ncol(weight_idx)
  n_triplets <- n * k * k
  rows <- integer(n_triplets)
  cols <- integer(n_triplets)
  vals <- numeric(n_triplets)
  rank_deficient_count <- 0L
  min_local_rank <- ndim

  for (i in seq_len(n)) {
    nni <- weight_idx[i, ]
    local <- ltsa_local_weights_r_reference(X, nni, ndim)
    if (local$rank < ndim) {
      rank_deficient_count <- rank_deficient_count + 1L
      min_local_rank <- min(min_local_rank, local$rank)
    }

    idx <- ((i - 1L) * k * k + 1L):(i * k * k)
    rows[idx] <- rep.int(nni, times = k)
    cols[idx] <- rep(nni, each = k)
    vals[idx] <- as.vector(local$Wi)
  }

  list(
    B = Matrix::sparseMatrix(
      i = rows,
      j = cols,
      x = vals,
      dims = c(n, n),
      giveCsparse = TRUE
    ),
    rank_deficient_count = rank_deficient_count,
    min_local_rank = min_local_rank
  )
}

test_that("assembly verbose logging suppresses non-actionable fallback reasons", {
  expect_false(flotsam:::ltsa_log_fallback_reason(""))
  expect_false(flotsam:::ltsa_log_fallback_reason("not_requested"))
  expect_false(flotsam:::ltsa_log_fallback_reason("not_applicable_svd_route"))
  expect_true(flotsam:::ltsa_log_fallback_reason("copy_size_exceeds_limit"))
})

test_that("default assembly behavior remains serial", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  default <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  explicit <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 1L
  )

  expect_sparse_equivalent(default$B, explicit$B, tolerance = 0)
  expect_identical(default$diagnostics$assembly_route, "serial_triangular")
  expect_identical(
    default$diagnostics$requested_assembly_threads,
    1L
  )
  expect_identical(
    default$diagnostics$effective_assembly_threads,
    1L
  )
  default_ltsa <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B"
  )
  explicit_ltsa <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_assembly_threads = 1L
  )
  expect_sparse_equivalent(default_ltsa, explicit_ltsa, tolerance = 0)
})

test_that("n_threads remains nearest-neighbor-only", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])

  serial_nn0 <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_threads = 0L,
    n_assembly_threads = 1L
  )
  serial_nn2 <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_threads = 2L,
    n_assembly_threads = 1L
  )

  expect_sparse_equivalent(serial_nn2, serial_nn0, tolerance = 0)
})

test_that("public output modes control B and detailed result contents", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])

  B <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B"
  )
  expect_s4_class(B, "dgCMatrix")

  embedding <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    eig_method = "eig",
    include_B = TRUE
  )
  expect_true(is.matrix(embedding))

  result_without_B <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    eig_method = "eig",
    eig_k = 4L,
    output = "result"
  )
  expect_true(is.list(result_without_B))
  expect_false("B" %in% names(result_without_B))
  expect_false(
    "duplicate_fallback_count" %in% names(result_without_B$assembly)
  )

  result_with_B <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    eig_method = "eig",
    eig_k = 4L,
    output = "result",
    include_B = TRUE
  )
  expect_sparse_equivalent(result_with_B$B, B, tolerance = 0)
  expect_named(result_with_B, c("embedding", "eigen", "assembly", "B"))
  expect_equal(dim(result_with_B$embedding), c(nrow(X), 2L))
  expect_identical(result_with_B$assembly$n_neighbors, 6L)
  expect_true(result_with_B$assembly$include_self)
})

test_that("parallel assembly matches serial on high-p Gram route", {
  set.seed(10)
  X <- matrix(rnorm(20L * 14L), nrow = 20L)

  for (include_self in c(TRUE, FALSE)) {
    nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = include_self)
    parallel <- expect_ltsa_assembly_parallel_matches(
      X,
      nn_idx,
      ndim = 2L,
      include_self = include_self,
      n_assembly_threads = 4L
    )
    expect_true(parallel$diagnostics$row_major_used)
  }
})

test_that("copy_max_mib controls optional row-major Gram copy", {
  set.seed(101)
  X <- matrix(rnorm(18L * 14L), nrow = 18L)
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  serial_default <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 1L
  )
  serial_disabled <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 1L,
    copy_max_mib = 0
  )
  parallel_default <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 3L
  )
  parallel_disabled <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 3L,
    copy_max_mib = 0
  )

  expect_true(serial_default$diagnostics$row_major_used)
  expect_false(serial_disabled$diagnostics$row_major_used)
  expect_identical(
    serial_disabled$diagnostics$row_major_fallback_reason,
    "copy_size_exceeds_limit"
  )
  expect_true(parallel_default$diagnostics$row_major_used)
  expect_false(parallel_disabled$diagnostics$row_major_used)
  expect_identical(
    parallel_disabled$diagnostics$row_major_fallback_reason,
    "copy_size_exceeds_limit"
  )
  expect_sparse_equivalent(
    serial_disabled$B,
    serial_default$B,
    tolerance = 1e-12
  )
  expect_sparse_equivalent(
    parallel_disabled$B,
    serial_default$B,
    tolerance = 1e-12
  )
})

test_that("assembly estimator reports checked route components without allocation", {
  estimates <- flotsam:::ltsa_assembly_memory_estimates(
    n = 10,
    k = 4,
    p = 3,
    ndim = 2L,
    requested_assembly_threads = 3L,
    include_self = TRUE,
    row_major_copy_max_bytes = 256 * 1024^2
  )

  expect_identical(estimates$serial$estimate_kind, "modeled_storage_bound")
  expect_identical(estimates$parallel$estimate_kind, "modeled_storage_bound")
  expect_equal(estimates$parallel$raw_entries, 100)
  expect_equal(
    estimates$parallel$components_bytes$sparse_slot_offsets,
    100 * estimates$sizeof_bytes$size_t
  )
  expect_equal(
    estimates$parallel$components_bytes$raw_row_staging,
    100 * estimates$sizeof_bytes$int
  )
  expect_equal(
    estimates$parallel$components_bytes$raw_value_staging,
    100 * estimates$sizeof_bytes$double
  )
  expect_equal(estimates$serial$components_bytes$sparse_slot_offsets, 0)
  expect_equal(estimates$serial$components_bytes$raw_row_staging, 0)
  expect_equal(estimates$serial$components_bytes$raw_value_staging, 0)
  expect_equal(estimates$serial$full_compact_slots_bound, 100)
  expect_equal(estimates$serial$final_sparse_slots_bound, 100)
  expect_equal(
    estimates$serial$components_bytes$final_sparse_output_bound,
    11 *
      estimates$sizeof_bytes$int +
      100 * (estimates$sizeof_bytes$int + estimates$sizeof_bytes$double)
  )
  expect_equal(
    estimates$serial$components_bytes$accepted_rank_diagnostics_bound,
    10 *
      estimates$sizeof_bytes$int +
      4 * estimates$sizeof_bytes$int +
      4 * estimates$sizeof_bytes$r_sexp_pointer
  )
  expect_equal(
    estimates$serial$components_bytes$rank_diagnostics_workspaces_bound,
    30 * estimates$sizeof_bytes$int
  )
  expect_equal(
    estimates$serial$components_bytes$accepted_component_diagnostics_bound,
    21 * estimates$sizeof_bytes$int
  )
  expect_equal(
    estimates$serial$components_bytes$component_discovery_workspaces_bound,
    40 * estimates$sizeof_bytes$int
  )
  expect_equal(
    estimates$serial$components_bytes$sparse_validation_staging_bound,
    100 * estimates$sizeof_bytes$int
  )
  expect_equal(
    estimates$serial$components_bytes$reverse_occurrence,
    10 * estimates$sizeof_bytes$int
  )
  expect_named(
    estimates$serial$phase_bytes_bound,
    c(
      "fill",
      "reduction",
      "expansion",
      "finalization",
      "return_copy",
      "r_rank_diagnostics",
      "r_component_diagnostics",
      "r_sparse_validation",
      "r_return"
    )
  )
  components <- estimates$serial$components_bytes
  r_base <- components$neighborhood_index_staging +
    components$local_rank_staging +
    components$final_sparse_output_bound
  expect_equal(
    estimates$serial$phase_bytes_bound$r_component_diagnostics,
    r_base +
      components$accepted_rank_diagnostics_bound +
      components$accepted_component_diagnostics_bound +
      components$component_discovery_workspaces_bound
  )
  expect_equal(
    estimates$serial$phase_bytes_bound$r_sparse_validation,
    r_base +
      components$accepted_rank_diagnostics_bound +
      components$accepted_component_diagnostics_bound +
      components$sparse_validation_staging_bound
  )
  expect_equal(
    estimates$serial$phase_bytes_bound$r_return,
    r_base +
      components$accepted_rank_diagnostics_bound +
      components$accepted_component_diagnostics_bound +
      components$reverse_occurrence
  )
  expect_equal(
    estimates$parallel$modeled_peak_bytes_bound,
    max(unlist(estimates$parallel$phase_bytes_bound))
  )
})

test_that("assembly estimator keeps compact and sparse slot bounds distinct", {
  args <- list(
    n = 50000,
    k = 250,
    p = 1,
    ndim = 1L,
    requested_assembly_threads = 2L,
    include_self = TRUE,
    row_major_copy_max_bytes = 0
  )

  if (.Machine$sizeof.pointer == 4L) {
    expect_error(
      do.call(flotsam:::ltsa_assembly_memory_estimates, args),
      "memory estimate overflowed"
    )
  } else {
    estimates <- do.call(flotsam:::ltsa_assembly_memory_estimates, args)
    expect_equal(estimates$serial$raw_entries, 1568750000)
    expect_equal(estimates$serial$full_compact_slots_bound, 2500000000)
    expect_gt(
      estimates$serial$full_compact_slots_bound,
      .Machine$integer.max
    )
    expect_equal(
      estimates$serial$final_sparse_slots_bound,
      .Machine$integer.max
    )
    expect_equal(
      estimates$serial$components_bytes$full_compact_staging_bound,
      2500000000 * estimates$sizeof_bytes$compact_entry
    )
  }
})

test_that("assembly estimator models the protected C++ to R return boundary", {
  estimates <- flotsam:::ltsa_assembly_memory_estimates(
    n = 37,
    k = 6,
    p = 4,
    ndim = 2L,
    requested_assembly_threads = 3L,
    include_self = FALSE,
    row_major_copy_max_bytes = 0
  )

  for (route_name in c("serial", "parallel")) {
    route <- estimates[[route_name]]
    components <- route$components_bytes
    expect_equal(
      components$cpp_to_r_output_copy_bound,
      components$final_sparse_output_bound
    )
    expect_equal(
      components$cpp_to_r_local_rank_copy,
      components$local_rank_staging
    )

    route_staging <- if (route_name == "serial") {
      "compact_column_containers"
    } else {
      c(
        "sparse_slot_offsets",
        "raw_row_staging",
        "raw_value_staging",
        "column_counters",
        "column_starts"
      )
    }
    return_components <- c(
      "neighborhood_index_staging",
      route_staging,
      "local_rank_staging",
      "worker_workspaces",
      "optional_row_major_copy",
      "final_sparse_output_bound",
      "cpp_to_r_output_copy_bound",
      "cpp_to_r_local_rank_copy",
      "control_objects"
    )
    expect_equal(
      route$phase_bytes_bound$return_copy,
      sum(unlist(components[return_components], use.names = FALSE))
    )
  }
})

test_that("assembly estimator accounts for optional row-major and neighbor copies", {
  with_copy <- flotsam:::ltsa_assembly_memory_estimates(
    n = 10,
    k = 4,
    p = 8,
    ndim = 2L,
    requested_assembly_threads = 3L,
    include_self = FALSE,
    row_major_copy_max_bytes = 256 * 1024^2
  )
  without_copy <- flotsam:::ltsa_assembly_memory_estimates(
    n = 10,
    k = 4,
    p = 8,
    ndim = 2L,
    requested_assembly_threads = 3L,
    include_self = FALSE,
    row_major_copy_max_bytes = 0
  )

  expect_true(with_copy$parallel$row_major_copy_included)
  expect_false(without_copy$parallel$row_major_copy_included)
  expect_equal(
    with_copy$parallel$components_bytes$optional_row_major_copy,
    10 * 8 * with_copy$sizeof_bytes$double
  )
  expect_equal(
    with_copy$parallel$components_bytes$neighborhood_index_staging,
    2 * 10 * 4 * with_copy$sizeof_bytes$int
  )
})

test_that("assembly estimator rejects overflow without large allocations", {
  expect_error(
    flotsam:::ltsa_assembly_memory_estimates(
      n = .Machine$integer.max - 1,
      k = .Machine$integer.max - 1,
      p = 1,
      ndim = 1L,
      requested_assembly_threads = 2L,
      include_self = TRUE,
      row_major_copy_max_bytes = 0
    ),
    "Too many LTSA triplets"
  )
  overflow_args <- if (.Machine$sizeof.pointer == 4L) {
    list(n = 2, k = 2, p = .Machine$integer.max)
  } else {
    list(n = 1, k = 2^30, p = 2^30)
  }
  expect_error(
    do.call(
      flotsam:::ltsa_assembly_memory_estimates,
      c(
        overflow_args,
        list(
          ndim = 1L,
          requested_assembly_threads = 2L,
          include_self = TRUE,
          row_major_copy_max_bytes = 0
        )
      )
    ),
    "memory estimate overflowed"
  )
})

test_that("public assembly cap falls back to serial with intact diagnostics", {
  set.seed(102)
  X <- matrix(rnorm(18L * 14L), nrow = 18L)
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)
  estimates <- flotsam:::ltsa_assembly_memory_estimates(
    n = nrow(X),
    k = ncol(nn_idx),
    p = ncol(X),
    ndim = 2L,
    requested_assembly_threads = 3L,
    include_self = TRUE,
    row_major_copy_max_bytes = 256 * 1024^2
  )
  fallback_cap_mib <- mean(c(
    estimates$serial$modeled_peak_mib_bound,
    estimates$parallel$modeled_peak_mib_bound
  ))

  serial <- ltsa(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    eig_method = "eig",
    eig_k = 4L,
    output = "result",
    include_B = TRUE,
    n_assembly_threads = 1L
  )
  expect_warning(
    fallback <- ltsa(
      X,
      ndim = 2L,
      nn_method = nn_idx,
      eig_method = "eig",
      eig_k = 4L,
      output = "result",
      include_B = TRUE,
      n_assembly_threads = 3L,
      assembly_max_mib = fallback_cap_mib
    ),
    "using serial assembly"
  )

  expect_sparse_equivalent(fallback$B, serial$B, tolerance = 0)
  expect_identical(
    fallback$assembly$local_rank_histogram,
    serial$assembly$local_rank_histogram
  )
  expect_identical(
    fallback$assembly$rank_deficient_neighborhood_indices,
    serial$assembly$rank_deficient_neighborhood_indices
  )
  expect_identical(fallback$assembly$assembly_route, "serial_triangular")
  expect_identical(fallback$assembly$requested_assembly_threads, 3L)
  expect_identical(fallback$assembly$effective_assembly_threads, 1L)
  expect_identical(fallback$assembly$memory$requested_route, "parallel")
  expect_identical(fallback$assembly$memory$selected_route, "serial")
  expect_identical(
    fallback$assembly$parallel_fallback_reason,
    "parallel_estimate_exceeds_assembly_cap"
  )
  expect_false("raw_bytes_estimate" %in% names(fallback$assembly))
})

test_that("public assembly cap fails before assembly staging", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  expect_error(
    ltsa(
      X,
      ndim = 2L,
      nn_method = nn_idx,
      output = "B",
      n_assembly_threads = 3L,
      assembly_max_mib = 0
    ),
    "memory preflight failed before assembly staging.*Largest modeled components"
  )
})

test_that("parallel assembly matches serial on low-p SVD route", {
  set.seed(11)
  X <- matrix(rnorm(22L * 3L), nrow = 22L)
  nn_idx <- exact_nn_idx(X, n_neighbors = 7L, include_self = TRUE)

  expect_ltsa_assembly_parallel_matches(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 4L
  )
})

test_that("parallel assembly preserves rank-deficient Gram and SVD metadata", {
  gram_X <- outer(seq_len(18L), seq_len(12L))
  gram_nn <- exact_nn_idx(gram_X, n_neighbors = 5L, include_self = TRUE)
  gram_parallel <- expect_ltsa_assembly_parallel_matches(
    gram_X,
    gram_nn,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 3L
  )
  expect_identical(gram_parallel$rank_deficient_count, nrow(gram_X))
  expect_identical(gram_parallel$min_local_rank, 1L)

  x <- seq_len(16L)
  svd_X <- cbind(x, 2 * x, -3 * x)
  svd_nn <- exact_nn_idx(svd_X, n_neighbors = 6L, include_self = TRUE)
  svd_parallel <- expect_ltsa_assembly_parallel_matches(
    svd_X,
    svd_nn,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 3L
  )
  expect_identical(svd_parallel$rank_deficient_count, nrow(svd_X))
  expect_identical(svd_parallel$min_local_rank, 1L)
})

test_that("public results expose serial and parallel local-rank diagnostics", {
  x <- seq_len(18L)
  fixtures <- list(
    svd = cbind(x, 2 * x, -3 * x),
    gram = outer(x, seq_len(12L))
  )
  neighborhood_sizes <- c(svd = 6L, gram = 5L)

  for (route in names(fixtures)) {
    X <- fixtures[[route]]
    nn_idx <- exact_nn_idx(
      X,
      n_neighbors = neighborhood_sizes[[route]],
      include_self = TRUE
    )
    expect_warning(
      serial <- ltsa(
        X,
        ndim = 2L,
        nn_method = nn_idx,
        include_self = TRUE,
        eig_method = "eig",
        eig_k = 4L,
        output = "result",
        n_assembly_threads = 1L
      ),
      "numerical rank"
    )
    expect_warning(
      parallel <- ltsa(
        X,
        ndim = 2L,
        nn_method = nn_idx,
        include_self = TRUE,
        eig_method = "eig",
        eig_k = 4L,
        output = "result",
        n_assembly_threads = 3L
      ),
      "numerical rank"
    )

    max_rank <- min(ncol(X), ncol(nn_idx))
    expected_histogram <- stats::setNames(
      integer(max_rank + 1L),
      0:max_rank
    )
    expected_histogram[["1"]] <- nrow(X)
    expect_identical(serial$assembly$local_solver_route, route)
    expect_identical(parallel$assembly$local_solver_route, route)
    expect_identical(
      serial$assembly$local_rank_histogram,
      expected_histogram
    )
    expect_identical(
      parallel$assembly$local_rank_histogram,
      expected_histogram
    )
    expect_identical(
      serial$assembly$rank_deficient_neighborhood_indices,
      seq_len(nrow(X))
    )
    expect_identical(
      parallel$assembly$rank_deficient_neighborhood_indices,
      seq_len(nrow(X))
    )
  }
})

test_that("parallel assembly keeps centered Gram stable on hostile offsets", {
  set.seed(1)
  n <- 200L
  p <- 100L
  X <- matrix(rnorm(n * p, sd = 1e-3), n, p) + 1e6
  nn_idx <- exact_nn_idx(X, n_neighbors = 8L, include_self = TRUE)

  expect_ltsa_assembly_parallel_matches(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 4L,
    tolerance = 1e-5
  )
})

test_that("parallel assembly handles more requested threads than observations", {
  set.seed(13)
  X <- matrix(rnorm(5L * 6L), nrow = 5L)
  nn_idx <- exact_nn_idx(X, n_neighbors = 4L, include_self = TRUE)

  parallel <- expect_ltsa_assembly_parallel_matches(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 8L
  )
  expect_identical(parallel$diagnostics$requested_assembly_threads, 8L)
  expect_identical(parallel$diagnostics$effective_assembly_threads, 5L)
})

test_that("ltsa B and embedding paths agree between serial and parallel assembly", {
  set.seed(14)
  X <- matrix(rnorm(16L * 8L), nrow = 16L)

  serial_B <- ltsa(
    X,
    n_neighbors = 5L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_assembly_threads = 1L
  )
  parallel_B <- ltsa(
    X,
    n_neighbors = 5L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_assembly_threads = 4L
  )
  expect_sparse_equivalent(parallel_B, serial_B, tolerance = 1e-11)

  for (normalize in c(FALSE, TRUE)) {
    serial_embedding <- ltsa(
      X,
      n_neighbors = 5L,
      ndim = 2L,
      nn_method = "exact",
      eig_method = "eig",
      include_self = TRUE,
      normalize = normalize,
      n_assembly_threads = 1L
    )
    parallel_embedding <- ltsa(
      X,
      n_neighbors = 5L,
      ndim = 2L,
      nn_method = "exact",
      eig_method = "eig",
      include_self = TRUE,
      normalize = normalize,
      n_assembly_threads = 4L
    )
    expect_embedding_subspace_equivalent(
      parallel_embedding,
      serial_embedding,
      tolerance = 1e-8
    )
  }
})

test_that("serial assembly matches R triplet reference", {
  X <- as.matrix(iris[seq_len(20L), seq_len(4L)])

  for (include_self in c(TRUE, FALSE)) {
    nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = include_self)
    reference <- assemble_ltsa_B_r_triplet_reference(
      X,
      nn_idx,
      ndim = 2L,
      include_self = include_self
    )
    candidate <- flotsam:::assemble_ltsa_B(
      X,
      nn_idx,
      ndim = 2L,
      include_self = include_self
    )

    expect_sparse_equivalent(candidate$B, reference$B)
    expect_identical(
      candidate$rank_deficient_count,
      reference$rank_deficient_count
    )
    expect_identical(candidate$min_local_rank, reference$min_local_rank)
    expect_true(Matrix::isSymmetric(candidate$B))
    expect_equal(sum(candidate$B@x == 0), 0)
  }
})

test_that("serial assembly Gram path matches R triplet reference", {
  set.seed(2)
  X <- matrix(rnorm(18L * 12L), nrow = 18L)
  nn_idx <- exact_nn_idx(X, n_neighbors = 5L, include_self = TRUE)
  expect_gt(ncol(X), ncol(nn_idx))

  reference <- assemble_ltsa_B_r_triplet_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )

  expect_sparse_equivalent(candidate$B, reference$B, tolerance = 1e-11)
  expect_identical(
    candidate$rank_deficient_count,
    reference$rank_deficient_count
  )
  expect_identical(candidate$min_local_rank, reference$min_local_rank)
  expect_true(Matrix::isSymmetric(candidate$B))
  expect_equal(sum(candidate$B@x == 0), 0)
})

test_that("triangular production assembly matches R reference on shuffled neighborhoods", {
  set.seed(3)
  X <- matrix(rnorm(8L * 10L), nrow = 8L)
  # fmt: skip
  nn_idx <- matrix(
    c(
      4L, 1L, 7L, 2L,
      6L, 2L, 8L, 3L,
      1L, 5L, 3L, 8L,
      7L, 4L, 2L, 6L,
      2L, 8L, 5L, 1L,
      3L, 6L, 4L, 7L,
      8L, 7L, 1L, 5L,
      5L, 3L, 6L, 4L
    ),
    nrow = 8L,
    byrow = TRUE
  )
  storage.mode(nn_idx) <- "integer"
  expect_gt(ncol(X), ncol(nn_idx))

  reference <- assemble_ltsa_B_r_triplet_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )

  expect_sparse_equivalent(candidate$B, reference$B, tolerance = 1e-11)
  expect_identical(
    candidate$rank_deficient_count,
    reference$rank_deficient_count
  )
  expect_identical(candidate$min_local_rank, reference$min_local_rank)
})

test_that("serial assembly Gram path preserves rank deficiency", {
  X <- outer(seq_len(18L), seq_len(12L))
  nn_idx <- exact_nn_idx(X, n_neighbors = 5L, include_self = TRUE)
  expect_gt(ncol(X), ncol(nn_idx))

  reference <- assemble_ltsa_B_r_triplet_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )

  expect_sparse_equivalent(candidate$B, reference$B, tolerance = 1e-11)
  expect_identical(
    candidate$rank_deficient_count,
    reference$rank_deficient_count
  )
  expect_identical(candidate$min_local_rank, reference$min_local_rank)
  expect_identical(candidate$rank_deficient_count, nrow(X))
  expect_identical(candidate$min_local_rank, 1L)
})

test_that("triangular assembly preserves low-p rank deficiency", {
  x <- seq_len(16L)
  X <- cbind(x, 2 * x, -3 * x)
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  reference <- assemble_ltsa_B_r_triplet_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )

  expect_sparse_equivalent(candidate$B, reference$B, tolerance = 1e-11)
  expect_identical(
    candidate$rank_deficient_count,
    reference$rank_deficient_count
  )
  expect_identical(candidate$min_local_rank, reference$min_local_rank)
  expect_identical(candidate$rank_deficient_count, nrow(X))
  expect_identical(candidate$min_local_rank, 1L)
})

test_that("centered high-p Gram path is stable on hostile large-offset data", {
  set.seed(1)
  n <- 200L
  p <- 100L
  X <- matrix(rnorm(n * p, sd = 1e-3), n, p) + 1e6
  nn_idx <- exact_nn_idx(X, n_neighbors = 8L, include_self = TRUE)
  expect_gt(ncol(X), ncol(nn_idx))

  reference <- assemble_ltsa_B_r_triplet_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )

  expect_sparse_equivalent(candidate$B, reference$B, tolerance = 1e-5)
  expect_identical(
    candidate$rank_deficient_count,
    reference$rank_deficient_count
  )
  expect_identical(candidate$min_local_rank, reference$min_local_rank)
})

test_that("Gram assembly handles near-machine SVD-rank boundary", {
  k <- 6L
  p <- 12L
  X <- matrix(0, nrow = k, ncol = p)
  X[, 1L] <- c(1, -1, 0, 0, 0, 0) / sqrt(2)
  # The second singular direction is above the direct-SVD tolerance but below
  # the Gram eigenvalue tolerance used to avoid promoting Gram roundoff.
  X[, 2L] <- 1e-8 * c(0, 0, 1, -1, 0, 0) / sqrt(2)
  nn_idx <- matrix(rep(seq_len(k), times = k), nrow = k, byrow = TRUE)
  storage.mode(nn_idx) <- "integer"

  reference <- assemble_ltsa_B_r_triplet_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )

  expect_gt(ncol(X), ncol(nn_idx))
  expect_identical(reference$rank_deficient_count, 0L)
  expect_identical(candidate$rank_deficient_count, k)
  expect_identical(candidate$min_local_rank, 1L)
  expect_gt(max(abs(as.matrix(candidate$B - reference$B))), 0.1)
})

test_that("ltsa output B uses production assembly", {
  X <- as.matrix(iris[seq_len(15L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 5L, include_self = FALSE)
  reference <- assemble_ltsa_B_r_triplet_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = FALSE
  )

  candidate <- ltsa(
    X,
    n_neighbors = 5L,
    ndim = 2L,
    nn_method = "exact",
    include_self = FALSE,
    output = "B"
  )

  expect_sparse_equivalent(candidate, reference$B)
  expect_equal(sum(candidate@x == 0), 0)
})

test_that("output B returns the unchanged full dgCMatrix class", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])

  candidate <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B"
  )

  expect_s4_class(candidate, "dgCMatrix")
  expect_false(methods::is(candidate, "dsCMatrix"))
  expect_true(Matrix::isSymmetric(candidate))
  expect_equal(sum(candidate@x == 0), 0)
})
