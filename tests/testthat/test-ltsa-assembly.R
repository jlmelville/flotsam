expect_parallel_assembly_matches <- function(
  X,
  nn_idx,
  ndim,
  include_self,
  n_assembly_threads = 3L,
  tolerance = 1e-11
) {
  serial <- flotsam:::assemble_alignment_matrix(
    X,
    nn_idx,
    ndim = ndim,
    include_self = include_self,
    n_assembly_threads = 1L
  )
  parallel <- flotsam:::assemble_alignment_matrix(
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
  removed_fields <- c(
    "memory",
    "raw_entries_estimate",
    "raw_bytes_estimate",
    "row_major_used",
    "row_major_fallback_reason",
    "parallel_fallback_reason"
  )
  expect_false(any(removed_fields %in% names(serial$diagnostics)))
  expect_false(any(removed_fields %in% names(parallel$diagnostics)))
  expect_true(Matrix::isSymmetric(parallel$B))
  expect_equal(sum(parallel$B@x == 0), 0)

  invisible(parallel)
}

test_that("compiled effective components match the pure-R oracle", {
  connected <- t(vapply(
    seq_len(8L),
    function(i) as.integer((i - 1L + 0:2) %% 8L + 1L),
    integer(3L)
  ))
  # Queries are intentionally absent from some rows. Components follow only
  # assembled neighborhood co-membership, including isolated observations.
  # fmt: skip
  isolated <- matrix(
    c(
      1L, 2L,
      1L, 2L,
      1L, 2L,
      4L, 5L,
      4L, 5L,
      4L, 5L
    ),
    nrow = 6L,
    byrow = TRUE
  )
  # fmt: skip
  disconnected <- matrix(
    c(
      3L, 1L, 2L,
      2L, 3L, 1L,
      1L, 2L, 3L,
      6L, 4L, 5L,
      5L, 6L, 4L,
      4L, 5L, 6L
    ),
    nrow = 6L,
    byrow = TRUE
  )

  for (indices in list(connected, isolated, disconnected)) {
    reference <- flotsam:::compute_effective_components(
      indices,
      nrow(indices)
    )
    compiled <- flotsam:::compute_effective_components_cpp(
      t(indices),
      nrow(indices),
      ncol(indices)
    )
    expect_identical(compiled, reference)
  }

  expect_identical(
    flotsam:::compute_effective_components_cpp(t(isolated), 6L, 2L),
    list(
      component_count = 4L,
      component_sizes = c(2L, 1L, 2L, 1L),
      component_membership = c(1L, 1L, 2L, 3L, 3L, 4L)
    )
  )
})

test_that("compiled effective components reject invalid native inputs", {
  expect_error(
    flotsam:::compute_effective_components_cpp(1:3, 2L, 2L),
    "Inconsistent component graph dimensions"
  )
  expect_error(
    flotsam:::compute_effective_components_cpp(c(1L, 3L), 2L, 1L),
    "outside the sparse matrix dimensions"
  )
  expect_error(
    flotsam:::compute_effective_components_cpp(integer(), 0L, 1L),
    "n_obs must be positive"
  )
  expect_error(
    flotsam:::compute_effective_components_cpp(integer(), 1L, 0L),
    "n_neighbors must be positive"
  )
})

test_that("production assembly bypasses the pure-R component oracle", {
  set.seed(20260828)
  X <- matrix(stats::rnorm(12L * 5L), nrow = 12L)
  nn_idx <- t(vapply(
    seq_len(nrow(X)),
    function(i) as.integer((i - 1L + 0:3) %% nrow(X) + 1L),
    integer(4L)
  ))
  local_mocked_bindings(
    compute_effective_components = function(...) {
      stop("pure-R component oracle reached", call. = FALSE)
    },
    .package = "flotsam"
  )

  expect_no_error(flotsam:::assemble_alignment_matrix(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  ))
})

raw_exact_neighbor_indices <- function(X, n_neighbors, include_self) {
  nn <- rnndescent::brute_force_knn(
    data = X,
    k = ifelse(include_self, n_neighbors, n_neighbors + 1L),
    n_threads = 0L,
    verbose = FALSE
  )
  mode(nn$idx) <- "integer"
  nn$idx
}

reference_local_weights <- function(X, neighbor_indices, ndim) {
  Xi <- scale(
    X[neighbor_indices, , drop = FALSE],
    center = TRUE,
    scale = FALSE
  )
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

  k <- length(neighbor_indices)
  Gi <- cbind(1 / sqrt(k), res$u[, keep, drop = FALSE])
  Wi <- -tcrossprod(Gi)
  diag(Wi) <- diag(Wi) + 1.0

  list(Wi = Wi, rank = rank)
}

reference_triplet_assembly <- function(X, nn_idx, ndim, include_self) {
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
    neighbor_indices <- weight_idx[i, ]
    local <- reference_local_weights(X, neighbor_indices, ndim)
    if (local$rank < ndim) {
      rank_deficient_count <- rank_deficient_count + 1L
      min_local_rank <- min(min_local_rank, local$rank)
    }

    idx <- ((i - 1L) * k * k + 1L):(i * k * k)
    rows[idx] <- rep.int(neighbor_indices, times = k)
    cols[idx] <- rep(neighbor_indices, each = k)
    vals[idx] <- as.vector(local$Wi)
  }

  list(
    B = Matrix::sparseMatrix(
      i = rows,
      j = cols,
      x = vals,
      dims = c(n, n),
      repr = "C"
    ),
    rank_deficient_count = rank_deficient_count,
    min_local_rank = min_local_rank
  )
}

test_that("default assembly behavior remains serial", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 6L, include_self = TRUE)

  default <- flotsam:::assemble_alignment_matrix(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  explicit <- flotsam:::assemble_alignment_matrix(
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
  default_B <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B"
  )
  explicit_B <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_assembly_threads = 1L
  )
  expect_sparse_equivalent(default_B, explicit_B, tolerance = 0)
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
    nn_idx <- raw_exact_neighbor_indices(
      X,
      n_neighbors = 6L,
      include_self = include_self
    )
    parallel <- expect_parallel_assembly_matches(
      X,
      nn_idx,
      ndim = 2L,
      include_self = include_self,
      n_assembly_threads = 4L
    )
    expect_identical(
      parallel$diagnostics$assembly_route,
      "parallel_triangular_two_pass"
    )
  }
})

test_that("public assembly thread control selects the requested route", {
  set.seed(102)
  X <- matrix(rnorm(18L * 14L), nrow = 18L)
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 6L, include_self = TRUE)
  common_args <- list(
    X = X,
    ndim = 2L,
    nn_method = nn_idx,
    eig_method = "eig",
    eig_k = 4L,
    output = "result",
    include_B = TRUE
  )

  serial <- expect_silent(
    do.call(ltsa, c(common_args, list(n_assembly_threads = 1L)))
  )
  parallel <- expect_silent(
    do.call(ltsa, c(common_args, list(n_assembly_threads = 3L)))
  )

  assembly_fields <- c(
    "n_neighbors",
    "include_self",
    "neighbor_source",
    "neighbor_elapsed",
    "rank_deficient_count",
    "min_local_rank",
    "assembly_route",
    "component_count",
    "component_sizes",
    "component_membership",
    "requested_assembly_threads",
    "effective_assembly_threads"
  )
  expect_named(serial$assembly, assembly_fields)
  expect_named(parallel$assembly, assembly_fields)
  expect_sparse_equivalent(parallel$B, serial$B, tolerance = 1e-11)
  expect_identical(serial$assembly$assembly_route, "serial_triangular")
  expect_identical(serial$assembly$requested_assembly_threads, 1L)
  expect_identical(serial$assembly$effective_assembly_threads, 1L)
  expect_identical(
    parallel$assembly$assembly_route,
    "parallel_triangular_two_pass"
  )
  expect_identical(parallel$assembly$requested_assembly_threads, 3L)
  expect_gte(parallel$assembly$effective_assembly_threads, 1L)
  expect_lte(parallel$assembly$effective_assembly_threads, 3L)

  removed_fields <- c(
    "memory",
    "raw_entries_estimate",
    "raw_bytes_estimate",
    "row_major_used",
    "row_major_fallback_reason",
    "parallel_fallback_reason"
  )
  expect_false(any(removed_fields %in% names(serial$assembly)))
  expect_false(any(removed_fields %in% names(parallel$assembly)))
})

test_that("parallel assembly matches serial on low-p SVD route", {
  set.seed(11)
  X <- matrix(rnorm(22L * 3L), nrow = 22L)
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 7L, include_self = TRUE)

  expect_parallel_assembly_matches(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 4L
  )
})

test_that("parallel assembly preserves rank-deficiency summaries", {
  gram_X <- outer(seq_len(18L), seq_len(12L))
  gram_nn <- raw_exact_neighbor_indices(
    gram_X,
    n_neighbors = 5L,
    include_self = TRUE
  )
  gram_parallel <- expect_parallel_assembly_matches(
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
  svd_nn <- raw_exact_neighbor_indices(
    svd_X,
    n_neighbors = 6L,
    include_self = TRUE
  )
  svd_parallel <- expect_parallel_assembly_matches(
    svd_X,
    svd_nn,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 3L
  )
  expect_identical(svd_parallel$rank_deficient_count, nrow(svd_X))
  expect_identical(svd_parallel$min_local_rank, 1L)
})

test_that("public results expose serial and parallel rank summaries", {
  x <- seq_len(18L)
  fixtures <- list(
    svd = cbind(x, 2 * x, -3 * x),
    gram = outer(x, seq_len(12L))
  )
  neighborhood_sizes <- c(svd = 6L, gram = 5L)

  for (route in names(fixtures)) {
    X <- fixtures[[route]]
    nn_idx <- raw_exact_neighbor_indices(
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

    expect_identical(serial$assembly$rank_deficient_count, nrow(X))
    expect_identical(parallel$assembly$rank_deficient_count, nrow(X))
    expect_identical(serial$assembly$min_local_rank, 1L)
    expect_identical(parallel$assembly$min_local_rank, 1L)
  }
})

test_that("parallel assembly keeps centered Gram stable on hostile offsets", {
  set.seed(1)
  n <- 200L
  p <- 100L
  X <- matrix(rnorm(n * p, sd = 1e-3), n, p) + 1e6
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 8L, include_self = TRUE)

  expect_parallel_assembly_matches(
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
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 4L, include_self = TRUE)

  parallel <- expect_parallel_assembly_matches(
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
    expect_same_subspace(
      parallel_embedding,
      serial_embedding,
      tolerance = 1e-8
    )
  }
})

test_that("serial assembly matches R triplet reference", {
  X <- as.matrix(iris[seq_len(20L), seq_len(4L)])

  for (include_self in c(TRUE, FALSE)) {
    nn_idx <- raw_exact_neighbor_indices(
      X,
      n_neighbors = 6L,
      include_self = include_self
    )
    reference <- reference_triplet_assembly(
      X,
      nn_idx,
      ndim = 2L,
      include_self = include_self
    )
    candidate <- flotsam:::assemble_alignment_matrix(
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
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 5L, include_self = TRUE)
  expect_gt(ncol(X), ncol(nn_idx))

  reference <- reference_triplet_assembly(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_alignment_matrix(
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

  reference <- reference_triplet_assembly(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_alignment_matrix(
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
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 5L, include_self = TRUE)
  expect_gt(ncol(X), ncol(nn_idx))

  reference <- reference_triplet_assembly(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_alignment_matrix(
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
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 6L, include_self = TRUE)

  reference <- reference_triplet_assembly(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_alignment_matrix(
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
  nn_idx <- raw_exact_neighbor_indices(X, n_neighbors = 8L, include_self = TRUE)
  expect_gt(ncol(X), ncol(nn_idx))

  reference <- reference_triplet_assembly(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_alignment_matrix(
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

  reference <- reference_triplet_assembly(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- flotsam:::assemble_alignment_matrix(
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
  nn_idx <- raw_exact_neighbor_indices(
    X,
    n_neighbors = 5L,
    include_self = FALSE
  )
  reference <- reference_triplet_assembly(
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
