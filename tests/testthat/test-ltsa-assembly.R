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
  expect_identical(serial$diagnostics$assembly_route, "serial_triangular")
  expect_identical(
    parallel$diagnostics$assembly_route,
    "parallel_triangular_two_pass"
  )
  expect_identical(
    parallel$diagnostics$requested_assembly_threads,
    as.integer(n_assembly_threads)
  )
  expect_true(Matrix::isSymmetric(parallel$B))
  expect_equal(sum(parallel$B@x == 0), 0)

  invisible(parallel)
}

expect_embedding_subspace_equivalent <- function(
  candidate,
  reference,
  tolerance = 1e-8
) {
  candidate_q <- qr.Q(qr(candidate))
  reference_q <- qr.Q(qr(reference))
  expect_equal(
    tcrossprod(candidate_q),
    tcrossprod(reference_q),
    tolerance = tolerance
  )
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
  weight_idx <- flotsam:::ltsa_weight_neighborhoods(nn_idx, include_self)
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
    ret_B = TRUE
  )
  explicit_ltsa <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    ret_B = TRUE,
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
    ret_B = TRUE,
    n_threads = 0L,
    n_assembly_threads = 1L
  )
  serial_nn2 <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    ret_B = TRUE,
    n_threads = 2L,
    n_assembly_threads = 1L
  )

  expect_sparse_equivalent(serial_nn2, serial_nn0, tolerance = 0)
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

test_that("parallel assembly preserves duplicate-neighborhood fallback", {
  set.seed(12)
  X <- matrix(rnorm(8L * 10L), nrow = 8L)
  # fmt: skip
  nn_idx <- matrix(
    c(
      1L, 1L, 2L, 3L,
      2L, 3L, 3L, 4L,
      3L, 4L, 5L, 5L,
      4L, 6L, 6L, 7L,
      5L, 7L, 8L, 8L,
      6L, 1L, 2L, 2L,
      7L, 3L, 4L, 4L,
      8L, 5L, 6L, 6L
    ),
    nrow = 8L,
    byrow = TRUE
  )
  storage.mode(nn_idx) <- "integer"

  parallel <- expect_ltsa_assembly_parallel_matches(
    X,
    nn_idx,
    ndim = 2L,
    include_self = TRUE,
    n_assembly_threads = 4L,
    tolerance = 1e-11
  )
  expect_identical(parallel$diagnostics$duplicate_fallback_count, nrow(X))
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

test_that("ltsa ret_B and embedding paths agree between serial and parallel assembly", {
  set.seed(14)
  X <- matrix(rnorm(16L * 8L), nrow = 16L)

  serial_B <- ltsa(
    X,
    n_neighbors = 5L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    ret_B = TRUE,
    n_assembly_threads = 1L
  )
  parallel_B <- ltsa(
    X,
    n_neighbors = 5L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    ret_B = TRUE,
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
      ret_B = FALSE,
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
      ret_B = FALSE,
      n_assembly_threads = 4L
    )
    expect_embedding_subspace_equivalent(
      parallel_embedding,
      serial_embedding,
      tolerance = 1e-8
    )
  }
})

test_that("triangular builder expands symmetric local weights to full CSC", {
  # fmt: skip
  neighborhoods <- matrix(
    c(
      3L, 1L, 4L, 2L,
      2L, 4L, 1L, 3L,
      4L, 3L, 2L, 1L,
      1L, 2L, 3L, 4L
    ),
    nrow = 4L,
    byrow = TRUE
  )
  k <- ncol(neighborhoods)
  W <- matrix(seq_len(k * k), k, k)
  W <- W + t(W)
  diag(W) <- seq_len(k)
  storage.mode(W) <- "double"

  builder <- flotsam:::ltsa_triplet_builder_create(
    value_nnt = as.integer(t(neighborhoods)),
    value_n_nbrs = k
  )
  rows <- integer(nrow(neighborhoods) * k * k)
  cols <- integer(nrow(neighborhoods) * k * k)
  vals <- numeric(nrow(neighborhoods) * k * k)

  for (obs in seq_len(nrow(neighborhoods))) {
    nni <- neighborhoods[obs, ]
    flotsam:::ltsa_triplet_builder_append(
      builder,
      nni = nni,
      weights = as.vector(W)
    )

    idx <- ((obs - 1L) * k * k + 1L):(obs * k * k)
    rows[idx] <- rep.int(nni, times = k)
    cols[idx] <- rep(nni, each = k)
    vals[idx] <- as.vector(W)
  }

  candidate <- flotsam:::ltsa_components_to_dgCMatrix(
    flotsam:::ltsa_triplet_builder_finalize(builder),
    n = nrow(neighborhoods)
  )
  reference <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = rep(nrow(neighborhoods), 2L),
    giveCsparse = TRUE
  )

  expect_sparse_equivalent(candidate, reference, tolerance = 0)
})

test_that("append/finalize assembly matches R triplet reference", {
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

test_that("append/finalize assembly Gram path matches R triplet reference", {
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

test_that("append/finalize assembly Gram path preserves rank deficiency", {
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

test_that("C++ local weights match exact R SVD reference", {
  set.seed(1)
  low_p <- matrix(rnorm(10 * 3), nrow = 10)
  low_nni <- c(1L, 3L, 4L, 6L, 8L, 10L)
  low_reference <- ltsa_local_weights_r_reference(low_p, low_nni, ndim = 2L)
  low_candidate <- flotsam:::ltsa_local_weights(low_p, low_nni, ndim = 2L)
  expect_equal(low_candidate$Wi, low_reference$Wi, tolerance = 1e-12)
  expect_identical(low_candidate$rank, low_reference$rank)

  high_p <- matrix(rnorm(9 * 12), nrow = 9)
  high_nni <- c(1L, 2L, 4L, 7L, 9L)
  high_reference <- ltsa_local_weights_r_reference(high_p, high_nni, ndim = 2L)
  high_candidate <- flotsam:::ltsa_local_weights(high_p, high_nni, ndim = 2L)
  expect_equal(high_candidate$Wi, high_reference$Wi, tolerance = 1e-12)
  expect_identical(high_candidate$rank, high_reference$rank)
})

test_that("ltsa ret_B uses append/finalize assembly", {
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
    ret_B = TRUE
  )

  expect_sparse_equivalent(candidate, reference$B)
  expect_equal(sum(candidate@x == 0), 0)
})

test_that("ret_B returns the unchanged full dgCMatrix class", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])

  candidate <- ltsa(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    ret_B = TRUE
  )

  expect_s4_class(candidate, "dgCMatrix")
  expect_false(methods::is(candidate, "dsCMatrix"))
  expect_true(Matrix::isSymmetric(candidate))
  expect_equal(sum(candidate@x == 0), 0)
})

test_that("append/finalize builder validates lifecycle", {
  builder <- flotsam:::ltsa_triplet_builder_create(
    value_nnt = c(1L, 2L, 1L, 2L),
    value_n_nbrs = 2L
  )

  expect_error(
    flotsam:::ltsa_triplet_builder_finalize(builder),
    "Not all LTSA neighborhoods were appended"
  )

  flotsam:::ltsa_triplet_builder_append(
    builder,
    nni = c(1L, 2L),
    weights = c(1, 0, 0, 1)
  )
  flotsam:::ltsa_triplet_builder_append(
    builder,
    nni = c(1L, 2L),
    weights = c(1, 0, 0, 1)
  )

  components <- flotsam:::ltsa_triplet_builder_finalize(builder)
  expect_named(components, c("i", "p", "x"))
  expect_error(
    flotsam:::ltsa_triplet_builder_finalize(builder),
    "already been finalized"
  )
  expect_error(
    flotsam:::ltsa_triplet_builder_append(
      builder,
      nni = c(1L, 2L),
      weights = c(1, 0, 0, 1)
    ),
    "already been finalized"
  )
})

test_that("append/finalize builder validates append dimensions and indices", {
  expect_error(
    flotsam:::ltsa_triplet_builder_create(
      value_nnt = integer(),
      value_n_nbrs = 2L
    ),
    "must not be empty"
  )
  expect_error(
    flotsam:::ltsa_triplet_builder_create(
      value_nnt = c(1L, 2L, 1L),
      value_n_nbrs = 2L
    ),
    "Inconsistent value neighborhood dimensions"
  )
  expect_error(
    flotsam:::ltsa_triplet_builder_create(
      value_nnt = c(1L, 3L, 1L, 2L),
      value_n_nbrs = 2L
    ),
    "outside the sparse matrix dimensions"
  )

  builder <- flotsam:::ltsa_triplet_builder_create(
    value_nnt = c(1L, 2L, 1L, 2L),
    value_n_nbrs = 2L
  )
  expect_error(
    flotsam:::ltsa_triplet_builder_append(
      builder,
      nni = c(1L),
      weights = c(1, 0, 0, 1)
    ),
    "Inconsistent value neighborhood dimensions"
  )
  expect_error(
    flotsam:::ltsa_triplet_builder_append(
      builder,
      nni = c(1L, 2L),
      weights = c(1, 0)
    ),
    "Inconsistent local weight dimensions"
  )
  expect_error(
    flotsam:::ltsa_triplet_builder_append(
      builder,
      nni = c(1L, 3L),
      weights = c(1, 0, 0, 1)
    ),
    "outside the sparse matrix dimensions"
  )
})

test_that("append/finalize builder rejects too many appended neighborhoods", {
  builder <- flotsam:::ltsa_triplet_builder_create(
    value_nnt = c(1L, 2L, 1L, 2L),
    value_n_nbrs = 2L
  )

  flotsam:::ltsa_triplet_builder_append(
    builder,
    nni = c(1L, 2L),
    weights = c(1, 0, 0, 1)
  )
  flotsam:::ltsa_triplet_builder_append(
    builder,
    nni = c(1L, 2L),
    weights = c(1, 0, 0, 1)
  )

  expect_error(
    flotsam:::ltsa_triplet_builder_append(
      builder,
      nni = c(1L, 2L),
      weights = c(1, 0, 0, 1)
    ),
    "Too many LTSA neighborhoods appended"
  )
})
