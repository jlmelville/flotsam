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
