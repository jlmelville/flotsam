expect_sparse_identical <- function(candidate, reference) {
  expect_s4_class(candidate, "dgCMatrix")
  expect_s4_class(reference, "dgCMatrix")
  expect_identical(candidate@Dim, reference@Dim)
  expect_identical(candidate@p, reference@p)
  expect_identical(candidate@i, reference@i)
  expect_identical(candidate@x, reference@x)
}

expect_sparse_equivalent <- function(candidate, reference) {
  expect_sparse_identical(Matrix::drop0(candidate), Matrix::drop0(reference))
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

test_that("C++ triplet assembly matches slot-search assembly with self", {
  X <- as.matrix(iris[seq_len(20L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  reference <- flotsam:::assemble_ltsa_B_slot_search_reference(
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

  expect_sparse_equivalent(candidate$B, reference$B)
  expect_identical(candidate$rank_deficient_count, reference$rank_deficient_count)
  expect_identical(candidate$min_local_rank, reference$min_local_rank)
  expect_true(Matrix::isSymmetric(candidate$B))
})

test_that("append/finalize assembly matches compact assembly", {
  X <- as.matrix(iris[seq_len(20L), seq_len(4L)])

  for (include_self in c(TRUE, FALSE)) {
    nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = include_self)

    reference <- flotsam:::assemble_ltsa_B_compact(
      X,
      nn_idx,
      ndim = 2L,
      include_self = include_self
    )
    candidate <- flotsam:::assemble_ltsa_B_append(
      X,
      nn_idx,
      ndim = 2L,
      include_self = include_self
    )

    expect_sparse_identical(candidate$B, reference$B)
    expect_identical(candidate$rank_deficient_count, reference$rank_deficient_count)
    expect_identical(candidate$min_local_rank, reference$min_local_rank)
    expect_true(Matrix::isSymmetric(candidate$B))
  }
})

test_that("C++ triplet assembly drops excluded-self zero-only pattern by default", {
  X <- as.matrix(iris[seq_len(20L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = FALSE)

  reference <- flotsam:::assemble_ltsa_B_slot_search_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = FALSE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = FALSE
  )

  expect_sparse_equivalent(candidate$B, reference$B)
  expect_lt(length(candidate$B@x), length(reference$B@x))
  expect_equal(sum(candidate$B@x == 0), 0)
  expect_true(Matrix::isSymmetric(candidate$B))
})

test_that("C++ triplet assembly can preserve excluded-self sparse pattern", {
  X <- as.matrix(iris[seq_len(20L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = FALSE)

  reference <- flotsam:::assemble_ltsa_B_slot_search_reference(
    X,
    nn_idx,
    ndim = 2L,
    include_self = FALSE
  )
  candidate <- flotsam:::assemble_ltsa_B(
    X,
    nn_idx,
    ndim = 2L,
    include_self = FALSE,
    preserve_pattern = TRUE
  )

  expect_sparse_identical(candidate$B, reference$B)
  expect_gt(sum(candidate$B@x == 0), 0)
  expect_true(Matrix::isSymmetric(candidate$B))
})

test_that("ltsa ret_B uses assembly matching the slot-search reference", {
  X <- as.matrix(iris[seq_len(15L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 5L, include_self = FALSE)
  reference <- flotsam:::assemble_ltsa_B_slot_search_reference(
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
})

test_that("triplet assembly component helper validates dimensions", {
  expect_error(
    flotsam:::ltsa_triplet_assembly_components(
      pattern_nnt = c(1L, 2L, 3L),
      pattern_n_nbrs = 2L,
      value_nnt = c(1L, 2L),
      weights = c(0, 0, 0, 0),
      value_n_nbrs = 2L,
      preserve_pattern = TRUE
    ),
    "Inconsistent pattern neighborhood dimensions"
  )

  expect_error(
    flotsam:::ltsa_triplet_assembly_components(
      pattern_nnt = c(1L, 2L),
      pattern_n_nbrs = 2L,
      value_nnt = c(1L, 2L),
      weights = c(0, 0),
      value_n_nbrs = 2L,
      preserve_pattern = FALSE
    ),
    "Inconsistent local weight dimensions"
  )
})

test_that("triplet assembly component helper validates neighbor indices", {
  expect_error(
    flotsam:::ltsa_triplet_assembly_components(
      pattern_nnt = c(1L, 3L),
      pattern_n_nbrs = 2L,
      value_nnt = c(1L, 2L),
      weights = c(0, 0, 0, 0),
      value_n_nbrs = 2L,
      preserve_pattern = TRUE
    ),
    "outside the sparse matrix dimensions"
  )

  expect_error(
    flotsam:::ltsa_triplet_assembly_components(
      pattern_nnt = c(1L, 2L),
      pattern_n_nbrs = 2L,
      value_nnt = c(1L, 3L),
      weights = c(0, 0, 0, 0),
      value_n_nbrs = 2L,
      preserve_pattern = FALSE
    ),
    "outside the sparse matrix dimensions"
  )
})
