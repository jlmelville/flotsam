expect_sparse_identical <- function(candidate, reference) {
  expect_s4_class(candidate, "dgCMatrix")
  expect_s4_class(reference, "dgCMatrix")
  expect_identical(candidate@Dim, reference@Dim)
  expect_identical(candidate@p, reference@p)
  expect_identical(candidate@i, reference@i)
  expect_identical(candidate@x, reference@x)
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

  expect_sparse_identical(candidate$B, reference$B)
  expect_identical(candidate$rank_deficient_count, reference$rank_deficient_count)
  expect_identical(candidate$min_local_rank, reference$min_local_rank)
  expect_true(Matrix::isSymmetric(candidate$B))
})

test_that("C++ triplet assembly preserves excluded-self sparse pattern", {
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

  expect_sparse_identical(candidate, reference$B)
})

test_that("triplet assembly component helper validates dimensions", {
  expect_error(
    flotsam:::ltsa_triplet_assembly_components(
      pattern_nnt = c(1L, 2L, 3L),
      pattern_n_nbrs = 2L,
      value_nnt = c(1L, 2L),
      weights = c(0, 0, 0, 0),
      value_n_nbrs = 2L
    ),
    "Inconsistent pattern neighborhood dimensions"
  )

  expect_error(
    flotsam:::ltsa_triplet_assembly_components(
      pattern_nnt = c(1L, 2L),
      pattern_n_nbrs = 2L,
      value_nnt = c(1L, 2L),
      weights = c(0, 0),
      value_n_nbrs = 2L
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
      value_n_nbrs = 2L
    ),
    "outside the sparse matrix dimensions"
  )
})
