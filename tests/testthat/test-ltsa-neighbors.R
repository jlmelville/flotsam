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

capture_ltsa_messages <- function(...) {
  capture.output(invisible(ltsa(...)), type = "message")
}

# fmt: skip
duplicate_nn_idx <- function() {
  matrix(
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
}

test_that("precomputed exact neighborhoods match computed exact LTSA B", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])

  for (include_self in c(TRUE, FALSE)) {
    nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = include_self)
    computed <- ltsa(
      X,
      n_neighbors = 6L,
      ndim = 2L,
      nn_method = "exact",
      include_self = include_self,
      output = "B"
    )
    precomputed <- ltsa(
      X,
      n_neighbors = NULL,
      ndim = 2L,
      nn_method = nn_idx,
      include_self = include_self,
      output = "B"
    )

    expect_sparse_equivalent(precomputed, computed, tolerance = 0)
  }
})

test_that("precomputed graph supplied as nn_method skips nearest-neighbor search", {
  set.seed(20)
  X <- matrix(rnorm(8L * 10L), nrow = 8L)
  nn_idx <- duplicate_nn_idx()
  storage.mode(nn_idx) <- "integer"

  reference <- flotsam:::assemble_ltsa_B(
    X = X,
    nn_idx = nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B"
  )

  expect_sparse_equivalent(candidate, reference$B, tolerance = 1e-11)
})

test_that("nn_method can carry a precomputed neighbor graph", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)
  reference <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B"
  )

  from_matrix <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B"
  )
  from_result <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = list(idx = nn_idx),
    include_self = TRUE,
    output = "B"
  )

  expect_sparse_equivalent(from_matrix, reference, tolerance = 0)
  expect_sparse_equivalent(from_result, reference, tolerance = 0)
})

test_that("precomputed duplicate neighborhoods work with serial and parallel assembly", {
  set.seed(21)
  X <- matrix(rnorm(8L * 10L), nrow = 8L)
  nn_idx <- duplicate_nn_idx()
  storage.mode(nn_idx) <- "integer"

  serial <- ltsa(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B",
    n_assembly_threads = 1L
  )
  parallel <- ltsa(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B",
    n_assembly_threads = 3L
  )

  expect_sparse_equivalent(parallel, serial, tolerance = 1e-11)
})

test_that("precomputed neighbor graph validation rejects invalid graphs", {
  X <- as.matrix(iris[seq_len(8L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 4L, include_self = TRUE)

  expect_error(
    ltsa(X, ndim = 2L, nn_method = nn_idx[-1L, ], output = "B"),
    "one row per observation"
  )

  bad <- nn_idx
  bad[1L, 2L] <- 9L
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, output = "B"),
    "between 1 and nrow"
  )

  bad <- nn_idx
  bad[1L, 2L] <- NA_integer_
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, output = "B"),
    "finite whole-number"
  )

  bad <- nn_idx
  storage.mode(bad) <- "double"
  bad[1L, 2L] <- bad[1L, 2L] + 0.5
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, output = "B"),
    "whole-number"
  )

  expect_error(
    ltsa(X, ndim = 2L, nn_method = as.vector(nn_idx), output = "B"),
    "nn_method"
  )

  expect_error(
    ltsa(X, ndim = 2L, nn_method = nn_idx[, 1:2], output = "B"),
    "greater than ndim"
  )

  expect_error(
    ltsa(
      X,
      n_neighbors = 5L,
      ndim = 2L,
      nn_method = nn_idx,
      output = "B"
    ),
    "ncol\\(nn_method\\)"
  )

  bad <- nn_idx
  bad[1L, ] <- c(2L, 3L, 4L, 5L)
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, include_self = TRUE, output = "B"),
    "own row index"
  )

  nn_idx_no_self <- exact_nn_idx(X, n_neighbors = 4L, include_self = FALSE)
  bad <- nn_idx_no_self
  bad[2L, 1L] <- 1L
  expect_error(
    ltsa(
      X,
      ndim = 2L,
      nn_method = bad,
      include_self = FALSE,
      output = "B"
    ),
    "first column"
  )
})

test_that("precomputed graph with n_threads does not warn", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  expect_warning(
    invisible(ltsa(
      X,
      ndim = 2L,
      nn_method = nn_idx,
      include_self = TRUE,
      output = "B",
      n_threads = 2L
    )),
    NA
  )
})

test_that("verbose output describes computed and precomputed neighbor handling", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  computed_messages <- capture_ltsa_messages(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_threads = 2L,
    verbose = TRUE
  )
  expect_true(any(grepl(
    "Finding nearest neighbors with method 'exact' using n_threads = 2",
    computed_messages,
    fixed = TRUE
  )))

  precomputed_messages <- capture_ltsa_messages(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B",
    n_threads = 2L,
    n_assembly_threads = 3L,
    verbose = TRUE
  )
  expect_true(any(grepl(
    "Using precomputed nearest-neighbor graph with k = 6",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Ignoring n_threads = 2 because precomputed nearest-neighbor graph was supplied",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "assembly workers requested/active: 3/",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "raw staging memory estimate: [0-9,.]+ [KMGT]?i?B",
    precomputed_messages
  )))
  expect_match(flotsam:::format_ltsa_count(2790000), "2[^0-9]790[^0-9]000")
})
