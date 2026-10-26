make_sparse_neighbor_fixture <- function(offsets = seq_len(4L)) {
  n <- 12L
  neighbors <- t(vapply(
    seq_len(n),
    function(i) as.integer((i - 1L + offsets) %% n + 1L),
    integer(length(offsets))
  ))
  set.seed(26)
  list(
    X = matrix(rnorm(n * 5L), nrow = n),
    idx = cbind(seq_len(n), neighbors),
    graph = Matrix::sparseMatrix(
      i = rep(seq_len(n), each = length(offsets)),
      j = as.vector(t(neighbors)),
      x = rep(c(-3, 0.2, 4, 100), length.out = length(neighbors)),
      dims = c(n, n)
    )
  )
}

test_that("sparse neighbor graphs reproduce assembly from index matrices", {
  fixture <- make_sparse_neighbor_fixture()
  graph <- fixture$graph
  full_diagonal <- mixed_diagonal <- graph
  Matrix::diag(full_diagonal) <- 1
  Matrix::diag(mixed_diagonal) <- rep(c(0, 2), nrow(graph) / 2L)
  edges <- summary(graph)
  explicit_zero <- Matrix::sparseMatrix(
    i = c(edges$i, 1L),
    j = c(edges$j, 8L),
    x = c(edges$x, 0),
    dims = dim(graph)
  )
  expect_true(any(explicit_zero@x == 0))

  duplicates <- Matrix::sparseMatrix(
    i = c(edges$i, edges$i, 1L, 1L),
    j = c(edges$j, edges$j, 8L, 8L),
    x = c(edges$x / 2, edges$x / 2, 2, -2),
    dims = dim(graph),
    repr = "T"
  )
  graphs <- list(
    column = graph,
    row = methods::as(graph, "RsparseMatrix"),
    triplet = methods::as(graph, "TsparseMatrix"),
    logical = graph != 0,
    pattern = methods::as(graph, "nsparseMatrix"),
    full_diagonal = full_diagonal,
    mixed_diagonal = mixed_diagonal,
    explicit_zero = explicit_zero,
    duplicates = duplicates
  )

  for (include_self in c(TRUE, FALSE)) {
    reference <- ltsa(
      fixture$X,
      nn_method = fixture$idx,
      include_self = include_self,
      output = "B"
    )
    for (candidate in graphs) {
      original <- candidate
      for (n_assembly_threads in c(1L, 2L)) {
        actual <- ltsa(
          fixture$X,
          nn_method = candidate,
          include_self = include_self,
          n_assembly_threads = n_assembly_threads,
          output = "B"
        )
        expect_sparse_equivalent(actual, reference, tolerance = 1e-11)
      }
      expect_identical(candidate, original)
    }
  }
})

test_that("symmetric sparse graphs include neighbors implied by symmetry", {
  fixture <- make_sparse_neighbor_fixture(c(-2L, -1L, 1L, 2L))
  graph <- Matrix::forceSymmetric(fixture$graph != 0)

  for (include_self in c(TRUE, FALSE)) {
    reference <- ltsa(
      fixture$X,
      nn_method = fixture$idx,
      include_self = include_self,
      output = "B"
    )
    actual <- ltsa(
      fixture$X,
      nn_method = graph,
      include_self = include_self,
      output = "B"
    )
    expect_sparse_equivalent(actual, reference, tolerance = 1e-11)
  }
})

test_that("sparse graphs infer neighborhood size and report precomputed input", {
  fixture <- make_sparse_neighbor_fixture()

  for (include_self in c(TRUE, FALSE)) {
    k <- 4L + as.integer(include_self)
    result <- ltsa(
      fixture$X,
      nn_method = fixture$graph,
      include_self = include_self,
      eig_method = "eig",
      output = "result",
      include_B = TRUE
    )
    expect_identical(result$assembly$n_neighbors, k)
    expect_identical(result$assembly$include_self, include_self)
    expect_identical(result$assembly$neighbor_source, "precomputed")
    expect_true(is.na(result$assembly$neighbor_elapsed))
    explicit <- ltsa(
      fixture$X,
      nn_method = fixture$graph,
      n_neighbors = k,
      include_self = include_self,
      output = "B"
    )
    expect_sparse_equivalent(explicit, result$B, tolerance = 0)
    expect_error(
      ltsa(
        fixture$X,
        nn_method = fixture$graph,
        n_neighbors = k + 1L,
        include_self = include_self,
        output = "B"
      ),
      "neighborhood size must match n_neighbors"
    )
  }
})

test_that("sparse graphs reject invalid dimensions, values, and neighborhood sizes", {
  fixture <- make_sparse_neighbor_fixture()
  X <- fixture$X
  graph <- fixture$graph

  expect_error(
    ltsa(X, nn_method = graph[, -1L], output = "B"),
    "must be square"
  )
  expect_error(
    ltsa(X[-1L, ], nn_method = graph, output = "B"),
    "one row per observation"
  )
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- graph
    bad[1L, 2L] <- value
    expect_error(
      ltsa(X, nn_method = bad, output = "B"),
      "must contain only finite values"
    )
  }

  # Preserve the total edge count so a premature reshape would silently reassign neighbors.
  ragged <- graph
  ragged[1L, 8L] <- 1
  ragged[2L, 3L] <- 0
  empty_row <- graph
  empty_row[6L, ] <- 0
  for (bad in list(ragged, empty_row)) {
    expect_error(
      ltsa(X, nn_method = bad, output = "B"),
      "same number of nonzero off-diagonal entries in every row"
    )
  }
  expect_error(
    ltsa(X, nn_method = graph * 0, output = "B"),
    "at least ndim \\+ 2"
  )
})

test_that("include_self affects whether sparse neighborhoods meet the minimum size", {
  fixture <- make_sparse_neighbor_fixture(seq_len(3L))

  expect_no_warning(ltsa(fixture$X, nn_method = fixture$graph, output = "B"))
  expect_error(
    ltsa(
      fixture$X,
      nn_method = fixture$graph,
      include_self = FALSE,
      output = "B"
    ),
    "at least ndim \\+ 2"
  )
})
