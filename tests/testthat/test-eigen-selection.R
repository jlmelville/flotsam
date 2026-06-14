selection_test_basis <- function() {
  list(
    u = c(1, 1, 1, 1) / 2,
    v1 = c(1, -1, 0, 0) / sqrt(2),
    v2 = c(0, 0, 1, -1) / sqrt(2),
    v3 = c(1, 1, -1, -1) / 2
  )
}

selection_test_matrix <- function(basis) {
  Q <- do.call(cbind, basis)
  Q %*% diag(c(0, 1, 2, 3)) %*% t(Q)
}

expect_selected_basis <- function(selected, expected) {
  expect_equal(abs(crossprod(selected, expected)), diag(ncol(expected)), tolerance = 1e-12)
}

test_that("embedding vector selection drops a returned trivial vector", {
  basis <- selection_test_basis()
  B <- selection_test_matrix(basis)
  candidates <- cbind(basis$v2, basis$u, basis$v3, basis$v1)

  selected <- flotsam:::select_ltsa_embedding_vectors(B, candidates, ndim = 2L)

  expect_selected_basis(selected, cbind(basis$v1, basis$v2))
})

test_that("embedding vector selection keeps smallest vectors when trivial vector is absent", {
  basis <- selection_test_basis()
  B <- selection_test_matrix(basis)
  candidates <- cbind(basis$v2, basis$v3, basis$v1)

  selected <- flotsam:::select_ltsa_embedding_vectors(B, candidates, ndim = 2L)

  expect_selected_basis(selected, cbind(basis$v1, basis$v2))
})
