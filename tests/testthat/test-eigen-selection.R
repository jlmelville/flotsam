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
  expect_equal(
    abs(crossprod(selected, expected)),
    diag(ncol(expected)),
    tolerance = 1e-12
  )
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

test_that("LTSA symmetrization returns a general sparse solver matrix", {
  B <- Matrix::sparseMatrix(
    i = c(1L, 2L, 2L, 3L),
    j = c(1L, 1L, 2L, 3L),
    x = c(2, 4, 6, 8),
    dims = c(3L, 3L)
  )

  B_sym <- flotsam:::symmetrize_ltsa_matrix(B)

  expect_s4_class(B_sym, "dgCMatrix")
  expect_false(methods::is(B_sym, "dsCMatrix"))
  expect_true(Matrix::isSymmetric(B_sym))
  expect_equal(as.matrix(B_sym), 0.5 * (as.matrix(B) + t(as.matrix(B))))
})

test_that("dense LTSA eigensolver fallback reports small residuals", {
  B <- Matrix::Diagonal(x = seq(0, 11))

  res <- flotsam:::rs_eig(B, k = 6L)

  expect_identical(res$backend, "dense_eigen")
  expect_identical(res$eig_k, 6L)
  expect_identical(res$nconv, 6L)
  expect_equal(res$values, seq(0, 5), tolerance = 1e-12)
  expect_lt(max(res$scaled_residuals), 1e-12)
})

test_that("RSpectra path uses shifted largest-algebraic solve with residual metadata", {
  B <- Matrix::Diagonal(x = seq(0, 29))

  res <- flotsam:::rs_eig(
    B,
    k = 6L,
    dense_n = 0L,
    tol = 1e-10,
    maxitr = 5000L
  )

  expect_identical(res$backend, "rspectra")
  expect_identical(res$solve_which, "LA")
  expect_identical(res$eig_k, 6L)
  expect_gte(res$nconv, 6L)
  expect_equal(res$values, seq(0, 5), tolerance = 1e-8)
  expect_lt(max(res$scaled_residuals), 1e-8)
  expect_s4_class(res$matrix, "dgCMatrix")
})

test_that("RSpectra partial convergence is a hard LTSA error", {
  set.seed(1)
  A <- crossprod(matrix(stats::rnorm(80L * 80L), nrow = 80L))
  B <- Matrix::Matrix(A, sparse = TRUE)
  lambda_max <- max(eigen(A, symmetric = TRUE, only.values = TRUE)$values)

  expect_error(
    suppressWarnings(
      flotsam:::rs_eig(
        B,
        k = 20L,
        lambda_max = lambda_max,
        dense_n = 0L,
        ncv = 21L,
        maxitr = 1L,
        tol = 1e-16
      )
    ),
    "RSpectra failed to converge enough LTSA candidate vectors: 0 / 20"
  )
})
