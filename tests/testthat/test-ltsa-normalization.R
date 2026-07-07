test_that("sparse normalization implements Jacobi scaling", {
  # Protects the algebra behind the public normalize = TRUE path.
  L <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 3L, 2L, 3L),
    x = c(4, 2, 2, 9, 3, 3, 16),
    dims = c(3L, 3L),
    giveCsparse = TRUE
  )

  normalized <- flotsam:::ltsa_normalize_sparse_operator(L)
  diagonal <- diag(L)
  Dinvs <- sqrt(1 / diagonal)
  reference <- Matrix::Diagonal(x = Dinvs) %*%
    L %*%
    Matrix::Diagonal(x = Dinvs)

  expect_equal(as.matrix(normalized$Lsym), as.matrix(reference))
  expect_equal(normalized$Dinvs, Dinvs)
  expect_equal(normalized$nullvec, sqrt(diagonal))
})
