expect_psd_sparse <- function(B, tolerance = 1e-8) {
  expect_true(Matrix::isSymmetric(B))
  expect_gte(min(diag(B)), -tolerance)
  vals <- eigen(as.matrix(B), symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(vals), -tolerance)
}

test_that("constant neighborhoods warn and use valid local projections", {
  expect_warning(
    B <- ltsa(
      matrix(1, nrow = 8, ncol = 3),
      nn_method = "exact",
      n_neighbors = 4,
      eig_method = "eig",
      ret_B = TRUE
    ),
    "local neighborhoods had numerical rank below ndim"
  )

  expect_psd_sparse(B)
})

test_that("collinear neighborhoods warn and use valid local projections", {
  x <- seq_len(8)
  expect_warning(
    B <- ltsa(
      cbind(x, 2 * x, 3 * x),
      nn_method = "exact",
      n_neighbors = 4,
      eig_method = "eig",
      ret_B = TRUE
    ),
    "local neighborhoods had numerical rank below ndim"
  )

  expect_psd_sparse(B)
})
