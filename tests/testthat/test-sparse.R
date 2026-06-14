test_that("find_in_spsq maps neighborhood pairs to sparse slots", {
  m <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = 1,
    dims = c(2, 2)
  )

  expect_equal(flotsam:::find_in_spsq(m, c(1L, 2L)), 1:4)
})

test_that("find_in_spsq rejects missing sparse entries", {
  m <- Matrix::sparseMatrix(
    i = c(1, 2),
    j = c(1, 2),
    x = 1,
    dims = c(2, 2)
  )

  expect_error(
    flotsam:::find_in_spsq(m, c(1L, 2L)),
    "does not contain requested neighborhood pair"
  )
})

test_that("find_in_spsq rejects out-of-range indices", {
  m <- Matrix::sparseMatrix(
    i = c(1, 2),
    j = c(1, 2),
    x = 1,
    dims = c(2, 2)
  )

  expect_error(
    flotsam:::find_in_spsq(m, c(1L, 3L)),
    "outside the sparse matrix dimensions"
  )
})
