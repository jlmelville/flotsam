test_that("invalid methods are rejected early", {
  expect_error(
    ltsa(iris[1:10, ], nn_method = "bad"),
    "nn_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = "bad"),
    "eig_method"
  )
})

test_that("eigen is accepted as an eig alias", {
  iris10_eig <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "eig"
  )
  iris10_eigen <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "eigen"
  )
  expect_equal(iris10_eigen, iris10_eig)
})

test_that("input data must be numeric and finite", {
  expect_error(
    ltsa(data.frame(group = letters[1:4])),
    "at least one numeric column"
  )
  expect_error(
    ltsa(matrix(letters[1:4], ncol = 2)),
    "numeric values"
  )
  expect_error(
    ltsa(matrix(c(1, 2, NA, 4), ncol = 2)),
    "finite"
  )
})

test_that("dimension and neighborhood arguments are validated", {
  expect_error(
    ltsa(iris[1:10, ], ndim = 0),
    "ndim"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 10),
    "ndim must be less than the number of observations"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 2, ndim = 2),
    "greater than ndim"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 11, include_self = TRUE),
    "too large"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 10, include_self = FALSE),
    "too large"
  )
  expect_error(
    ltsa(iris[1:10, ], n_threads = -1),
    "n_threads"
  )
})

test_that("logical arguments must be scalar TRUE or FALSE", {
  expect_error(
    ltsa(iris[1:10, ], include_self = NA),
    "include_self"
  )
  expect_error(
    ltsa(iris[1:10, ], normalize = c(TRUE, FALSE)),
    "normalize"
  )
  expect_error(
    ltsa(iris[1:10, ], ret_B = 1),
    "ret_B"
  )
  expect_error(
    ltsa(iris[1:10, ], verbose = NA),
    "verbose"
  )
})

test_that("eigenanalysis errors are not swallowed", {
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = "svdr",
      not_an_argument = TRUE
    ),
    "unused argument"
  )
})
