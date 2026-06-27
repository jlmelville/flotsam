test_that("invalid methods are rejected early", {
  expect_error(
    ltsa(iris[1:10, ], nn_method = "bad"),
    "nn_method"
  )
  expect_error(
    ltsa(iris[1:10, ], nn_method = c("exact", "nnd")),
    "nn_method"
  )
  expect_error(
    ltsa(iris[1:10, ], nn_method = NA_character_),
    "nn_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = "bad"),
    "eig_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = c("eig", "eigen")),
    "eig_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = NA_character_),
    "eig_method"
  )
  expect_error(
    ltsa(iris[1:10, ], output = "bad"),
    "should be one of"
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
    ltsa(seq_len(4)),
    "matrix or data frame"
  )
  expect_error(
    ltsa(data.frame(group = letters[1:4])),
    "at least one numeric column"
  )
  mixed_df <- flotsam:::x2m(data.frame(x = 1:3, group = letters[1:3]))
  expect_equal(unname(mixed_df), matrix(as.double(1:3), ncol = 1))
  expect_identical(colnames(mixed_df), "x")
  expect_error(
    ltsa(matrix(letters[1:4], ncol = 2)),
    "numeric values"
  )
  expect_error(
    ltsa(matrix(c(1, 2, NA, 4), ncol = 2)),
    "finite"
  )
  expect_error(
    ltsa(matrix(numeric(), nrow = 2, ncol = 0)),
    "at least one column"
  )
  expect_error(
    ltsa(matrix(1, nrow = 1, ncol = 2)),
    "at least two observations"
  )
})

test_that("dimension and neighborhood arguments are validated", {
  expect_error(
    ltsa(iris[1:10, ], ndim = 0),
    "ndim"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 1.5),
    "ndim"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 10),
    "ndim must be less than the number of observations"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 2, eig_k = 2L),
    "eig_k"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 2, eig_k = 10L),
    "eig_k"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 2.5),
    "n_neighbors"
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
  expect_error(
    ltsa(iris[1:10, ], n_threads = 1.5),
    "n_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = 0),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = -1),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = 1.5),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = NA_real_),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = c(1, 2)),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], copy_max_mib = -1),
    "copy_max_mib"
  )
  expect_error(
    ltsa(iris[1:10, ], copy_max_mib = NA_real_),
    "copy_max_mib"
  )
  expect_error(
    ltsa(iris[1:10, ], copy_max_mib = Inf),
    "copy_max_mib"
  )
  expect_error(
    ltsa(iris[1:10, ], copy_max_mib = c(128, 256)),
    "copy_max_mib"
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
    ltsa(iris[1:10, ], include_B = 1),
    "include_B"
  )
  expect_error(
    ltsa(iris[1:10, ], verbose = NA),
    "verbose"
  )
})

test_that("eigen diagnostic tolerances must be finite positive numbers", {
  expect_error(
    ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8L, resid_tol = 0),
    "resid_tol must be a finite positive number"
  )
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8L,
      resid_tol = c(1e-5, 1e-4)
    ),
    "resid_tol must be a finite positive number"
  )
  expect_error(
    ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8L, gap_tol = Inf),
    "gap_tol must be a finite positive number"
  )
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8L,
      gap_tol = NA_real_
    ),
    "gap_tol must be a finite positive number"
  )
})

test_that("normalized LTSA rejects matrices with non-positive diagonal", {
  expect_error(
    ltsa(
      iris[1:4, ],
      nn_method = "exact",
      n_neighbors = 3L,
      include_self = TRUE,
      ndim = 2L,
      normalize = TRUE,
      eig_method = "eig"
    ),
    "Cannot normalize the LTSA matrix"
  )
})

test_that("eigenanalysis errors include ltsa context and solver details", {
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = "svdr",
      not_an_argument = TRUE
    ),
    "Eigenanalysis failed:.*unused argument"
  )
})
