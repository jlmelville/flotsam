test_that("iris10", {
  iris10_ltsa_expected <- matrix(
    byrow = TRUE, ncol = 2,
    c(
      0.23438576, 0.09178531,
      -0.14825134, 0.53317022,
      -0.17014801, -0.13276554,
      -0.2669889, -0.06793736,
      0.22719947, -0.2243461,
      0.69774771, -0.05138839,
      -0.11459417, -0.6165803,
      0.12743673, 0.11721779,
      -0.49692296, -0.12175932,
      -0.08986429, 0.47260369
    )
  )

  iris10_ltsa <- ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8, include_self = FALSE, eig_method = "rspectra")
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <- ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8, include_self = FALSE, eig_method = "irlba")
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <- ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8, include_self = FALSE, eig_method = "svdr")
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <- ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8, include_self = FALSE, eig_method = "eig")
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <- ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8, include_self = FALSE, eig_method = "fullsvd")
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)
})

test_that("normalized", {
  iris10_ltsa_expected <- matrix(
    c(
      -0.06988, -0.1593,
      0.2234, -0.04539,
      -0.01874, 0.1182,
      0.0399, 0.1608,
      -0.185, -0.08534,
      -0.2591, -0.4034,
      -0.2188, 0.1961,
      -0.02322, -0.1096,
      0.09676, 0.3019,
      0.1833, -0.06714
    ),
    byrow = TRUE, ncol = 2
  )
  iris10_ltsa <-
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = TRUE,
      normalize = TRUE,
      eig_method = "rspectra"
    )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <-
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = TRUE,
      normalize = TRUE,
      eig_method = "irlba"
    )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <-
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = TRUE,
      normalize = TRUE,
      eig_method = "svdr"
    )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <-
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = TRUE,
      normalize = TRUE,
      eig_method = "eig"
    )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <-
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = TRUE,
      normalize = TRUE,
      eig_method = "fullsvd"
    )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)
})
