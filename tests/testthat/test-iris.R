test_that("iris", {
  irisnn50_ltsa <-
    ltsa(
      iris,
      nn_method = "exact",
      n_neighbors = 50,
      include_self = FALSE,
      eig_method = "eig"
    )
  expect_equal(abs(irisnn50_ltsa), abs(irisnn50_ltsa_expected), tolerance = 1e-2)
})
