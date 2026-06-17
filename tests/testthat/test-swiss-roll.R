test_that("default RSpectra path unfolds a swiss roll", {
  set.seed(42)
  n <- 500L
  phi <- stats::runif(n, min = 1.5 * pi, max = 4.5 * pi)
  swiss_roll <- data.frame(
    x = phi * cos(phi),
    y = phi * sin(phi),
    z = stats::runif(n, max = 10)
  )

  swiss_ltsa <- NULL
  expect_warning(
    swiss_ltsa <- ltsa(
      swiss_roll,
      nn_method = "exact",
      eig_method = "rspectra",
      n_threads = 0
    ),
    NA
  )

  truth <- scale(cbind(phi = phi, z = swiss_roll$z))
  agreement <- stats::cancor(scale(swiss_ltsa), truth)$cor

  expect_gt(min(agreement), 0.95)
})
