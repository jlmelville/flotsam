test_that("explicit RSpectra path unfolds a swiss roll", {
  set.seed(42)
  n <- 500L
  phi <- stats::runif(n, min = 1.5 * pi, max = 4.5 * pi)
  swiss_roll <- data.frame(
    x = phi * cos(phi),
    y = phi * sin(phi),
    z = stats::runif(n, max = 10)
  )

  swiss_embedding <- NULL
  expect_warning(
    swiss_embedding <- ltsa(
      swiss_roll,
      nn_method = "exact",
      eig_method = "rspectra",
      n_threads = 0
    ),
    "LTSA eigenanalysis status is warning:.*Weak Ritz boundary gap"
  )

  truth <- scale(cbind(phi = phi, z = swiss_roll$z))
  agreement <- stats::cancor(scale(swiss_embedding), truth)$cor

  expect_gt(min(agreement), 0.95)
})
