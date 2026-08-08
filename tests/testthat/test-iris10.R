test_that("iris10", {
  # fmt: skip
  iris10_ltsa_expected <- matrix(
    c(
       0.23438576,  0.09178531,
      -0.14825134,  0.53317022,
      -0.17014801, -0.13276554,
      -0.26698890, -0.06793736,
       0.22719947, -0.22434610,
       0.69774771, -0.05138839,
      -0.11459417, -0.61658030,
       0.12743673,  0.11721779,
      -0.49692296, -0.12175932,
      -0.08986429,  0.47260369
    ),
    ncol = 2,
    byrow = TRUE
  )

  iris10_ltsa <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "rspectra"
  )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "irlba"
  )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "svdr"
  )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)

  iris10_ltsa <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "eig"
  )
  expect_equal(abs(iris10_ltsa), abs(iris10_ltsa_expected), tolerance = 1e-2)
})

test_that("normalized", {
  # A spare exact-search candidate changes include-self boundary ties before
  # deterministic self canonicalization.
  # fmt: skip
  iris10_ltsa_expected <- matrix(
    c(
      -0.07766,  0.15450,
       0.22300,  0.04879,
      -0.01936, -0.11990,
       0.04053, -0.16160,
      -0.19330,  0.07628,
      -0.27560,  0.39040,
      -0.22190, -0.20310,
      -0.02942,  0.10530,
       0.10190, -0.29940,
       0.18200,  0.06882
    ),
    ncol = 2,
    byrow = TRUE
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
})

expect_ltsa_public_result <- function(
  result,
  method,
  normalized = FALSE,
  n = 10L,
  ndim = 2L,
  eig_k = 4L
) {
  expect_named(result, c("embedding", "eigen", "assembly"))
  expect_true(is.matrix(result$embedding))
  expect_equal(dim(result$embedding), c(n, ndim))
  eigen_fields <- c(
    "method",
    "normalized",
    "eig_k",
    "values",
    "ritz_values",
    "residuals",
    "rank",
    "lambda_max",
    "status",
    "messages",
    "backend",
    "diagnostics"
  )
  if (normalized) {
    eigen_fields <- c(eigen_fields, "normalized_details")
  }
  expect_named(result$eigen, eigen_fields)
  expect_identical(result$eigen$method, method)
  expect_identical(result$eigen$normalized, normalized)
  expect_identical(result$eigen$eig_k, eig_k)
  expect_length(result$eigen$values, ndim)
  expect_length(result$eigen$residuals, ndim)
  expect_gte(length(result$eigen$ritz_values), ndim)
  expect_gte(result$eigen$rank, ndim)
  expect_true(result$eigen$status %in% c("ok", "warning", "invalid"))
  expect_true(is.list(result$eigen$backend))
  expect_named(
    result$eigen$diagnostics,
    c(
      "selected_boundary_gap",
      "global_gap",
      "local_gap",
      "candidate_span_size",
      "near_zero_nonconstant_count",
      "near_zero_nonconstant_counts",
      "near_zero_threshold",
      "near_zero_thresholds",
      "near_zero_boundary_gap",
      "near_zero_boundary_gaps",
      "near_zero_boundary_observed",
      "near_zero_boundaries_observed",
      "near_zero_block_truncated"
    )
  )
  expect_type(
    result$eigen$diagnostics$near_zero_nonconstant_count,
    "integer"
  )
  expect_type(result$eigen$diagnostics$near_zero_threshold, "double")
  expect_type(result$eigen$diagnostics$near_zero_block_truncated, "logical")
  expect_false("B" %in% names(result))
  expect_identical(result$assembly$n_neighbors, 8L)
  expect_identical(result$assembly$neighbor_source, "exact")
  expect_type(result$assembly$neighbor_elapsed, "double")
  expect_length(result$assembly$neighbor_elapsed, 1L)
  expect_true(is.finite(result$assembly$neighbor_elapsed))
  expect_gte(result$assembly$neighbor_elapsed, 0)
}

test_that("default public return remains an embedding matrix", {
  iris10_ltsa <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "rspectra"
  )

  expect_true(is.matrix(iris10_ltsa))
  expect_equal(dim(iris10_ltsa), c(10L, 2L))
})

test_that("public embedding status is actionable and result status is inspectable", {
  warning_embedding <- NULL
  expect_warning(
    warning_embedding <- ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = "eig",
      eig_k = 3L
    ),
    "LTSA eigenanalysis status is warning:.*no spare boundary"
  )
  expect_true(is.matrix(warning_embedding))

  warning_result <- expect_silent(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = "eig",
      eig_k = 3L,
      output = "result"
    )
  )
  expect_identical(warning_result$eigen$status, "warning")

  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = "eig",
      resid_tol = 1e-20
    ),
    "LTSA eigenanalysis status is invalid"
  )

  invalid_result <- expect_silent(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = "eig",
      resid_tol = 1e-20,
      output = "result"
    )
  )
  expect_identical(invalid_result$eigen$status, "invalid")
})

test_that("all solver methods support detailed public results", {
  methods <- c("rspectra", "irlba", "svdr", "eig", "eigen")

  for (method in methods) {
    set.seed(7)
    result <- ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = method,
      eig_k = 4L,
      output = "result"
    )

    expect_ltsa_public_result(
      result,
      method = if (identical(method, "eigen")) "eig" else method
    )
  }
})

test_that("normalized iterative detailed results have consistent diagnostics", {
  methods <- c("rspectra", "irlba", "svdr")

  for (method in methods) {
    set.seed(8)
    result <- ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = TRUE,
      normalize = TRUE,
      eig_method = method,
      eig_k = 4L,
      output = "result",
      dense_n = 0L
    )

    expect_ltsa_public_result(
      result,
      method = method,
      normalized = TRUE
    )
    expect_length(result$eigen$normalized_details$generalized_residuals, 2L)
    expect_false("component_embedding_overlap" %in% names(result$assembly))
  }
})
