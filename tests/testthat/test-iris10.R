expect_iris10_backends_match_dense <- function(normalize, include_self) {
  common_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = include_self,
    normalize = normalize,
    eig_k = 4L,
    output = "result"
  )
  dense_reference <- do.call(
    ltsa,
    c(common_args, list(eig_method = "eig"))
  )$embedding

  for (method in c("rspectra", "irlba", "svdr")) {
    set.seed(7)
    result <- do.call(
      ltsa,
      c(
        common_args,
        list(eig_method = method)
      )
    )

    expect_identical(result$eigen$backend$name, method)
    expect_same_subspace(
      result$embedding,
      dense_reference,
      tolerance = 1e-5
    )
  }
}

test_that("iris10 iterative backends match the dense reference subspace", {
  expect_iris10_backends_match_dense(
    normalize = FALSE,
    include_self = FALSE
  )
})

test_that("normalized iris10 backends match the dense reference subspace", {
  expect_iris10_backends_match_dense(
    normalize = TRUE,
    include_self = TRUE
  )
})

expect_public_result <- function(
  result,
  method,
  normalized = FALSE,
  n = 10L,
  ndim = 2L,
  eig_k = 4L,
  row_names = as.character(seq_len(n))
) {
  expect_named(
    result,
    c("embedding", "eigen", "assembly"),
    ignore.order = TRUE
  )
  expect_true(is.matrix(result$embedding))
  expect_equal(dim(result$embedding), c(n, ndim))
  coordinate_dimnames <- list(
    row_names,
    paste0("LTSA", seq_len(ndim))
  )
  expect_identical(dimnames(result$embedding), coordinate_dimnames)
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
  expect_named(result$eigen, eigen_fields, ignore.order = TRUE)
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
    "near_zero_block_truncated",
    ignore.order = TRUE
  )
  expect_type(result$eigen$diagnostics$near_zero_block_truncated, "logical")
  expect_false("B" %in% names(result))
  assembly_fields <- c(
    "n_neighbors",
    "include_self",
    "neighbor_source",
    "neighbor_elapsed",
    "rank_deficient_count",
    "min_local_rank",
    "assembly_route",
    "component_count",
    "component_sizes",
    "component_membership",
    "requested_assembly_threads",
    "effective_assembly_threads"
  )
  expect_named(result$assembly, assembly_fields, ignore.order = TRUE)
  expect_identical(result$assembly$n_neighbors, 8L)
  expect_identical(result$assembly$neighbor_source, "exact")
  expect_type(result$assembly$neighbor_elapsed, "double")
  expect_length(result$assembly$neighbor_elapsed, 1L)
  expect_true(is.finite(result$assembly$neighbor_elapsed))
  expect_gte(result$assembly$neighbor_elapsed, 0)
}

test_that("default public return remains an embedding matrix", {
  iris10_embedding <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE
  )

  expect_true(is.matrix(iris10_embedding))
  expect_equal(dim(iris10_embedding), c(10L, 2L))
  expect_identical(
    dimnames(iris10_embedding),
    list(as.character(seq_len(10L)), c("LTSA1", "LTSA2"))
  )
})

test_that("public coordinate matrices preserve row names and name dimensions", {
  named_iris <- as.matrix(iris[1:10, 1:4])
  rownames(named_iris) <- paste0("sample-", seq_len(nrow(named_iris)))
  expected_dimnames <- list(
    rownames(named_iris),
    c("LTSA1", "LTSA2")
  )

  embedding <- ltsa(
    named_iris,
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_method = "eig"
  )
  expect_identical(dimnames(embedding), expected_dimnames)

  normalized_result <- ltsa(
    named_iris,
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = TRUE,
    normalize = TRUE,
    eig_method = "eig",
    eig_k = 4L,
    output = "result"
  )
  expect_identical(dimnames(normalized_result$embedding), expected_dimnames)

  unnamed_iris <- named_iris
  rownames(unnamed_iris) <- NULL
  unnamed_embedding <- ltsa(
    unnamed_iris,
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_method = "eig"
  )
  expect_identical(
    dimnames(unnamed_embedding),
    list(NULL, c("LTSA1", "LTSA2"))
  )

  B <- ltsa(
    named_iris,
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = TRUE,
    output = "B"
  )
  expect_identical(dimnames(B), list(NULL, NULL))
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
  methods <- c("auto", "rspectra", "irlba", "svdr", "eig", "eigen")

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

    expect_public_result(
      result,
      method = if (identical(method, "eigen")) "eig" else method
    )
  }
})

test_that("normalized iterative detailed results have consistent shape", {
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
      output = "result"
    )

    expect_public_result(
      result,
      method = method,
      normalized = TRUE
    )
  }
})
