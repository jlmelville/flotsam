test_that("sparse normalization implements Jacobi scaling", {
  # Protects the algebra behind the public normalize = TRUE path.
  L <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 3L, 2L, 3L),
    x = c(4, 2, 2, 9, 3, 3, 16),
    dims = c(3L, 3L),
    giveCsparse = TRUE
  )

  normalized <- flotsam:::ltsa_normalize_sparse_operator(L)
  diagonal <- diag(L)
  Dinvs <- sqrt(1 / diagonal)
  reference <- Matrix::Diagonal(x = Dinvs) %*%
    L %*%
    Matrix::Diagonal(x = Dinvs)

  expect_equal(as.matrix(normalized$Lsym), as.matrix(reference))
  expect_equal(normalized$Dinvs, Dinvs)
  expect_equal(normalized$nullvec, sqrt(diagonal))
})

test_that("normalized public results expose the generalized algebra", {
  result <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = TRUE,
    normalize = TRUE,
    eig_method = "eig",
    eig_k = 4L,
    output = "result",
    include_B = TRUE
  )

  details <- result$eigen$normalized_details
  expect_named(
    details,
    c(
      "mass",
      "mass_summary",
      "symmetric_embedding",
      "generalized_absolute_residuals",
      "generalized_residual_scale",
      "generalized_residuals",
      "weighted_orthogonality_error",
      "weighted_centering_error",
      "map_back_error",
      "reverse_occurrence",
      "component_embedding_overlap"
    )
  )
  expect_false(any(c("embedding", "mapped_embedding") %in% names(details)))

  B <- result$B
  V <- result$embedding
  U <- details$symmetric_embedding
  D <- details$mass
  values <- result$eigen$values
  generalized_residual <- as.matrix(B %*% V) -
    D * sweep(V, 2L, values, "*")
  absolute_residuals <- sqrt(colSums(generalized_residual^2))
  residual_scale <- sqrt(max(D)) * max(result$eigen$lambda_max, 1)

  expect_equal(D, diag(B), tolerance = 0)
  expect_equal(
    as.matrix(B %*% V),
    D * sweep(V, 2L, values, "*"),
    tolerance = 1e-10
  )
  expect_equal(crossprod(V, D * V), diag(ncol(V)), tolerance = 1e-10)
  expect_equal(as.numeric(crossprod(V, D)), rep(0, ncol(V)), tolerance = 1e-10)
  expect_equal(U, sqrt(D) * V, tolerance = 1e-12)
  expect_equal(details$generalized_absolute_residuals, absolute_residuals)
  expect_equal(details$generalized_residual_scale, residual_scale)
  expect_equal(
    details$generalized_residuals,
    absolute_residuals / residual_scale
  )
  expect_lt(details$weighted_orthogonality_error, 1e-10)
  expect_lt(details$weighted_centering_error, 1e-10)
  expect_lt(details$map_back_error, 1e-12)
  expect_equal(
    details$reverse_occurrence$correlation_with_mass,
    stats::cor(details$reverse_occurrence$counts, D)
  )

  summary <- details$mass_summary
  expect_named(
    summary,
    c(
      "quantiles",
      "max_to_min_ratio",
      "min_to_median_ratio",
      "log10_max_to_min_ratio",
      "log10_min_to_median_ratio",
      "index_limit",
      "smallest_mass_indices",
      "largest_mass_indices"
    )
  )
  expect_length(summary$quantiles, 5L)
  expect_equal(summary$max_to_min_ratio, max(D) / min(D))
  expect_equal(summary$min_to_median_ratio, min(D) / stats::median(D))
  expect_equal(
    summary$log10_max_to_min_ratio,
    log10(max(D)) - log10(min(D))
  )
  expect_equal(
    summary$log10_min_to_median_ratio,
    log10(min(D)) - log10(stats::median(D))
  )
  index_limit <- min(length(D), summary$index_limit)
  expect_identical(
    summary$smallest_mass_indices,
    head(order(D, seq_along(D)), index_limit)
  )
  expect_identical(
    summary$largest_mass_indices,
    head(order(-D, seq_along(D)), index_limit)
  )
})

test_that("normalized component overlap uses symmetric weighted contrasts", {
  set.seed(20260807)
  X <- matrix(rnorm(12L * 3L), nrow = 12L)
  # fmt: skip
  nn_idx <- matrix(
    c(
       1L,  2L,  3L,  4L,
       2L,  3L,  4L,  5L,
       3L,  4L,  5L,  6L,
       4L,  5L,  6L,  1L,
       5L,  6L,  1L,  2L,
       6L,  1L,  2L,  3L,
       7L,  8L,  9L, 10L,
       8L,  9L, 10L, 11L,
       9L, 10L, 11L, 12L,
      10L, 11L, 12L,  7L,
      11L, 12L,  7L,  8L,
      12L,  7L,  8L,  9L
    ),
    nrow = 12L,
    byrow = TRUE
  )

  result <- ltsa(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    normalize = TRUE,
    eig_method = "eig",
    output = "result"
  )

  details <- result$eigen$normalized_details
  overlap <- details$component_embedding_overlap
  membership <- result$assembly$component_membership
  D <- details$mass
  U <- details$symmetric_embedding
  indicator <- model.matrix(~ factor(membership) - 1L)
  component_mass <- colSums(D * indicator)
  weighted_indicators <- sqrt(D) *
    sweep(
      indicator,
      2L,
      sqrt(component_mass),
      "/"
    )
  global_direction <- sqrt(D / sum(D))
  contrasts <- weighted_indicators -
    global_direction %*% crossprod(global_direction, weighted_indicators)
  contrast_decomposition <- qr(contrasts)
  contrast_basis <- qr.Q(contrast_decomposition)[,
    seq_len(contrast_decomposition$rank),
    drop = FALSE
  ]
  embedding_decomposition <- qr(U)
  embedding_basis <- qr.Q(embedding_decomposition)[,
    seq_len(embedding_decomposition$rank),
    drop = FALSE
  ]
  cosines <- svd(crossprod(contrast_basis, embedding_basis), nu = 0, nv = 0)$d

  expect_identical(result$assembly$component_count, 2L)
  expect_false("component_embedding_overlap" %in% names(result$assembly))
  expect_equal(overlap$principal_angle_cosines, cosines, tolerance = 1e-10)
  expect_equal(overlap$projection_energy, sum(cosines^2), tolerance = 1e-10)
  expect_equal(details$reverse_occurrence$counts, rep.int(4L, 12L))
  expect_true(is.na(details$reverse_occurrence$correlation_with_mass))
})

test_that("reverse occurrence uses effective assembly members", {
  set.seed(20260808)
  X <- matrix(rnorm(8L * 3L), nrow = 8L)
  # fmt: skip
  nn_idx <- matrix(
    c(
      1L, 2L, 3L, 4L, 5L,
      2L, 1L, 3L, 4L, 5L,
      3L, 1L, 2L, 4L, 5L,
      4L, 1L, 2L, 3L, 5L,
      5L, 1L, 2L, 3L, 4L,
      6L, 1L, 2L, 7L, 8L,
      7L, 1L, 2L, 6L, 8L,
      8L, 1L, 2L, 6L, 7L
    ),
    nrow = 8L,
    byrow = TRUE
  )

  result <- ltsa(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = FALSE,
    normalize = TRUE,
    eig_method = "eig",
    output = "result"
  )

  details <- result$eigen$normalized_details
  expected_counts <- tabulate(nn_idx[, -1L], nbins = nrow(X))
  expect_identical(details$reverse_occurrence$counts, expected_counts)
  expect_identical(sum(details$reverse_occurrence$counts), 8L * 4L)
  expect_equal(
    details$reverse_occurrence$correlation_with_mass,
    stats::cor(expected_counts, details$mass)
  )
})

test_that("output B remains the unnormalized alignment operator", {
  args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = TRUE,
    output = "B"
  )

  unnormalized <- do.call(ltsa, c(args, list(normalize = FALSE)))
  normalized_request <- do.call(ltsa, c(args, list(normalize = TRUE)))

  expect_equal(
    as.matrix(normalized_request),
    as.matrix(unnormalized),
    tolerance = 0
  )
})
