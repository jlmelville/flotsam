selection_test_basis <- function() {
  list(
    u = c(1, 1, 1, 1) / 2,
    v1 = c(1, -1, 0, 0) / sqrt(2),
    v2 = c(0, 0, 1, -1) / sqrt(2),
    v3 = c(1, 1, -1, -1) / 2
  )
}

selection_test_matrix <- function(basis) {
  Q <- do.call(cbind, basis)
  Q %*% diag(c(0, 1, 2, 3)) %*% t(Q)
}

expect_selected_basis <- function(selected, expected) {
  expect_equal(
    abs(crossprod(selected, expected)),
    diag(ncol(expected)),
    tolerance = 1e-12
  )
}

orthonormalize_test_basis <- function(V) {
  qrV <- qr(as.matrix(V))
  qr.Q(qrV)[, seq_len(qrV$rank), drop = FALSE]
}

expect_same_subspace <- function(actual, expected, tolerance = 1e-10) {
  q_actual <- orthonormalize_test_basis(actual)
  q_expected <- orthonormalize_test_basis(expected)
  sv <- svd(crossprod(q_actual, q_expected), nu = 0, nv = 0)$d
  sv <- pmax(0, pmin(1, sv))
  projection_distance <- sqrt(max(
    0,
    ncol(q_actual) + ncol(q_expected) - 2 * sum(sv^2)
  ))

  expect_equal(ncol(q_actual), ncol(q_expected))
  expect_lt(projection_distance, tolerance)
}

centered_test_basis <- function(n, rank) {
  set.seed(123)
  nullvec <- rep(1, n) / sqrt(n)
  Z <- matrix(stats::rnorm(n * rank * 2L), nrow = n)
  Z <- Z - nullvec %*% crossprod(nullvec, Z)
  Q <- orthonormalize_test_basis(Z)
  Q[, seq_len(rank), drop = FALSE]
}

synthetic_ltsa_problem <- function(values) {
  n <- length(values)
  basis <- cbind(rep(1, n) / sqrt(n), centered_test_basis(n, n - 1L))
  list(
    matrix = basis %*% diag(values) %*% t(basis),
    basis = basis
  )
}

synthetic_ltsa_matrix <- function(values) {
  synthetic_ltsa_problem(values)$matrix
}

fixed_width_provider_factory <- function(
  problem,
  backend = "synthetic",
  lambda_max = NULL,
  convergence_known = TRUE,
  nconv = function(eig_k) NA_integer_,
  converged_columns = function(eig_k) {
    if (isTRUE(convergence_known)) eig_k else NA_integer_
  },
  cols = function(eig_k) seq_len(eig_k),
  niter = NA_integer_,
  nops = NA_integer_,
  mprod = NA_integer_
) {
  if (is.null(lambda_max)) {
    lambda_max <- max(eigen(
      problem$matrix,
      symmetric = TRUE,
      only.values = TRUE
    )$values)
  }
  lambda_max_value <- lambda_max
  calls <- data.frame(eig_k = integer())
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    calls <<- rbind(calls, data.frame(eig_k = eig_k))
    cols_i <- cols(eig_k)
    nconv_i <- if (is.function(nconv)) nconv(eig_k) else nconv
    converged_i <- if (is.function(converged_columns)) {
      converged_columns(eig_k)
    } else {
      converged_columns
    }

    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, cols_i, drop = FALSE],
      backend = backend,
      eig_k = eig_k,
      matrix = B,
      lambda_max = lambda_max_value,
      nconv = nconv_i,
      niter = niter,
      nops = nops,
      mprod = mprod,
      convergence_known = convergence_known,
      returned_columns = length(cols_i),
      converged_columns = converged_i
    )
  }

  list(
    provider = provider,
    calls = function() calls
  )
}

calibrated_synthetic_ritz <- function(values, ndim, lambda_max = max(values)) {
  problem <- synthetic_ltsa_problem(values)
  rr <- flotsam:::ltsa_ritz_select(
    problem$matrix,
    problem$basis,
    ndim = ndim,
    lambda_max = lambda_max
  )
  acceptance <- flotsam:::accept_ltsa_ritz(
    rr,
    ndim = ndim,
    lambda_max = lambda_max
  )
  flotsam:::ltsa_mark_partial_near_zero_block(
    list(
      values = rr$values,
      near_zero_nonconstant_count = rr$near_zero_nonconstant_count,
      near_zero_nonconstant_counts = rr$near_zero_nonconstant_counts,
      near_zero_tol = rr$near_zero_tol,
      acceptance = acceptance
    ),
    ndim = ndim
  )
}

synthetic_energy_rescue_result <- function(
  values,
  scaled_residuals,
  ndim,
  lambda_max = 1,
  rank_ok = TRUE,
  resid_ok = TRUE
) {
  residual_scale <- flotsam:::ltsa_residual_scale(lambda_max)
  values <- as.numeric(values)
  scaled_residuals <- as.numeric(scaled_residuals)
  if (length(scaled_residuals) == 1L) {
    scaled_residuals <- rep(scaled_residuals, length(values))
  }

  list(
    values = values,
    scaled_residuals = scaled_residuals,
    absolute_residuals = scaled_residuals * residual_scale,
    lambda_max = lambda_max,
    zero_tol = flotsam:::ltsa_gap_zero_tol(lambda_max),
    near_zero_tol = flotsam:::ltsa_near_zero_tol(lambda_max),
    near_zero_nonconstant_count = as.integer(sum(
      abs(values[seq_len(ndim)]) <= flotsam:::ltsa_near_zero_tol(lambda_max)
    )),
    acceptance = list(
      rank_ok = isTRUE(rank_ok),
      resid_ok = isTRUE(resid_ok),
      lambda_max = lambda_max
    )
  )
}

width_rescue_problem <- function(ndim, n = 96L) {
  low <- seq_len(ndim) * 1e-10
  high_count <- n - ndim - 1L
  high <- 5e-4 + seq_len(high_count) * 1e-8
  synthetic_ltsa_problem(c(0, low, high))
}

width_rescue_cols <- function(eig_k, ndim, n, mode) {
  low <- seq.int(2L, ndim + 1L)
  high <- seq.int(ndim + 2L, n)
  cols <- switch(
    mode,
    partial = c(1L, low[seq_len(max(1L, ndim - 1L))], high),
    complete = c(1L, low, high),
    high = c(1L, high),
    stop("unknown width rescue mode", call. = FALSE)
  )
  cols[seq_len(eig_k)]
}

width_rescue_provider_factory <- function(
  problem,
  ndim,
  ordinary_mode,
  strict_mode = function(eig_k) "complete",
  lambda_max = 1
) {
  calls <- data.frame(eig_k = integer(), strict = logical())
  provider <- function(
    B,
    eig_k,
    lambda_max = NULL,
    verbose = FALSE,
    strict = FALSE
  ) {
    calls <<- rbind(calls, data.frame(eig_k = eig_k, strict = strict))
    mode <- if (strict) strict_mode(eig_k) else ordinary_mode(eig_k)
    cols <- width_rescue_cols(
      eig_k = eig_k,
      ndim = ndim,
      n = nrow(problem$basis),
      mode = mode
    )
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, cols, drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = lambda_max,
      convergence_known = TRUE,
      returned_columns = length(cols),
      converged_columns = length(cols),
      nconv = length(cols)
    )
  }

  list(
    provider = provider,
    calls = function() calls
  )
}

run_width_rescue_case <- function(
  ndim = 2L,
  n = 96L,
  ordinary_mode,
  strict_mode = function(eig_k) "complete",
  ...
) {
  problem <- width_rescue_problem(ndim = ndim, n = n)
  fixture <- width_rescue_provider_factory(
    problem = problem,
    ndim = ndim,
    ordinary_mode = ordinary_mode,
    strict_mode = strict_mode
  )
  res <- flotsam:::ltsa_adaptive_ritz_eig(
    problem$matrix,
    ndim = ndim,
    provider = fixture$provider,
    gap_expansion_steps = 0L,
    width_first_rescue = TRUE,
    strict_rescue_arg_mapper = function(provider_args, ...) {
      c(provider_args, list(strict = TRUE))
    },
    ...
  )

  list(result = res, calls = fixture$calls(), problem = problem)
}

test_that("embedding vector selection drops a returned trivial vector", {
  basis <- selection_test_basis()
  B <- selection_test_matrix(basis)
  candidates <- cbind(basis$v2, basis$u, basis$v3, basis$v1)

  selected <- flotsam:::select_ltsa_embedding_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 3
  )

  expect_selected_basis(selected, cbind(basis$v1, basis$v2))
})

test_that("embedding vector selection keeps smallest vectors when trivial vector is absent", {
  basis <- selection_test_basis()
  B <- selection_test_matrix(basis)
  candidates <- cbind(basis$v2, basis$v3, basis$v1)

  selected <- flotsam:::select_ltsa_embedding_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 3
  )

  expect_selected_basis(selected, cbind(basis$v1, basis$v2))
})

test_that("Ritz selection is invariant to rotations in the candidate subspace", {
  basis <- selection_test_basis()
  B <- selection_test_matrix(basis)
  target <- cbind(basis$v1, basis$v2)
  candidate_space <- cbind(basis$u, target, basis$v3)
  # fmt: skip
  rotation <- qr.Q(qr(matrix(c(
    1, 2, 0, 1,
    -1, 1, 2, 0,
    0, 1, 1, -2,
    2, 0, -1, 1
  ), nrow = 4L)))
  candidates <- candidate_space %*% rotation

  rr <- flotsam:::ltsa_ritz_select(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 3
  )

  expect_same_subspace(rr$vectors, target)
  expect_equal(rr$values, c(1, 2), tolerance = 1e-12)
  expect_lt(max(rr$scaled_residuals), 1e-12)
  expect_identical(rr$rank_after_null, 3L)
  expect_equal(rr$boundary_gap, 1, tolerance = 1e-12)
})

test_that("clustered low-eigenvalue subspaces are compared by projector", {
  basis <- selection_test_basis()
  Q <- do.call(cbind, basis)
  B <- Q %*% diag(c(0, 1, 1, 4)) %*% t(Q)
  target <- cbind(basis$v1, basis$v2)
  candidates <- cbind(
    basis$u,
    (basis$v1 + basis$v2) / sqrt(2),
    (basis$v1 - basis$v2) / sqrt(2),
    basis$v3
  )

  rr <- flotsam:::ltsa_ritz_select(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 4
  )

  expect_same_subspace(rr$vectors, target)
  expect_equal(rr$values, c(1, 1), tolerance = 1e-12)
  expect_lt(max(rr$scaled_residuals), 1e-12)
  expect_equal(rr$boundary_gap, 3, tolerance = 1e-12)
})

test_that("rank loss after null projection fails clearly", {
  basis <- selection_test_basis()
  B <- selection_test_matrix(basis)
  candidates <- cbind(basis$u, 2 * basis$u, basis$v1)

  expect_error(
    flotsam:::ltsa_ritz_select(
      B,
      candidates,
      ndim = 2L,
      lambda_max = 3
    ),
    "rank after null projection is 1, less than ndim = 2"
  )
})

test_that("Ritz residual diagnostics use the scaled residual convention", {
  B <- synthetic_ltsa_matrix(c(0, 0.25, 0.5, 2, 8))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  candidates <- dense$vectors[, ord[seq_len(4L)], drop = FALSE]

  rr <- flotsam:::ltsa_ritz_select(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 8
  )

  expect_equal(rr$values, c(0.25, 0.5), tolerance = 1e-12)
  expect_lt(max(rr$absolute_residuals), 1e-12)
  expect_equal(
    rr$scaled_residuals,
    rr$absolute_residuals / 8,
    tolerance = 1e-15
  )
})

test_that("Ritz gap diagnostics report global and local scales", {
  problem <- synthetic_ltsa_problem(c(0, 5e-9, 5e-8, 5e-7, 5e-6, 1))

  rr <- flotsam:::ltsa_ritz_select(
    problem$matrix,
    problem$basis,
    ndim = 2L,
    lambda_max = 1
  )

  expect_equal(unname(rr$boundary_gap), 4.5e-7, tolerance = 1e-8)
  expect_equal(
    unname(rr$global_gap),
    unname(rr$boundary_gap),
    tolerance = 1e-12
  )
  expect_equal(
    unname(rr$boundary_gap_relative),
    unname(rr$global_gap),
    tolerance = 1e-12
  )
  expect_equal(
    unname(rr$zero_tol),
    sqrt(.Machine$double.eps),
    tolerance = 1e-15
  )
  expect_equal(unname(rr$local_gap), 0.9, tolerance = 1e-8)
  expect_gt(abs(rr$local_gap - rr$global_gap), 0.1)
})

test_that("near-zero Ritz counts are reported at multiple thresholds", {
  problem <- synthetic_ltsa_problem(c(0, 5e-9, 5e-8, 5e-7, 5e-6, 1))

  rr <- flotsam:::ltsa_ritz_select(
    problem$matrix,
    problem$basis,
    ndim = 2L,
    lambda_max = 1
  )

  expect_equal(
    rr$near_zero_nonconstant_counts,
    c("1e-08" = 1L, "1e-07" = 2L, "1e-06" = 3L, "1e-05" = 4L)
  )
  expected_ritz_values <- c(5e-9, 5e-8, 5e-7, 5e-6, 1)
  expect_length(rr$reported_ritz_values, length(expected_ritz_values))
  expect_lt(max(abs(rr$reported_ritz_values - expected_ritz_values)), 1e-14)
})

test_that("near-zero trigger calibration documents synthetic spectra", {
  case_a <- calibrated_synthetic_ritz(
    c(0, 1e-10, 2e-10, 1e-3, 1),
    ndim = 3L
  )
  expect_identical(case_a$near_zero_nonconstant_count, 2L)
  expect_true(isTRUE(case_a$partial_near_zero_block))
  expect_true(flotsam:::ltsa_strict_rescue_needed(case_a, ndim = 3L))

  case_b <- calibrated_synthetic_ritz(
    c(0, 1e-8, 1e-4, 1e-3, 1),
    ndim = 2L
  )
  expect_identical(case_b$near_zero_nonconstant_count, 1L)
  expect_true(isTRUE(case_b$partial_near_zero_block))
  expect_warning(
    flotsam:::ltsa_maybe_warn_partial_near_zero_block(
      case_b,
      ndim = 2L,
      strict_rescue = FALSE
    ),
    "missing near-zero LTSA coordinate"
  )

  case_c_unit_scale <- calibrated_synthetic_ritz(
    c(0, 1e-6, 2e-6, 3e-6, 1),
    ndim = 2L
  )
  expect_identical(case_c_unit_scale$near_zero_nonconstant_count, 0L)
  expect_false(isTRUE(case_c_unit_scale$partial_near_zero_block))
  expect_identical(
    unname(case_c_unit_scale$near_zero_nonconstant_counts[["1e-05"]]),
    3L
  )

  case_c_large_scale <- calibrated_synthetic_ritz(
    c(0, 1e-6, 2e-6, 3e-6, 100),
    ndim = 2L
  )
  expect_gt(case_c_large_scale$near_zero_tol, 1e-6)
  expect_lt(case_c_large_scale$near_zero_tol, 2e-6)
  expect_identical(case_c_large_scale$near_zero_nonconstant_count, 1L)
  expect_true(isTRUE(case_c_large_scale$partial_near_zero_block))
})

test_that("multiple exact zeros use ambiguity diagnostics, not partial-block rescue", {
  exact_equal_ndim <- calibrated_synthetic_ritz(
    c(0, 0, 0, 1e-3, 1),
    ndim = 2L
  )
  expect_identical(exact_equal_ndim$near_zero_nonconstant_count, 2L)
  expect_false(isTRUE(exact_equal_ndim$partial_near_zero_block))
  expect_false(flotsam:::ltsa_strict_rescue_needed(
    exact_equal_ndim,
    ndim = 2L
  ))
  expect_length(
    flotsam:::ltsa_spectral_ambiguity_issues(exact_equal_ndim, ndim = 2L),
    0L
  )

  exact_extra <- calibrated_synthetic_ritz(
    c(0, 0, 0, 1e-3, 1),
    ndim = 1L
  )
  expect_identical(exact_extra$near_zero_nonconstant_count, 2L)
  expect_false(isTRUE(exact_extra$partial_near_zero_block))
  expect_false(flotsam:::ltsa_strict_rescue_needed(exact_extra, ndim = 1L))
  issues <- flotsam:::ltsa_spectral_ambiguity_issues(
    exact_extra,
    ndim = 1L
  )
  expect_true(any(grepl(
    "2 near-zero nonconstant modes",
    issues,
    fixed = TRUE
  )))
})

test_that("small Ritz-selected cases agree with dense eigen reference subspaces", {
  B <- synthetic_ltsa_matrix(c(0, 0.1, 0.2, 0.8, 1.5, 3))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  lambda_max <- max(dense$values)
  reference <- dense$vectors[, ord[2:3], drop = FALSE]
  candidates <- dense$vectors[, ord[seq_len(5L)], drop = FALSE]
  # fmt: skip
  candidates <- candidates %*% qr.Q(qr(matrix(c(
    1, 0, 2, 0, -1,
    0, 1, 1, -1, 0,
    1, -1, 0, 1, 2,
    2, 0, -1, 0, 1,
    0, 2, 0, 1, -1
  ), nrow = 5L)))

  rr <- flotsam:::ltsa_ritz_select(
    B,
    candidates,
    ndim = 2L,
    lambda_max = lambda_max
  )

  expect_equal(rr$values, dense$values[ord[2:3]], tolerance = 1e-12)
  expect_same_subspace(rr$vectors, reference, tolerance = 1e-7)
  expect_lt(max(rr$scaled_residuals), 1e-12)
  expect_gt(rr$boundary_gap_relative, 0)
})

test_that("fixed-width default eig_k follows the public rule", {
  expect_identical(flotsam:::ltsa_default_eig_k(ndim = 2L, n = 50L), 12L)
  expect_identical(flotsam:::ltsa_default_eig_k(ndim = 15L, n = 50L), 17L)
  expect_identical(flotsam:::ltsa_default_eig_k(ndim = 2L, n = 8L), 7L)
  expect_identical(flotsam:::ltsa_validate_eig_k(NULL, ndim = 2L, n = 50L), 12L)
})

test_that("fixed-width driver accepts minimum eig_k equal to ndim plus one", {
  problem <- synthetic_ltsa_problem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- fixed_width_provider_factory(problem)

  res <- flotsam:::ltsa_fixed_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 3L
  )

  expect_equal(fixture$calls()$eig_k, 3L)
  expect_equal(res$values, c(0.1, 0.2), tolerance = 1e-12)
  expect_same_subspace(res$vectors, problem$basis[, 2:3], tolerance = 1e-7)
  expect_identical(res$eigen$eig_k, 3L)
  expect_identical(res$eigen$rank, 2L)
  expect_identical(res$eigen$status, "warning")
  expect_true(any(grepl("no spare boundary", res$eigen$messages)))
})

test_that("fixed-width eig_k validation rejects values below ndim plus one", {
  expect_error(
    flotsam:::ltsa_validate_eig_k(2L, ndim = 2L, n = 6L),
    "ndim \\+ 1 <= eig_k < n"
  )
})

test_that("fixed-width eig_k validation rejects values at least n", {
  expect_error(
    flotsam:::ltsa_validate_eig_k(6L, ndim = 2L, n = 6L),
    "ndim \\+ 1 <= eig_k < n"
  )
})

test_that("fixed-width driver calls the provider exactly once", {
  problem <- synthetic_ltsa_problem(c(
    0, 0.1, 0.2, 1, 2, 3, 5, 8,
    13, 21, 34, 55, 89, 144, 233, 377
  ))
  fixture <- fixed_width_provider_factory(problem)

  res <- flotsam:::ltsa_fixed_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider
  )

  expect_equal(nrow(fixture$calls()), 1L)
  expect_equal(fixture$calls()$eig_k, 12L)
  expect_identical(res$eigen$eig_k, 12L)
  expect_identical(res$eigen$status, "ok")
})

test_that("fixed-width diagnostics use compact solver-neutral shape", {
  problem <- synthetic_ltsa_problem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- fixed_width_provider_factory(problem)

  res <- flotsam:::ltsa_fixed_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 4L
  )

  expect_named(
    res$eigen,
    c(
      "method",
      "eig_k",
      "values",
      "ritz_values",
      "residuals",
      "rank",
      "lambda_max",
      "status",
      "messages",
      "backend"
    )
  )
  expect_equal(res$eigen$values, res$values, tolerance = 1e-12)
  expect_equal(res$eigen$residuals, rep(0, 2L), tolerance = 1e-12)
  expect_identical(res$eigen$method, "synthetic")
  expect_identical(res$eigen$backend$name, "synthetic")
  expect_false(any(c(
    "attempts",
    "acceptance",
    "boundary_gap",
    "global_gap",
    "local_gap",
    "zero_tol",
    "near_zero_tol"
  ) %in% names(res$eigen)))
  expect_false(any(c("attempts", "acceptance", "ritz") %in% names(res)))
})

test_that("fixed-width diagnostics mark RSpectra nconv shortfall invalid", {
  problem <- synthetic_ltsa_problem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- fixed_width_provider_factory(
    problem,
    backend = "rspectra",
    convergence_known = TRUE,
    nconv = function(eig_k) eig_k - 1L,
    converged_columns = function(eig_k) eig_k - 1L,
    niter = 11L,
    nops = 101L
  )

  res <- flotsam:::ltsa_fixed_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 4L
  )

  expect_identical(res$eigen$status, "invalid")
  expect_identical(res$eigen$backend$nconv, 3L)
  expect_identical(res$eigen$backend$niter, 11L)
  expect_identical(res$eigen$backend$nops, 101L)
  expect_true(any(grepl("3 / 4", res$eigen$messages, fixed = TRUE)))
})

test_that("fixed-width diagnostics do not invent non-RSpectra convergence", {
  problem <- synthetic_ltsa_problem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- fixed_width_provider_factory(
    problem,
    backend = "irlba",
    convergence_known = FALSE,
    converged_columns = function(eig_k) NA_integer_,
    niter = 7L,
    mprod = 31L
  )

  res <- flotsam:::ltsa_fixed_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 4L
  )

  expect_identical(res$eigen$status, "warning")
  expect_identical(res$eigen$backend$name, "irlba")
  expect_false(res$eigen$backend$convergence_known)
  expect_identical(res$eigen$backend$iter, 7L)
  expect_identical(res$eigen$backend$mprod, 31L)
  expect_false("nconv" %in% names(res$eigen$backend))
  expect_true(any(grepl(
    "native convergence certificate",
    res$eigen$messages,
    fixed = TRUE
  )))
})

test_that("fixed-width driver handles arbitrary normalized-style null vectors", {
  nullvec <- c(1, 2, 3, 2, 1, 4)
  nullvec <- nullvec / sqrt(sum(nullvec * nullvec))
  Z <- cbind(
    c(1, -1, 0, 0, 0, 0),
    c(0, 1, -1, 0, 0, 0),
    c(0, 0, 1, -1, 0, 0),
    c(0, 0, 0, 1, -1, 0),
    c(0, 0, 0, 0, 1, -1)
  )
  Z <- Z - nullvec %*% crossprod(nullvec, Z)
  Q <- qr.Q(qr(cbind(nullvec, Z)))
  basis <- Q[, seq_len(6L), drop = FALSE]
  B <- basis %*% diag(c(0, 0.1, 0.2, 1, 3, 5)) %*% t(basis)
  candidates <- basis[, c(3L, 1L, 4L, 2L, 5L), drop = FALSE]
  # fmt: skip
  candidates <- candidates %*% qr.Q(qr(matrix(c(
    1, 0, 2, 0, -1,
    0, 1, 1, -1, 0,
    1, -1, 0, 1, 2,
    2, 0, -1, 0, 1,
    0, 2, 0, 1, -1
  ), nrow = 5L)))
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    flotsam:::ltsa_candidate_result(
      vectors = candidates[, seq_len(eig_k), drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 5,
      convergence_known = TRUE,
      returned_columns = eig_k,
      converged_columns = eig_k
    )
  }

  res <- flotsam:::ltsa_fixed_ritz_eig(
    B,
    ndim = 2L,
    provider = provider,
    nullvec = nullvec,
    eig_k = 5L
  )

  expect_equal(res$values, c(0.1, 0.2), tolerance = 1e-12)
  expect_same_subspace(res$vectors, basis[, 2:3], tolerance = 1e-7)
  expect_lt(max(res$eigen$residuals), 1e-12)
  expect_identical(res$eigen$status, "ok")
})

test_that("LTSA symmetrization returns a general sparse solver matrix", {
  B <- Matrix::sparseMatrix(
    i = c(1L, 2L, 2L, 3L),
    j = c(1L, 1L, 2L, 3L),
    x = c(2, 4, 6, 8),
    dims = c(3L, 3L)
  )

  B_sym <- flotsam:::symmetrize_ltsa_matrix(B)

  expect_s4_class(B_sym, "dgCMatrix")
  expect_false(methods::is(B_sym, "dsCMatrix"))
  expect_true(Matrix::isSymmetric(B_sym))
  expect_equal(as.matrix(B_sym), 0.5 * (as.matrix(B) + t(as.matrix(B))))
})

test_that("dense LTSA eigensolver fallback reports small residuals", {
  B <- Matrix::Diagonal(x = seq(0, 11))

  res <- flotsam:::rs_eig(B, k = 6L)

  expect_identical(res$backend, "dense_eigen")
  expect_identical(res$eig_k, 6L)
  expect_identical(res$nconv, 6L)
  expect_equal(res$values, seq(0, 5), tolerance = 1e-12)
  expect_lt(max(res$scaled_residuals), 1e-12)
})

test_that("RSpectra path uses shifted largest-algebraic solve with residual metadata", {
  B <- Matrix::Diagonal(x = seq(0, 29))

  res <- flotsam:::rs_eig(
    B,
    k = 6L,
    dense_n = 0L,
    tol = 1e-10,
    maxitr = 5000L
  )

  expect_identical(res$backend, "rspectra")
  expect_identical(res$solve_which, "LA")
  expect_identical(res$eig_k, 6L)
  expect_gte(res$nconv, 6L)
  expect_equal(res$values, seq(0, 5), tolerance = 1e-8)
  expect_lt(max(res$scaled_residuals), 1e-8)
  expect_s4_class(res$matrix, "dgCMatrix")
})

test_that("RSpectra candidate provider returns backend-neutral fields", {
  B <- Matrix::Diagonal(x = seq(0, 29))

  res <- flotsam:::ltsa_rspectra_candidate_provider(
    B,
    eig_k = 6L,
    dense_n = 0L,
    tol = 1e-10,
    maxitr = 5000L
  )

  expect_true(all(
    c(
      "vectors",
      "values",
      "shifted_values",
      "backend",
      "eig_k",
      "matrix",
      "lambda_max",
      "lambda_probe",
      "nconv",
      "niter",
      "nops",
      "mprod",
      "opts",
      "convergence_known",
      "returned_columns",
      "converged_columns"
    ) %in%
      names(res)
  ))
  expect_identical(res$backend, "rspectra")
  expect_identical(res$eig_k, 6L)
  expect_true(res$convergence_known)
  expect_gte(res$nconv, 6L)
  expect_equal(res$values, seq(0, 5), tolerance = 1e-8)
  expect_equal(ncol(res$vectors), 6L)
  expect_s4_class(res$matrix, "dgCMatrix")
})

test_that("RSpectra Ritz wrapper forwards backend tolerance and ncv controls", {
  B <- synthetic_ltsa_matrix(c(
    0, 0.1, 0.2, 1, 3, 5, 8, 13, 21, 34,
    55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181
  ))

  res <- flotsam:::ltsa_rspectra_ritz_eig(
    Matrix::Matrix(B, sparse = TRUE),
    ndim = 2L,
    dense_n = 0L,
    tol = 1e-8,
    ncv = 18L,
    maxitr = 5000L,
    strict_rescue = FALSE
  )

  expect_identical(res$opts$tol, 1e-8)
  expect_identical(res$opts$ncv, 18L)
  expect_identical(res$attempts[[1L]]$backend, "rspectra")
  expect_identical(res$attempts[[1L]]$ncv, 18L)
  expect_equal(res$attempts[[1L]]$eig_k, 8L)
})

test_that("irlba and svdr candidate providers return backend-neutral fields", {
  B <- Matrix::Diagonal(x = seq(0, 29))
  providers <- list(
    irlba = list(
      provider = flotsam:::ltsa_irlba_candidate_provider,
      args = list(dense_n = 0L, tol = 1e-10, maxit = 1000L)
    ),
    svdr = list(
      provider = flotsam:::ltsa_svdr_candidate_provider,
      args = list(dense_n = 0L, tol = 1e-10, it = 1000L)
    )
  )

  for (backend in names(providers)) {
    set.seed(11)
    res <- do.call(
      providers[[backend]]$provider,
      c(list(B = B, eig_k = 6L), providers[[backend]]$args)
    )

    expect_identical(res$backend, backend)
    expect_identical(res$eig_k, 6L)
    expect_false(res$convergence_known)
    expect_true(is.na(res$nconv))
    expect_true(is.na(res$converged_columns))
    expect_true(is.finite(res$mprod))
    expect_equal(res$values, seq(0, 5), tolerance = 1e-8)
    expect_equal(ncol(res$vectors), 6L)
    expect_s4_class(res$matrix, "dgCMatrix")
  }
})

test_that("generic adaptive Ritz driver consumes the RSpectra provider", {
  B <- synthetic_ltsa_matrix(c(0, 0.1, 0.2, 1, 3, 5, 8, 13))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  reference <- dense$vectors[, ord[2:3], drop = FALSE]

  res <- flotsam:::ltsa_adaptive_ritz_eig(
    Matrix::Matrix(B, sparse = TRUE),
    ndim = 2L,
    provider = flotsam:::ltsa_rspectra_candidate_provider,
    provider_args = list(dense_n = 0L, tol = 1e-10, maxitr = 5000L),
    strict_rescue_arg_mapper = flotsam:::ltsa_strict_rescue_args,
    strict_rescue_controls = list(
      strict_rescue_tol = 1e-10,
      strict_rescue_maxitr = 5000L
    )
  )

  expect_identical(res$backend, "rspectra")
  expect_true(!is.null(res$candidate_vectors))
  expect_equal(res$values, dense$values[ord[2:3]], tolerance = 1e-8)
  expect_same_subspace(res$vectors, reference, tolerance = 1e-7)
  expect_lt(max(res$scaled_residuals), 1e-8)
})

test_that("adaptive irlba and svdr paths agree with dense reference subspaces", {
  # fmt: skip
  B <- synthetic_ltsa_matrix(c(
    0, 0.1, 0.2, 1, 3, 5, 8, 13, 21, 34,
    55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181
  ))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  reference <- dense$vectors[, ord[2:3], drop = FALSE]
  backends <- list(
    irlba = list(
      eig = flotsam:::ltsa_irlba_ritz_eig,
      args = list(dense_n = 0L, tol = 1e-10, maxit = 5000L)
    ),
    svdr = list(
      eig = flotsam:::ltsa_svdr_ritz_eig,
      args = list(dense_n = 0L, tol = 1e-10, it = 5000L, extra = 12L)
    )
  )

  for (backend in names(backends)) {
    set.seed(42)
    res <- do.call(
      backends[[backend]]$eig,
      c(
        list(
          B = Matrix::Matrix(B, sparse = TRUE),
          ndim = 2L,
          strict_rescue = FALSE
        ),
        backends[[backend]]$args
      )
    )

    expect_identical(res$backend, backend)
    expect_equal(res$values, dense$values[ord[2:3]], tolerance = 1e-6)
    expect_same_subspace(res$vectors, reference, tolerance = 1e-5)
    expect_lt(max(res$scaled_residuals), 1e-6)
  }
})

test_that("adaptive Ritz driver handles rotated normalized null vectors", {
  nullvec <- c(1, 2, 3, 2, 1, 4)
  nullvec <- nullvec / sqrt(sum(nullvec * nullvec))
  Z <- cbind(
    c(1, -1, 0, 0, 0, 0),
    c(0, 1, -1, 0, 0, 0),
    c(0, 0, 1, -1, 0, 0),
    c(0, 0, 0, 1, -1, 0),
    c(0, 0, 0, 0, 1, -1)
  )
  Z <- Z - nullvec %*% crossprod(nullvec, Z)
  Q <- qr.Q(qr(cbind(nullvec, Z)))
  basis <- Q[, seq_len(6L), drop = FALSE]
  B <- basis %*% diag(c(0, 0.1, 0.2, 1, 3, 5)) %*% t(basis)
  candidates <- basis[, c(3L, 1L, 4L, 2L, 5L), drop = FALSE]
  # fmt: skip
  candidates <- candidates %*% qr.Q(qr(matrix(c(
    1, 0, 2, 0, -1,
    0, 1, 1, -1, 0,
    1, -1, 0, 1, 2,
    2, 0, -1, 0, 1,
    0, 2, 0, 1, -1
  ), nrow = 5L)))
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    flotsam:::ltsa_candidate_result(
      vectors = candidates[, seq_len(eig_k), drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 5,
      convergence_known = FALSE,
      returned_columns = eig_k
    )
  }

  res <- flotsam:::ltsa_adaptive_ritz_eig(
    B,
    ndim = 2L,
    provider = provider,
    nullvec = nullvec,
    strict_rescue = FALSE
  )

  expect_equal(res$values, c(0.1, 0.2), tolerance = 1e-12)
  expect_same_subspace(res$vectors, basis[, 2:3], tolerance = 1e-7)
  expect_lt(max(res$scaled_residuals), 1e-12)
  expect_identical(res$backend, "synthetic")
})

test_that("adaptive attempt summaries include backend work and coverage metadata", {
  problem <- synthetic_ltsa_problem(c(0, 0.1, 0.2, 1, 2, 3, 4, 5, 6, 7))
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, seq_len(eig_k), drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 7,
      nconv = eig_k,
      niter = 3L,
      nops = 101L,
      mprod = 202L,
      opts = list(ncv = 13L),
      convergence_known = TRUE,
      returned_columns = eig_k,
      converged_columns = eig_k
    )
  }

  res <- flotsam:::ltsa_adaptive_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = provider,
    strict_rescue = FALSE
  )
  attempt <- res$attempts[[1L]]

  expect_identical(attempt$backend, "synthetic")
  expect_identical(attempt$niter, 3L)
  expect_identical(attempt$nops, 101L)
  expect_identical(attempt$mprod, 202L)
  expect_identical(attempt$ncv, 13L)
  expect_identical(attempt$returned_columns, 8L)
  expect_identical(attempt$converged_columns, 8L)
  expect_true(attempt$convergence_known)
  expect_true(is.finite(attempt$candidate_elapsed))
  expect_true(attempt$candidate_elapsed >= 0)
  expect_true(is.finite(attempt$ritz_elapsed))
  expect_true(attempt$ritz_elapsed >= 0)
  expect_null(attempt$.candidate_vectors)
  expect_equal(attempt$reference_projection_norms, c(1, 1), tolerance = 1e-12)
  expect_equal(attempt$reference_projection_min_norm, 1, tolerance = 1e-12)
  expect_identical(attempt$reference_candidate_space_rank, 7L)
  expect_identical(attempt$reference_space_rank, 2L)
  expect_equal(attempt$reference_overlap_singular_values, c(1, 1), tolerance = 1e-12)
  expect_equal(attempt$reference_overlap_min_singular_value, 1, tolerance = 1e-12)
})

test_that("candidate reference projection pads missing reference directions", {
  basis <- selection_test_basis()
  candidates <- cbind(basis$u, basis$v1)
  reference <- cbind(basis$v1, basis$v2)

  coverage <- flotsam:::ltsa_candidate_reference_projection(
    candidate_vectors = candidates,
    reference_vectors = reference,
    nullvec = basis$u
  )

  expect_identical(coverage$reference_candidate_space_rank, 1L)
  expect_identical(coverage$reference_space_rank, 2L)
  expect_equal(coverage$reference_projection_norms, c(1, 0), tolerance = 1e-12)
  expect_equal(coverage$reference_overlap_singular_values, c(1, 0), tolerance = 1e-12)
  expect_equal(coverage$reference_overlap_min_singular_value, 0, tolerance = 1e-12)
})

test_that("adaptive Ritz driver can retain attempts and use a fixed reference", {
  problem <- synthetic_ltsa_problem(c(0, 0.1, 0.2, 1, 2, 3, 4, 5, 6, 7))
  fixed_reference <- problem$basis[, 3:4, drop = FALSE]
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, seq_len(eig_k), drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 7,
      convergence_known = TRUE,
      returned_columns = eig_k,
      converged_columns = eig_k
    )
  }

  res <- flotsam:::ltsa_adaptive_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = provider,
    strict_rescue = FALSE,
    attempt_reference_vectors = fixed_reference,
    retain_attempt_candidate_spaces = TRUE
  )
  attempt <- res$attempts[[1L]]

  expect_false(is.null(attempt$.candidate_vectors))
  expect_equal(attempt$reference_overlap_singular_values, c(1, 1), tolerance = 1e-12)
  expect_equal(attempt$reference_projection_norms, c(1, 1), tolerance = 1e-12)
})

test_that("adaptive RSpectra path returns Ritz-polished vectors", {
  B <- synthetic_ltsa_matrix(c(0, 0.1, 0.2, 1, 3, 5, 8, 13))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  reference <- dense$vectors[, ord[2:3], drop = FALSE]

  res <- flotsam:::ltsa_rspectra_ritz_eig(
    Matrix::Matrix(B, sparse = TRUE),
    ndim = 2L,
    dense_n = 0L,
    tol = 1e-10,
    maxitr = 5000L
  )

  expect_identical(res$backend, "rspectra")
  expect_true(!is.null(res$candidate_vectors))
  expect_equal(res$values, dense$values[ord[2:3]], tolerance = 1e-8)
  expect_same_subspace(res$vectors, reference, tolerance = 1e-7)
  expect_lt(max(res$scaled_residuals), 1e-8)
  expect_gte(res$rank_after_null, 3L)
  expect_true(is.finite(res$boundary_gap_relative))
})

test_that("strict rescue trigger is limited to partial near-zero selected blocks", {
  fake_result <- function(near_zero_count, rank_ok = TRUE, resid_ok = TRUE) {
    list(
      values = c(1e-9, 8e-7),
      near_zero_nonconstant_count = near_zero_count,
      acceptance = list(rank_ok = rank_ok, resid_ok = resid_ok)
    )
  }

  expect_true(flotsam:::ltsa_partial_near_zero_block(fake_result(1L), ndim = 2L))
  expect_false(flotsam:::ltsa_partial_near_zero_block(fake_result(0L), ndim = 2L))
  expect_false(flotsam:::ltsa_partial_near_zero_block(fake_result(2L), ndim = 2L))
  expect_true(flotsam:::ltsa_strict_rescue_needed(fake_result(1L), ndim = 2L))
  expect_false(flotsam:::ltsa_strict_rescue_needed(fake_result(0L), ndim = 2L))
  expect_false(flotsam:::ltsa_strict_rescue_needed(fake_result(2L), ndim = 2L))
  expect_false(flotsam:::ltsa_strict_rescue_needed(
    fake_result(1L, resid_ok = FALSE),
    ndim = 2L
  ))
  expect_false(flotsam:::ltsa_strict_rescue_needed(
    fake_result(1L, rank_ok = FALSE),
    ndim = 2L
  ))
  expect_false(flotsam:::ltsa_partial_near_zero_block(
    fake_result(1L, resid_ok = FALSE),
    ndim = 2L
  ))
})

test_that("partial near-zero helpers mark and warn only no-rescue partial blocks", {
  fake_result <- function(near_zero_count, rank_ok = TRUE, resid_ok = TRUE) {
    list(
      values = c(1e-8, 1e-4),
      near_zero_nonconstant_count = near_zero_count,
      acceptance = list(rank_ok = rank_ok, resid_ok = resid_ok)
    )
  }

  marked <- flotsam:::ltsa_mark_partial_near_zero_block(
    fake_result(1L),
    ndim = 2L
  )
  expect_true(isTRUE(marked$partial_near_zero_block))
  expect_true(isTRUE(marked$acceptance$partial_near_zero_block))

  expect_warning(
    flotsam:::ltsa_maybe_warn_partial_near_zero_block(
      marked,
      ndim = 2L,
      strict_rescue = FALSE
    ),
    "missing near-zero LTSA coordinate"
  )
  expect_warning(
    flotsam:::ltsa_maybe_warn_partial_near_zero_block(
      marked,
      ndim = 2L,
      strict_rescue = TRUE
    ),
    NA
  )

  unmarked <- flotsam:::ltsa_mark_partial_near_zero_block(
    fake_result(1L, resid_ok = FALSE),
    ndim = 2L
  )
  expect_false(isTRUE(unmarked$partial_near_zero_block))
  expect_warning(
    flotsam:::ltsa_maybe_warn_partial_near_zero_block(
      unmarked,
      ndim = 2L,
      strict_rescue = FALSE
    ),
    NA
  )
})

test_that("strict rescue result ranking prefers lower selected block energy", {
  fake_result <- function(
    values,
    near_zero_count,
    rank_ok = TRUE,
    resid_ok = TRUE
  ) {
    list(
      values = values,
      near_zero_nonconstant_count = near_zero_count,
      acceptance = list(rank_ok = rank_ok, resid_ok = resid_ok)
    )
  }
  default <- fake_result(c(2e-9, 8e-7), near_zero_count = 1L)
  strict <- fake_result(c(4e-10, 2e-9), near_zero_count = 2L)
  bad_residual <- fake_result(
    c(4e-10, 2e-9),
    near_zero_count = 2L,
    resid_ok = FALSE
  )

  expect_true(flotsam:::ltsa_energy_better(strict, default, ndim = 2L))
  expect_identical(
    flotsam:::ltsa_rescue_candidate(strict, default, ndim = 2L),
    strict
  )
  expect_identical(
    flotsam:::ltsa_rescue_candidate(bad_residual, default, ndim = 2L),
    default
  )
})

test_that("energy rescue predicate handles synthetic acceptance fixtures", {
  cases <- list(
    list(
      name = "accept_better_lambda_and_trace",
      ndim = 2L,
      incumbent_values = c(1e-10, 7e-5),
      candidate_values = c(1e-10, 2e-10),
      incumbent_scaled_residuals = c(2e-17, 2e-17),
      candidate_scaled_residuals = c(2e-17, 2e-17),
      expected_accept = TRUE,
      expected_decision = "accept_near_zero_lambda_ndim_shortcut",
      expected_warning = ""
    ),
    list(
      name = "reject_similar_lambda_worse_trace",
      ndim = 3L,
      incumbent_values = c(1e-10, 2e-10, 6e-5),
      candidate_values = c(2e-5, 2e-5, 6e-5),
      incumbent_scaled_residuals = c(1e-17, 1e-17, 1e-17),
      candidate_scaled_residuals = c(1e-17, 1e-17, 1e-17),
      expected_accept = FALSE,
      expected_decision = "reject_trace_worse_beyond_residual_allowance",
      expected_warning = ""
    ),
    list(
      name = "reject_worse_lambda_beyond_tolerance",
      ndim = 2L,
      incumbent_values = c(1e-10, 7e-5),
      candidate_values = c(7.2e-5, 7.3e-5),
      incumbent_scaled_residuals = c(1e-17, 1e-17),
      candidate_scaled_residuals = c(1e-17, 1e-17),
      expected_accept = FALSE,
      expected_decision = "reject_lambda_ndim_worse_beyond_residual_allowance",
      expected_warning = ""
    ),
    list(
      name = "accept_energy_tie_nonpartial_replacement",
      ndim = 2L,
      incumbent_values = c(1e-10, 7e-5),
      candidate_values = c(2e-8, 7.00005e-5),
      incumbent_scaled_residuals = c(2e-7, 2e-7),
      candidate_scaled_residuals = c(2e-7, 2e-7),
      expected_accept = TRUE,
      expected_decision = "accept_energy_tie_nonpartial_replacement",
      expected_warning = ""
    ),
    list(
      name = "tiny_negative_clamped_no_spurious_trace_win",
      ndim = 2L,
      incumbent_values = c(1e-10, 5e-5),
      candidate_values = c(-1e-8, 5.0005e-5),
      incumbent_scaled_residuals = c(1e-17, 1e-17),
      candidate_scaled_residuals = c(1e-17, 1e-17),
      expected_accept = FALSE,
      expected_decision = "reject_ineligible;candidate_partial",
      expected_warning = ""
    ),
    list(
      name = "reject_material_negative_below_zero_tol",
      ndim = 2L,
      incumbent_values = c(1e-10, 7e-5),
      candidate_values = c(-2.1e-8, -2.0e-8),
      incumbent_scaled_residuals = c(1e-17, 1e-17),
      candidate_scaled_residuals = c(1e-17, 1e-17),
      expected_accept = FALSE,
      expected_decision = "reject_ineligible;candidate_material_negative",
      expected_warning = "selected eigenvalue below -zero_tol"
    )
  )

  for (case in cases) {
    incumbent <- synthetic_energy_rescue_result(
      values = case$incumbent_values,
      scaled_residuals = case$incumbent_scaled_residuals,
      ndim = case$ndim
    )
    candidate <- synthetic_energy_rescue_result(
      values = case$candidate_values,
      scaled_residuals = case$candidate_scaled_residuals,
      ndim = case$ndim
    )
    decision <- flotsam:::ltsa_energy_rescue_accept_replacement(
      incumbent = incumbent,
      candidate = candidate,
      ndim = case$ndim
    )

    expect_identical(decision$accept, case$expected_accept, info = case$name)
    expect_identical(
      decision$decision,
      case$expected_decision,
      info = case$name
    )
    expect_identical(decision$warning, case$expected_warning, info = case$name)
  }

  tiny_negative <- synthetic_energy_rescue_result(
    values = c(-1e-8, 5.0005e-5),
    scaled_residuals = c(1e-17, 1e-17),
    ndim = 2L
  )
  key <- flotsam:::ltsa_energy_rescue_key(tiny_negative, ndim = 2L)
  expect_identical(unname(key$clamped_values[[1L]]), 0)
  expect_true(isTRUE(key$partial))
  expect_false(isTRUE(key$material_negative))
})

test_that("strict rescue arguments tighten tolerance and preserve stricter user opts", {
  expect_equal(
    flotsam:::ltsa_strict_rescue_args(
      list(),
      strict_rescue_tol = 1e-10,
      strict_rescue_maxitr = 5000L
    ),
    list(tol = 1e-10, maxitr = 5000L)
  )
  expect_equal(
    flotsam:::ltsa_strict_rescue_args(
      list(tol = 1e-12, maxitr = 8000L),
      strict_rescue_tol = 1e-10,
      strict_rescue_maxitr = 5000L
    ),
    list(tol = 1e-12, maxitr = 8000L)
  )
  expect_equal(
    flotsam:::ltsa_strict_rescue_args(
      list(tol = 1e-6, maxitr = 100L),
      strict_rescue_tol = 1e-10,
      strict_rescue_maxitr = 5000L
    ),
    list(tol = 1e-10, maxitr = 5000L)
  )
})

test_that("strict rescue candidate count includes expanded diagnostic attempts", {
  selected <- list(
    eig_k = 8L,
    acceptance = list(diagnostic_final_eig_k = 18L),
    attempts = list(
      list(eig_k = 8L),
      list(eig_k = 12L),
      list(eig_k = 18L)
    )
  )

  expect_equal(
    flotsam:::ltsa_strict_rescue_eig_k(
      selected,
      ndim = 2L,
      n = 100L,
      strict_rescue_extra = 5L
    ),
    18L
  )

  selected$acceptance$diagnostic_final_eig_k <- NULL
  selected$attempts <- list(list(eig_k = 8L), list(eig_k = 12L))
  expect_equal(
    flotsam:::ltsa_strict_rescue_eig_k(
      selected,
      ndim = 2L,
      n = 100L,
      strict_rescue_extra = 5L
    ),
    12L
  )

  selected$eig_k <- 30L
  expect_equal(
    flotsam:::ltsa_strict_rescue_eig_k(
      selected,
      ndim = 2L,
      n = 20L,
      strict_rescue_extra = 5L
    ),
    19L
  )
})

test_that("adaptive Ritz driver keeps strict rescue enabled by default", {
  # fmt: skip
  problem <- synthetic_ltsa_problem(c(
    0, 1e-7, 2e-7, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3,
    2e-3, 5e-3, 1e-2, 2e-2
  ))
  calls <- data.frame(eig_k = integer(), strict = logical())
  provider <- function(
    B,
    eig_k,
    lambda_max = NULL,
    verbose = FALSE,
    strict = FALSE
  ) {
    calls <<- rbind(calls, data.frame(eig_k = eig_k, strict = strict))
    cols <- if (strict) {
      seq_len(eig_k)
    } else {
      c(1L, 2L, 4L:9L)
    }
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, cols, drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 100,
      convergence_known = FALSE,
      returned_columns = length(cols)
    )
  }

  res <- NULL
  expect_warning(
    res <- flotsam:::ltsa_adaptive_ritz_eig(
      problem$matrix,
      ndim = 2L,
      provider = provider,
      gap_expansion_steps = 0L,
      strict_rescue_arg_mapper = function(provider_args, ...) {
        c(provider_args, list(strict = TRUE))
      }
    ),
    NA
  )

  expect_equal(calls$eig_k, c(8L, 8L))
  expect_equal(calls$strict, c(FALSE, TRUE))
  expect_true(isTRUE(res$attempts[[1L]]$partial_near_zero_block))
  expect_false(isTRUE(res$attempts[[2L]]$partial_near_zero_block))
  expect_false(isTRUE(res$partial_near_zero_block))
  expect_false(isTRUE(res$acceptance$partial_near_zero_block))
  expect_true(isTRUE(res$acceptance$strict_rescue_used))
  expect_identical(res$acceptance$return_reason, "strict_rescue")
  expect_lt(max(abs(res$values - c(1e-7, 2e-7))), 1e-14)
})

test_that("width-first rescue adaptive widths cover ndim and small-n caps", {
  adaptive_widths <- function(ndim, n, count) {
    max_k <- flotsam:::ltsa_max_ritz_candidate_k(ndim, n)
    widths <- flotsam:::ltsa_initial_ritz_candidate_k(ndim, n)
    while (length(widths) < count && tail(widths, 1L) < max_k) {
      widths <- c(
        widths,
        flotsam:::ltsa_next_ritz_candidate_k(tail(widths, 1L), max_k)
      )
    }
    widths
  }

  expect_equal(adaptive_widths(2L, 96L, 3L), c(8L, 12L, 18L))
  expect_equal(adaptive_widths(3L, 96L, 3L), c(9L, 14L, 21L))
  expect_equal(adaptive_widths(5L, 96L, 3L), c(11L, 17L, 26L))
  expect_equal(adaptive_widths(5L, 13L, 3L), c(11L, 12L))
})

test_that("width-first rescue accepts ordinary widening after a partial block", {
  out <- NULL
  expect_warning(
    out <- run_width_rescue_case(
      ordinary_mode = function(eig_k) {
        if (eig_k < 12L) "partial" else "complete"
      }
    ),
    NA
  )

  calls <- out$calls
  res <- out$result
  expect_equal(calls$eig_k, c(8L, 12L))
  expect_equal(calls$strict, c(FALSE, FALSE))
  expect_equal(
    vapply(res$attempts, `[[`, integer(1), "eig_k"),
    c(8L, 12L)
  )
  expect_identical(res$acceptance$return_reason, "width_first_ordinary_rescue")
  expect_true(isTRUE(res$acceptance$width_first_rescue_ordinary))
  expect_false(isTRUE(res$partial_near_zero_block))
  expect_lt(max(abs(res$values - c(1e-10, 2e-10))), 1e-14)
})

test_that("width-first rescue invokes strict fallback after ordinary exhaustion", {
  out <- NULL
  expect_warning(
    out <- run_width_rescue_case(
      ordinary_mode = function(eig_k) "partial",
      strict_mode = function(eig_k) "complete"
    ),
    NA
  )

  calls <- out$calls
  ordinary_calls <- calls$eig_k[!calls$strict]
  strict_calls <- calls$eig_k[calls$strict]
  res <- out$result
  expect_equal(ordinary_calls, c(8L, 12L, 18L))
  expect_equal(anyDuplicated(ordinary_calls), 0L)
  expect_equal(strict_calls, 12L)
  expect_identical(res$acceptance$return_reason, "width_first_strict_rescue")
  expect_true(isTRUE(res$acceptance$strict_rescue_used))
  expect_identical(
    res$acceptance$strict_rescue_stage,
    "previous_ordinary_width"
  )
  expect_identical(res$acceptance$strict_rescue_width_cap, 18L)
  expect_message(
    flotsam:::ltsa_maybe_message_width_first_rescue_decision(
      res,
      ndim = 2L,
      verbose = TRUE
    ),
    "ordinary widths 8, 12, 18 remained partial; strict width 12 accepted"
  )
})

test_that("width-first staged strict fallback tries widest after previous remains partial", {
  out <- NULL
  expect_warning(
    out <- run_width_rescue_case(
      ordinary_mode = function(eig_k) "partial",
      strict_mode = function(eig_k) {
        if (eig_k == 12L) "partial" else "complete"
      }
    ),
    NA
  )

  calls <- out$calls
  res <- out$result
  expect_equal(calls$eig_k, c(8L, 12L, 18L, 12L, 18L))
  expect_equal(calls$strict, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(
    res$acceptance$strict_rescue_stage,
    "widest_ordinary_width"
  )
  strict_attempts <- Filter(function(x) isTRUE(x$strict_rescue), res$attempts)
  expect_equal(
    vapply(strict_attempts, `[[`, character(1), "strict_rescue_stage"),
    c("previous_ordinary_width", "widest_ordinary_width")
  )
  expect_equal(
    vapply(strict_attempts, `[[`, integer(1), "strict_rescue_width_cap"),
    c(18L, 18L)
  )
  expect_lt(max(abs(res$values - c(1e-10, 2e-10))), 1e-14)
})

test_that("width-first rescue rejects energy-worse ordinary nonpartial candidates", {
  out <- NULL
  expect_warning(
    out <- run_width_rescue_case(
      ordinary_mode = function(eig_k) {
        if (eig_k == 8L) {
          "partial"
        } else if (eig_k == 12L) {
          "high"
        } else {
          "complete"
        }
      }
    ),
    NA
  )

  calls <- out$calls
  res <- out$result
  expect_equal(calls$eig_k, c(8L, 12L, 18L))
  expect_equal(calls$strict, c(FALSE, FALSE, FALSE))
  expect_identical(
    res$attempts[[2L]]$width_first_rescue_decision,
    "reject_lambda_ndim_worse_beyond_residual_allowance"
  )
  expect_false(isTRUE(res$attempts[[2L]]$width_first_rescue_accept))
  expect_identical(res$acceptance$return_reason, "width_first_ordinary_rescue")
  expect_lt(max(abs(res$values - c(1e-10, 2e-10))), 1e-14)
})

test_that("width-first rescue honors small-n cap without duplicate ordinary attempts", {
  out <- NULL
  expect_warning(
    out <- run_width_rescue_case(
      ndim = 5L,
      n = 13L,
      ordinary_mode = function(eig_k) "partial",
      strict_mode = function(eig_k) "complete"
    ),
    NA
  )

  calls <- out$calls
  ordinary_calls <- calls$eig_k[!calls$strict]
  strict_calls <- calls$eig_k[calls$strict]
  expect_equal(ordinary_calls, c(11L, 12L))
  expect_equal(anyDuplicated(ordinary_calls), 0L)
  expect_equal(strict_calls, 11L)
  expect_identical(
    out$result$acceptance$strict_rescue_stage,
    "previous_ordinary_width"
  )
  expect_identical(out$result$acceptance$strict_rescue_width_cap, 12L)
})

test_that("width-first rescue returns warned partial incumbent after capped failures", {
  out <- NULL
  expect_warning(
    expect_warning(
      out <- run_width_rescue_case(
        ordinary_mode = function(eig_k) "partial",
        strict_mode = function(eig_k) "partial"
      ),
      "width-first rescue exhausted capped ordinary and strict attempts"
    ),
    "ambiguous low-energy eigenspace"
  )

  calls <- out$calls
  res <- out$result
  expect_equal(calls$eig_k, c(8L, 12L, 18L, 12L, 18L))
  expect_equal(calls$strict, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_true(isTRUE(res$partial_near_zero_block))
  expect_true(isTRUE(res$acceptance$width_first_rescue_unresolved))
  expect_false(isTRUE(res$acceptance$strict_rescue_used))
  expect_identical(
    res$acceptance$return_reason,
    "width_first_rescue_unresolved_partial"
  )
})

test_that("disabled strict rescue warns on returned partial near-zero blocks", {
  # fmt: skip
  problem <- synthetic_ltsa_problem(c(
    0, 1e-7, 2e-7, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3,
    2e-3, 5e-3, 1e-2, 2e-2
  ))
  calls <- integer()
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    calls <<- c(calls, eig_k)
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, c(1L, 2L, 4L:9L), drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 100,
      convergence_known = FALSE,
      returned_columns = eig_k
    )
  }

  res <- NULL
  expect_warning(
    res <- flotsam:::ltsa_adaptive_ritz_eig(
      problem$matrix,
      ndim = 2L,
      provider = provider,
      gap_expansion_steps = 0L,
      strict_rescue = FALSE
    ),
    "missing near-zero LTSA coordinate"
  )

  expect_equal(calls, 8L)
  expect_true(res$acceptance$rank_ok)
  expect_true(res$acceptance$resid_ok)
  expect_identical(res$near_zero_nonconstant_count, 1L)
  expect_true(isTRUE(res$partial_near_zero_block))
  expect_true(isTRUE(res$acceptance$partial_near_zero_block))
  expect_true(isTRUE(res$attempts[[1L]]$partial_near_zero_block))
  expect_identical(
    res$acceptance$return_reason,
    "weak_lowest_energy_residual_rank_good"
  )
})

test_that("partial near-zero warning preserves gap-ok return reason", {
  # fmt: skip
  problem <- synthetic_ltsa_problem(c(
    0, 1e-7, 2e-7, 1, 2, 3, 4, 5, 6, 7, 8, 9
  ))
  calls <- integer()
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    calls <<- c(calls, eig_k)
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, c(1L, 2L, 4L:9L), drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 100,
      convergence_known = FALSE,
      returned_columns = eig_k
    )
  }

  res <- NULL
  expect_warning(
    res <- flotsam:::ltsa_adaptive_ritz_eig(
      problem$matrix,
      ndim = 2L,
      provider = provider,
      strict_rescue = FALSE
    ),
    "missing near-zero LTSA coordinate"
  )

  expect_equal(calls, 8L)
  expect_true(res$acceptance$gap_ok)
  expect_true(isTRUE(res$partial_near_zero_block))
  expect_identical(res$acceptance$return_reason, "residual_rank_gap_ok")
})

test_that("high-energy partial near-zero blocks do not trigger rescue", {
  # fmt: skip
  problem <- synthetic_ltsa_problem(c(
    0, 1e-7, 2e-7, 1, 2, 3, 4, 5, 6, 7, 8, 9
  ))
  calls <- data.frame(eig_k = integer(), strict = logical())
  provider <- function(
    B,
    eig_k,
    lambda_max = NULL,
    verbose = FALSE,
    strict = FALSE
  ) {
    calls <<- rbind(calls, data.frame(eig_k = eig_k, strict = strict))
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, c(1L, 2L, 4L:9L), drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 100,
      convergence_known = FALSE,
      returned_columns = eig_k
    )
  }

  res <- NULL
  expect_warning(
    res <- flotsam:::ltsa_adaptive_ritz_eig(
      problem$matrix,
      ndim = 2L,
      provider = provider,
      strict_rescue_arg_mapper = function(provider_args, ...) {
        c(provider_args, list(strict = TRUE))
      }
    ),
    NA
  )

  expect_equal(calls$eig_k, 8L)
  expect_false(any(calls$strict))
  expect_true(res$acceptance$gap_ok)
  expect_true(isTRUE(res$partial_near_zero_block))
  expect_false(isTRUE(res$acceptance$strict_rescue_used))
  expect_identical(res$acceptance$return_reason, "residual_rank_gap_ok")
})

test_that("extra near-zero nonconstant modes warn about spectral ambiguity", {
  B <- synthetic_ltsa_matrix(c(0, 1e-9, 2e-9, 3e-9, 4e-9, 1))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    flotsam:::ltsa_candidate_result(
      vectors = dense$vectors[, ord[seq_len(eig_k)], drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 1,
      convergence_known = FALSE,
      returned_columns = eig_k
    )
  }

  res <- NULL
  expect_warning(
    res <- flotsam:::ltsa_adaptive_ritz_eig(
      B,
      ndim = 2L,
      provider = provider,
      strict_rescue = TRUE
    ),
    "near-zero nonconstant modes"
  )

  expect_true(res$acceptance$spectral_ambiguity_warning)
  expect_gt(res$near_zero_nonconstant_count, 2L)
  expect_true(any(grepl(
    "near-zero nonconstant modes",
    res$acceptance$spectral_ambiguity_issues
  )))
})

test_that("adaptive weak-gap expansion keeps lower-energy candidates", {
  # fmt: skip
  problem <- synthetic_ltsa_problem(c(
    0, 1e-7, 2e-7, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3,
    2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.4,
    0.8, 1, 2, 4, 8, 16, 32, 100
  ))
  calls <- integer()
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    calls <<- c(calls, eig_k)
    cols <- if (eig_k == 8L) {
      c(1L, 2L, 4L:9L)
    } else {
      seq_len(eig_k)
    }
    flotsam:::ltsa_candidate_result(
      vectors = problem$basis[, cols, drop = FALSE],
      backend = "synthetic",
      eig_k = eig_k,
      matrix = B,
      lambda_max = 100,
      convergence_known = FALSE,
      returned_columns = length(cols)
    )
  }

  res <- NULL
  expect_warning(
    res <- flotsam:::ltsa_adaptive_ritz_eig(
      problem$matrix,
      ndim = 2L,
      provider = provider,
      strict_rescue = TRUE
    ),
    NA
  )

  expect_equal(calls, c(8L, 12L))
  expect_equal(res$eig_k, 12L)
  expect_lt(max(abs(res$values - c(1e-7, 2e-7))), 1e-14)
  expect_equal(
    res$attempts[[1L]]$reference_projection_norms,
    c(1, 0),
    tolerance = 1e-8
  )
  expect_equal(
    res$attempts[[2L]]$reference_projection_norms,
    c(1, 1),
    tolerance = 1e-12
  )
  expect_lt(res$attempts[[1L]]$reference_projection_min_norm, 1e-8)
  expect_equal(
    res$attempts[[2L]]$reference_projection_min_norm,
    1,
    tolerance = 1e-12
  )
  expect_false(isTRUE(res$acceptance$strict_rescue_used))
  expect_false(isTRUE(res$acceptance$spectral_ambiguity_warning))
  expect_identical(
    res$acceptance$return_reason,
    "weak_lowest_energy_residual_rank_good"
  )
  expect_equal(res$acceptance$diagnostic_final_eig_k, 12L)
})

test_that("weak-gap-only wrapper default confirms once", {
  # fmt: skip
  B <- synthetic_ltsa_matrix(c(
    0, 0.1, 0.2, 0.200001, 1, 2, 3, 5, 8, 13,
    21, 34, 55, 89, 144, 233, 377, 610, 987, 1597
  ))

  messages <- character()
  res <- NULL
  expect_warning(
    res <- withCallingHandlers(
      flotsam:::ltsa_rspectra_ritz_eig(
        Matrix::Matrix(B, sparse = TRUE),
        ndim = 2L,
        dense_n = 0L,
        tol = 1e-10,
        maxitr = 5000L,
        verbose = TRUE
      ),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    "boundary gap after ndim is weak"
  )
  expect_true(any(grepl(
    paste0(
      "LTSA Ritz boundary gap remains below tolerance ",
      "after diagnostic expansion to 12 candidates; ",
      "returning residual-good, rank-good Ritz vectors"
    ),
    messages
  )))

  requested <- vapply(res$attempts, `[[`, integer(1), "eig_k")
  expect_equal(requested, c(8L, 12L))
  expect_equal(res$acceptance$diagnostic_final_eig_k, 12L)
  expect_true(res$acceptance$resid_ok)
  expect_true(res$acceptance$rank_ok)
  expect_false(res$acceptance$gap_ok)
  expect_true(res$acceptance$spectral_ambiguity_warning)
  expect_identical(
    res$acceptance$return_reason,
    "weak_lowest_energy_residual_rank_good"
  )
})

test_that("explicit two-step weak-gap expansion stops before max_extra", {
  # fmt: skip
  B <- synthetic_ltsa_matrix(c(
    0, 0.1, 0.2, 0.200001, 1, 2, 3, 5, 8, 13,
    21, 34, 55, 89, 144, 233, 377, 610, 987, 1597
  ))

  messages <- character()
  res <- NULL
  expect_warning(
    res <- withCallingHandlers(
      flotsam:::ltsa_rspectra_ritz_eig(
        Matrix::Matrix(B, sparse = TRUE),
        ndim = 2L,
        dense_n = 0L,
        tol = 1e-10,
        maxitr = 5000L,
        gap_expansion_steps = 2L,
        verbose = TRUE
      ),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    "boundary gap after ndim is weak"
  )
  expect_true(any(grepl(
    paste0(
      "LTSA Ritz boundary gap remains below tolerance ",
      "after diagnostic expansion to 18 candidates; ",
      "returning residual-good, rank-good Ritz vectors"
    ),
    messages
  )))

  requested <- vapply(res$attempts, `[[`, integer(1), "eig_k")
  expect_equal(requested, c(8L, 12L, 18L))
  expect_equal(res$eig_k, 8L)
  expect_lt(max(requested), flotsam:::ltsa_max_ritz_candidate_k(2L, nrow(B)))
  expect_true(res$acceptance$resid_ok)
  expect_true(res$acceptance$rank_ok)
  expect_true(res$acceptance$ok)
  expect_false(res$acceptance$gap_ok)
  expect_true(res$acceptance$spectral_ambiguity_warning)
  expect_identical(
    res$acceptance$return_reason,
    "weak_lowest_energy_residual_rank_good"
  )
  expect_equal(res$acceptance$diagnostic_final_eig_k, 18L)
})

test_that("high selected Ritz residual keeps expanding to the hard cap", {
  B <- synthetic_ltsa_matrix(c(0, 0.1, 0.2, 1, 3, 5, 8, 13, 21, 34, 55, 89))

  res <- NULL
  expect_warning(
    res <- flotsam:::ltsa_rspectra_ritz_eig(
      Matrix::Matrix(B, sparse = TRUE),
      ndim = 2L,
      dense_n = 0L,
      tol = 1e-10,
      maxitr = 5000L,
      resid_tol = 1e-20,
      max_extra = 8L,
      gap_expansion_steps = 0L
    ),
    "residuals remain above tolerance"
  )

  requested <- vapply(res$attempts, `[[`, integer(1), "eig_k")
  expect_equal(
    max(requested),
    flotsam:::ltsa_max_ritz_candidate_k(
      2L,
      nrow(B),
      max_extra = 8L
    )
  )
  expect_false(res$acceptance$resid_ok)
  expect_false(tail(res$attempts, 1L)[[1L]]$accepted)
})

test_that("RSpectra partial convergence is a hard LTSA error", {
  set.seed(1)
  A <- crossprod(matrix(stats::rnorm(80L * 80L), nrow = 80L))
  B <- Matrix::Matrix(A, sparse = TRUE)
  lambda_max <- max(eigen(A, symmetric = TRUE, only.values = TRUE)$values)

  expect_error(
    suppressWarnings(
      flotsam:::rs_eig(
        B,
        k = 20L,
        lambda_max = lambda_max,
        dense_n = 0L,
        ncv = 21L,
        maxitr = 1L,
        tol = 1e-16
      )
    ),
    "RSpectra failed to converge enough LTSA candidate vectors: 0 / 20"
  )
})
