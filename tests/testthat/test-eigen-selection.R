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

synthetic_ltsa_matrix <- function(values) {
  n <- length(values)
  Q <- cbind(rep(1, n) / sqrt(n), centered_test_basis(n, n - 1L))
  Q %*% diag(values) %*% t(Q)
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
  B <- synthetic_ltsa_matrix(c(0, 5e-9, 5e-8, 5e-7, 5e-6, 1))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)

  rr <- flotsam:::ltsa_ritz_select(
    B,
    dense$vectors[, ord, drop = FALSE],
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
  B <- synthetic_ltsa_matrix(c(0, 5e-9, 5e-8, 5e-7, 5e-6, 1))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)

  rr <- flotsam:::ltsa_ritz_select(
    B,
    dense$vectors[, ord, drop = FALSE],
    ndim = 2L,
    lambda_max = 1
  )

  expect_equal(
    rr$near_zero_nonconstant_counts,
    c("1e-08" = 1L, "1e-07" = 2L, "1e-06" = 3L, "1e-05" = 4L)
  )
  expect_equal(
    rr$reported_ritz_values,
    c(5e-9, 5e-8, 5e-7, 5e-6, 1),
    tolerance = 1e-12
  )
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

  expect_true(all(c(
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
  ) %in% names(res)))
  expect_identical(res$backend, "rspectra")
  expect_identical(res$eig_k, 6L)
  expect_true(res$convergence_known)
  expect_gte(res$nconv, 6L)
  expect_equal(res$values, seq(0, 5), tolerance = 1e-8)
  expect_equal(ncol(res$vectors), 6L)
  expect_s4_class(res$matrix, "dgCMatrix")
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

test_that("weak-gap-only adaptive expansion stops before max_extra", {
  # fmt: skip
  B <- synthetic_ltsa_matrix(c(
    0, 0.1, 0.2, 0.200001, 1, 2, 3, 5, 8, 13,
    21, 34, 55, 89, 144, 233, 377, 610, 987, 1597
  ))

  messages <- character()
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
  expect_identical(
    res$acceptance$return_reason,
    "weak_first_residual_rank_good"
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
