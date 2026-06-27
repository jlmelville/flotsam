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

test_that("public iterative LTSA honors explicit fixed eig_k", {
  methods <- c("rspectra", "irlba", "svdr")

  for (method in methods) {
    set.seed(17)
    result <- ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = method,
      eig_k = 4L,
      output = "result",
      dense_n = 0L
    )

    expect_identical(result$eigen$method, method)
    expect_identical(result$eigen$eig_k, 4L)
    expect_false("attempts" %in% names(result))
    expect_false("acceptance" %in% names(result$eigen))
    expect_false("boundary_gap" %in% names(result$eigen))
  }
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

test_that("fixed-width diagnostics give eig_k and backend-setting guidance", {
  problem <- synthetic_ltsa_problem(c(0, 1e-10, 1e-3, 2e-3, 1))
  fixture <- fixed_width_provider_factory(problem)

  res <- flotsam:::ltsa_fixed_ritz_eig(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 4L
  )

  expect_equal(nrow(fixture$calls()), 1L)
  expect_identical(res$eigen$status, "warning")
  expect_true(any(grepl(
    "Only part of a near-zero selected Ritz block",
    res$eigen$messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "consider increasing eig_k",
    res$eigen$messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "stricter backend settings",
    res$eigen$messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "not completeness certificates",
    res$eigen$messages,
    fixed = TRUE
  )))
})

test_that("runtime eigensolver surfaces reject rescue-policy arguments", {
  B <- synthetic_ltsa_matrix(c(0, 0.1, 0.2, 1, 2, 3))

  expect_error(
    flotsam:::ltsa_rspectra_ritz_eig(
      Matrix::Matrix(B, sparse = TRUE),
      ndim = 2L,
      strict_rescue = FALSE
    ),
    "no longer supports rescue-policy"
  )
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_k = 4L,
      initial_extra = 4L
    ),
    "no longer supports rescue-policy"
  )
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

test_that("RSpectra Ritz wrapper returns fixed-width diagnostics", {
  B <- synthetic_ltsa_matrix(c(
    0, 0.1, 0.2, 1, 3, 5, 8, 13, 21, 34,
    55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181
  ))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  reference <- dense$vectors[, ord[2:3], drop = FALSE]

  res <- flotsam:::ltsa_rspectra_ritz_eig(
    Matrix::Matrix(B, sparse = TRUE),
    ndim = 2L,
    eig_k = 6L,
    dense_n = 0L,
    tol = 1e-8,
    ncv = 18L,
    maxitr = 5000L
  )

  expect_identical(res$eigen$method, "rspectra")
  expect_identical(res$eigen$eig_k, 6L)
  expect_identical(res$eigen$backend$name, "rspectra")
  expect_gte(res$eigen$backend$nconv, 6L)
  expect_false("attempts" %in% names(res))
  expect_equal(res$values, dense$values[ord[2:3]], tolerance = 1e-8)
  expect_same_subspace(res$vectors, reference, tolerance = 1e-7)
  expect_lt(max(res$eigen$residuals), 1e-8)
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

test_that("fixed-width irlba and svdr wrappers agree with dense reference subspaces", {
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
          eig_k = 8L
        ),
        backends[[backend]]$args
      )
    )

    expect_identical(res$eigen$backend$name, backend)
    expect_identical(res$eigen$eig_k, 8L)
    expect_false("attempts" %in% names(res))
    expect_equal(res$values, dense$values[ord[2:3]], tolerance = 1e-6)
    expect_same_subspace(res$vectors, reference, tolerance = 1e-5)
    expect_lt(max(res$eigen$residuals), 1e-6)
  }
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
