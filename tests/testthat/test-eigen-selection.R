make_selection_test_basis <- function() {
  list(
    u = c(1, 1, 1, 1) / 2,
    v1 = c(1, -1, 0, 0) / sqrt(2),
    v2 = c(0, 0, 1, -1) / sqrt(2),
    v3 = c(1, 1, -1, -1) / 2
  )
}

make_selection_test_matrix <- function(basis) {
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

make_centered_test_basis <- function(n, rank) {
  set.seed(123)
  null_vector <- rep(1, n) / sqrt(n)
  Z <- matrix(stats::rnorm(n * rank * 2L), nrow = n)
  Z <- Z - null_vector %*% crossprod(null_vector, Z)
  Q <- orthonormalize_test_basis(Z)
  Q[, seq_len(rank), drop = FALSE]
}

synthetic_eigenproblem <- function(values) {
  n <- length(values)
  basis <- cbind(rep(1, n) / sqrt(n), make_centered_test_basis(n, n - 1L))
  list(
    matrix = basis %*% diag(values) %*% t(basis),
    basis = basis
  )
}

synthetic_operator <- function(values) {
  synthetic_eigenproblem(values)$matrix
}

synthetic_candidate_provider_factory <- function(
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
    lambda_max <- max(
      eigen(
        problem$matrix,
        symmetric = TRUE,
        only.values = TRUE
      )$values
    )
  }
  lambda_max_value <- lambda_max
  state <- new.env(parent = emptyenv())
  state$calls <- data.frame(eig_k = integer())
  provider <- function(B, eig_k, lambda_max = NULL, verbose = FALSE) {
    state$calls <- rbind(state$calls, data.frame(eig_k = eig_k))
    cols_i <- cols(eig_k)
    nconv_i <- if (is.function(nconv)) nconv(eig_k) else nconv
    converged_i <- if (is.function(converged_columns)) {
      converged_columns(eig_k)
    } else {
      converged_columns
    }

    flotsam:::new_eigen_candidates(
      vectors = problem$basis[, cols_i, drop = FALSE],
      backend = backend,
      lambda_max = lambda_max_value,
      nconv = nconv_i,
      niter = niter,
      nops = nops,
      mprod = mprod,
      convergence_known = convergence_known,
      converged_columns = converged_i
    )
  }

  list(
    provider = provider,
    calls = function() state$calls
  )
}

test_that("embedding vector selection drops a returned trivial vector", {
  basis <- make_selection_test_basis()
  B <- make_selection_test_matrix(basis)
  candidates <- cbind(basis$v2, basis$u, basis$v3, basis$v1)

  selected <- flotsam:::select_ritz_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 3
  )$vectors

  expect_selected_basis(selected, cbind(basis$v1, basis$v2))
})

test_that("embedding vector selection keeps smallest vectors when trivial vector is absent", {
  basis <- make_selection_test_basis()
  B <- make_selection_test_matrix(basis)
  candidates <- cbind(basis$v2, basis$v3, basis$v1)

  selected <- flotsam:::select_ritz_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 3
  )$vectors

  expect_selected_basis(selected, cbind(basis$v1, basis$v2))
})

test_that("Ritz selection is invariant to rotations in the candidate subspace", {
  basis <- make_selection_test_basis()
  B <- make_selection_test_matrix(basis)
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

  ritz_result <- flotsam:::select_ritz_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 3
  )

  expect_same_subspace(ritz_result$vectors, target)
  expect_equal(ritz_result$values, c(1, 2), tolerance = 1e-12)
  expect_lt(max(ritz_result$residuals), 1e-12)
  expect_identical(ritz_result$rank_after_null, 3L)
  expect_equal(
    ritz_result$ritz_values[[3L]] - ritz_result$ritz_values[[2L]],
    1,
    tolerance = 1e-12
  )
})

test_that("clustered low-eigenvalue subspaces are compared by projector", {
  basis <- make_selection_test_basis()
  Q <- do.call(cbind, basis)
  B <- Q %*% diag(c(0, 1, 1, 4)) %*% t(Q)
  target <- cbind(basis$v1, basis$v2)
  candidates <- cbind(
    basis$u,
    (basis$v1 + basis$v2) / sqrt(2),
    (basis$v1 - basis$v2) / sqrt(2),
    basis$v3
  )

  ritz_result <- flotsam:::select_ritz_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 4
  )

  expect_same_subspace(ritz_result$vectors, target)
  expect_equal(ritz_result$values, c(1, 1), tolerance = 1e-12)
  expect_lt(max(ritz_result$residuals), 1e-12)
  expect_equal(
    ritz_result$ritz_values[[3L]] - ritz_result$ritz_values[[2L]],
    3,
    tolerance = 1e-12
  )
})

test_that("rank loss after null projection fails clearly", {
  basis <- make_selection_test_basis()
  B <- make_selection_test_matrix(basis)
  candidates <- cbind(basis$u, 2 * basis$u, basis$v1)

  expect_error(
    flotsam:::select_ritz_vectors(
      B,
      candidates,
      ndim = 2L,
      lambda_max = 3
    ),
    "rank after null projection is 1, less than ndim = 2"
  )
})

test_that("Ritz residual diagnostics use the scaled residual convention", {
  B <- synthetic_operator(c(0, 0.25, 0.5, 2, 8))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  candidates <- dense$vectors[, ord[seq_len(4L)], drop = FALSE]

  ritz_result <- flotsam:::select_ritz_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = 8
  )

  expect_equal(ritz_result$values, c(0.25, 0.5), tolerance = 1e-12)
  residuals <- flotsam:::ritz_residuals(
    B,
    ritz_result$vectors,
    ritz_result$values,
    8
  )
  expect_lt(max(residuals$absolute_residuals), 1e-12)
  expect_equal(
    ritz_result$residuals,
    residuals$absolute_residuals / 8,
    tolerance = 1e-15
  )
})

test_that("near-zero truncation uses the observed candidate span", {
  ndim <- 2L
  threshold <- flotsam:::near_zero_tolerance(4)
  cases <- list(
    below = list(
      values = c(0, threshold / 2, 0.1, 1, 2, 3, 4),
      count = 1L
    ),
    equal = list(
      values = c(0, threshold / 2, threshold / 2, 0.1, 1, 2, 4),
      count = 2L
    ),
    above = list(
      values = c(
        0,
        threshold / 2,
        threshold / 2,
        threshold / 2,
        0.1,
        1,
        4
      ),
      count = 3L
    ),
    fills_candidate_span = list(
      values = c(
        0,
        threshold / 2,
        threshold / 2,
        threshold / 2,
        threshold / 2,
        threshold / 2,
        4
      ),
      count = 5L
    )
  )

  for (case in cases) {
    problem <- synthetic_eigenproblem(case$values)
    fixture <- synthetic_candidate_provider_factory(problem)
    result <- flotsam:::solve_with_ritz(
      problem$matrix,
      ndim = ndim,
      provider = fixture$provider,
      eig_k = 6L
    )
    diagnostics <- result$eigen$diagnostics

    expect_identical(
      diagnostics$near_zero_block_truncated,
      case$count > ndim
    )
  }
})

test_that("small Ritz-selected cases agree with dense eigen reference subspaces", {
  B <- synthetic_operator(c(0, 0.1, 0.2, 0.8, 1.5, 3))
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

  ritz_result <- flotsam:::select_ritz_vectors(
    B,
    candidates,
    ndim = 2L,
    lambda_max = lambda_max
  )

  expect_equal(ritz_result$values, dense$values[ord[2:3]], tolerance = 1e-12)
  expect_same_subspace(ritz_result$vectors, reference, tolerance = 1e-7)
  expect_lt(max(ritz_result$residuals), 1e-12)
})

test_that("default eig_k follows the public rule", {
  expect_identical(flotsam:::default_eig_k(ndim = 2L, n = 50L), 12L)
  expect_identical(flotsam:::default_eig_k(ndim = 15L, n = 50L), 17L)
  expect_identical(flotsam:::default_eig_k(ndim = 2L, n = 8L), 7L)
  expect_identical(flotsam:::validate_eig_k(NULL, ndim = 2L, n = 50L), 12L)
})

test_that("Ritz driver accepts minimum eig_k equal to ndim plus one", {
  problem <- synthetic_eigenproblem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- synthetic_candidate_provider_factory(problem)

  res <- flotsam:::solve_with_ritz(
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

test_that("eig_k validation rejects values below ndim plus one", {
  expect_error(
    flotsam:::validate_eig_k(2L, ndim = 2L, n = 6L),
    "ndim \\+ 1 <= eig_k < n"
  )
})

test_that("eig_k validation rejects values at least n", {
  expect_error(
    flotsam:::validate_eig_k(6L, ndim = 2L, n = 6L),
    "ndim \\+ 1 <= eig_k < n"
  )
})

test_that("Ritz driver calls the provider exactly once", {
  problem <- synthetic_eigenproblem(c(
    0,
    0.1,
    0.2,
    1,
    2,
    3,
    5,
    8,
    13,
    21,
    34,
    55,
    89,
    144,
    233,
    377
  ))
  fixture <- synthetic_candidate_provider_factory(problem)

  res <- flotsam:::solve_with_ritz(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider
  )

  expect_equal(nrow(fixture$calls()), 1L)
  expect_equal(fixture$calls()$eig_k, 12L)
  expect_identical(res$eigen$eig_k, 12L)
  expect_identical(res$eigen$status, "ok")
})

test_that("public iterative LTSA honors explicit eig_k", {
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
      output = "result"
    )

    expect_identical(result$eigen$method, method)
    expect_identical(result$eigen$backend$name, method)
    expect_identical(result$eigen$eig_k, 4L)
    expect_false("attempts" %in% names(result))
    expect_false("acceptance" %in% names(result$eigen))
    expect_false("boundary_gap" %in% names(result$eigen))
  }
})

test_that("Ritz diagnostics use compact solver-neutral shape", {
  problem <- synthetic_eigenproblem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- synthetic_candidate_provider_factory(problem)

  res <- flotsam:::solve_with_ritz(
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
      "backend",
      "diagnostics"
    )
  )
  expect_equal(res$eigen$values, res$values, tolerance = 1e-12)
  expect_equal(res$eigen$residuals, rep(0, 2L), tolerance = 1e-12)
  expect_identical(res$eigen$method, "synthetic")
  expect_identical(res$eigen$backend$name, "synthetic")
  expect_named(
    res$eigen$diagnostics,
    "near_zero_block_truncated"
  )
  expect_false(any(c("attempts", "acceptance", "ritz") %in% names(res)))
})

test_that("Ritz diagnostics mark RSpectra nconv shortfall invalid", {
  problem <- synthetic_eigenproblem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- synthetic_candidate_provider_factory(
    problem,
    backend = "rspectra",
    convergence_known = TRUE,
    nconv = function(eig_k) eig_k - 1L,
    converged_columns = function(eig_k) eig_k - 1L,
    niter = 11L,
    nops = 101L
  )

  res <- flotsam:::solve_with_ritz(
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

test_that("Ritz diagnostics report generic shortfalls and weak gaps", {
  problem <- synthetic_eigenproblem(c(0, 0.1, 0.2, 0.20001, 1, 2))
  fixture <- synthetic_candidate_provider_factory(
    problem,
    backend = "synthetic_backend",
    convergence_known = TRUE,
    converged_columns = function(eig_k) eig_k - 1L
  )

  res <- flotsam:::solve_with_ritz(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 4L
  )

  expect_identical(res$eigen$status, "invalid")
  expect_true(any(grepl(
    "Backend reported fewer converged LTSA candidate vectors than requested",
    res$eigen$messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Weak Ritz boundary gap after the selected block",
    res$eigen$messages,
    fixed = TRUE
  )))
})

test_that("Ritz diagnostics do not invent non-RSpectra convergence", {
  problem <- synthetic_eigenproblem(c(0, 0.1, 0.2, 1, 2, 3))
  fixture <- synthetic_candidate_provider_factory(
    problem,
    backend = "irlba",
    convergence_known = FALSE,
    converged_columns = function(eig_k) NA_integer_,
    niter = 7L,
    mprod = 31L
  )

  res <- flotsam:::solve_with_ritz(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 4L
  )

  expect_identical(res$eigen$status, "ok")
  expect_identical(res$eigen$backend$name, "irlba")
  expect_false(res$eigen$backend$convergence_known)
  expect_identical(res$eigen$backend$iter, 7L)
  expect_identical(res$eigen$backend$mprod, 31L)
  expect_false("nconv" %in% names(res$eigen$backend))
  expect_length(res$eigen$messages, 0L)
})

test_that("Ritz diagnostics give near-zero truncation guidance", {
  problem <- synthetic_eigenproblem(c(0, 1e-10, 1e-10, 1e-10, 0.002, 1))
  fixture <- synthetic_candidate_provider_factory(problem)

  res <- flotsam:::solve_with_ritz(
    problem$matrix,
    ndim = 2L,
    provider = fixture$provider,
    eig_k = 4L
  )

  expect_equal(nrow(fixture$calls()), 1L)
  expect_identical(res$eigen$status, "warning")
  expect_true(any(grepl(
    "individual coordinates are not uniquely identifiable",
    res$eigen$messages,
    fixed = TRUE
  )))
  expect_true(any(grepl("Increase eig_k", res$eigen$messages, fixed = TRUE)))
  expect_true(any(grepl(
    "does not resolve a genuine repeated eigenspace",
    res$eigen$messages,
    fixed = TRUE
  )))
})

test_that("Ritz driver handles arbitrary normalized-style null vectors", {
  null_vector <- c(1, 2, 3, 2, 1, 4)
  null_vector <- null_vector / sqrt(sum(null_vector * null_vector))
  Z <- cbind(
    c(1, -1, 0, 0, 0, 0),
    c(0, 1, -1, 0, 0, 0),
    c(0, 0, 1, -1, 0, 0),
    c(0, 0, 0, 1, -1, 0),
    c(0, 0, 0, 0, 1, -1)
  )
  Z <- Z - null_vector %*% crossprod(null_vector, Z)
  Q <- qr.Q(qr(cbind(null_vector, Z)))
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
    flotsam:::new_eigen_candidates(
      vectors = candidates[, seq_len(eig_k), drop = FALSE],
      backend = "synthetic",
      lambda_max = 5,
      convergence_known = TRUE,
      converged_columns = eig_k
    )
  }

  res <- flotsam:::solve_with_ritz(
    B,
    ndim = 2L,
    provider = provider,
    null_vector = null_vector,
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

  B_sym <- flotsam:::symmetrize_operator(B)

  expect_s4_class(B_sym, "dgCMatrix")
  expect_false(methods::is(B_sym, "dsCMatrix"))
  expect_true(Matrix::isSymmetric(B_sym))
  expect_equal(as.matrix(B_sym), 0.5 * (as.matrix(B) + t(as.matrix(B))))
})

test_that("automatic dense LTSA route returns the lowest candidate subspace", {
  B <- flotsam:::symmetrize_operator(Matrix::Diagonal(x = seq(0, 11)))
  reference <- diag(12L)[, seq_len(6L), drop = FALSE]

  res <- flotsam:::auto_candidate_provider(B, eig_k = 6L)

  expect_identical(res$backend, "dense_eigen")
  expect_identical(res$nconv, 6L)
  expect_equal(ncol(res$vectors), 6L)
  expect_same_subspace(res$vectors, reference, tolerance = 1e-12)
})

test_that("RSpectra provider reports progress and returns low candidates", {
  B <- flotsam:::symmetrize_operator(Matrix::Diagonal(x = seq(0, 29)))
  reference <- diag(30L)[, seq_len(6L), drop = FALSE]

  messages <- capture.output(
    res <- flotsam:::rspectra_candidate_provider(
      B,
      eig_k = 6L,
      tol = 1e-10,
      maxitr = 5000L,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_identical(res$backend, "rspectra")
  expect_gte(res$nconv, 6L)
  expect_same_subspace(res$vectors, reference, tolerance = 1e-7)
  expect_true(any(grepl("Decomposing shifted matrix", messages, fixed = TRUE)))
  expect_true(any(grepl("max scaled residual", messages, fixed = TRUE)))
})

test_that("RSpectra candidate provider returns candidate bundle", {
  B <- flotsam:::symmetrize_operator(Matrix::Diagonal(x = seq(0, 29)))

  res <- flotsam:::rspectra_candidate_provider(
    B,
    eig_k = 6L,
    tol = 1e-10,
    maxitr = 5000L
  )

  expect_identical(res$backend, "rspectra")
  expect_true(res$convergence_known)
  expect_gte(res$nconv, 6L)
  expect_true(is.finite(res$lambda_max))
  expect_equal(ncol(res$vectors), 6L)
})

test_that("RSpectra Ritz driver returns diagnostics", {
  B <- synthetic_operator(c(
    0,
    0.1,
    0.2,
    1,
    3,
    5,
    8,
    13,
    21,
    34,
    55,
    89,
    144,
    233,
    377,
    610,
    987,
    1597,
    2584,
    4181
  ))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  reference <- dense$vectors[, ord[2:3], drop = FALSE]

  res <- flotsam:::solve_with_ritz(
    Matrix::Matrix(B, sparse = TRUE),
    ndim = 2L,
    provider = flotsam:::rspectra_candidate_provider,
    provider_args = list(
      tol = 1e-8,
      ncv = 18L,
      maxitr = 5000L
    ),
    eig_k = 6L
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

test_that("irlba and svdr candidate providers return candidate bundles", {
  B <- flotsam:::symmetrize_operator(Matrix::Diagonal(x = seq(0, 29)))
  reference <- diag(30L)[, seq_len(6L), drop = FALSE]
  providers <- list(
    irlba = list(
      provider = flotsam:::irlba_candidate_provider,
      args = list(tol = 1e-10, maxit = 1000L)
    ),
    svdr = list(
      provider = flotsam:::svdr_candidate_provider,
      args = list(tol = 1e-10, it = 1000L)
    )
  )

  for (backend in names(providers)) {
    set.seed(11)
    res <- do.call(
      providers[[backend]]$provider,
      c(list(B = B, eig_k = 6L), providers[[backend]]$args)
    )

    expect_identical(res$backend, backend)
    expect_false(res$convergence_known)
    expect_true(is.na(res$nconv))
    expect_true(is.finite(res$mprod))
    expect_true(is.finite(res$lambda_max))
    expect_equal(ncol(res$vectors), 6L)
    expect_same_subspace(res$vectors, reference, tolerance = 1e-6)
  }
})

test_that("irlba and svdr Ritz drivers agree with dense reference subspaces", {
  # fmt: skip
  B <- synthetic_operator(c(
    0, 0.1, 0.2, 1, 3, 5, 8, 13, 21, 34,
    55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181
  ))
  dense <- eigen(B, symmetric = TRUE)
  ord <- order(dense$values)
  reference <- dense$vectors[, ord[2:3], drop = FALSE]
  backends <- list(
    irlba = list(
      provider = flotsam:::irlba_candidate_provider,
      args = list(tol = 1e-10, maxit = 5000L)
    ),
    svdr = list(
      provider = flotsam:::svdr_candidate_provider,
      args = list(tol = 1e-10, it = 5000L, extra = 12L)
    )
  )

  for (backend in names(backends)) {
    set.seed(42)
    res <- flotsam:::solve_with_ritz(
      B = Matrix::Matrix(B, sparse = TRUE),
      ndim = 2L,
      provider = backends[[backend]]$provider,
      provider_args = backends[[backend]]$args,
      eig_k = 8L
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
  B <- flotsam:::symmetrize_operator(Matrix::Matrix(A, sparse = TRUE))
  lambda_max <- max(eigen(A, symmetric = TRUE, only.values = TRUE)$values)

  expect_error(
    suppressWarnings(
      flotsam:::rspectra_candidate_provider(
        B,
        eig_k = 20L,
        lambda_max = lambda_max,
        ncv = 21L,
        maxitr = 1L,
        tol = 1e-16
      )
    ),
    "RSpectra failed to converge enough LTSA candidate vectors: 0 / 20"
  )
})

test_that("verbose public iterative backends report candidate residual summaries", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result",
    verbose = TRUE
  )
  backend_arguments <- list(
    rspectra = list(tol = 1e-8, maxitr = 5000L, ncv = 8L),
    irlba = list(tol = 1e-8, maxit = 5000L, reorth = TRUE),
    svdr = list(tol = 1e-8, it = 500L, extra = 2L)
  )
  backend_labels <- c(rspectra = "RSpectra", irlba = "irlba", svdr = "svdr")

  for (method in names(backend_arguments)) {
    set.seed(42)
    messages <- capture.output(
      invisible(do.call(
        ltsa,
        c(
          base_args,
          list(eig_method = method),
          backend_arguments[[method]]
        )
      )),
      type = "message"
    )

    expect_true(any(grepl(
      paste0(
        backend_labels[[method]],
        ".* LTSA candidate vectors; max scaled residual = "
      ),
      messages
    )))
  }
})
