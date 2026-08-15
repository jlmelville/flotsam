# Eigensolver backends return candidate vectors for the final Rayleigh-Ritz
# selection step, plus compact backend metadata for diagnostics.
new_eigen_candidates <- function(
  vectors,
  backend,
  lambda_max = NULL,
  nconv = NA_integer_,
  niter = NA_integer_,
  nops = NA_integer_,
  mprod = NA_integer_,
  convergence_known = FALSE,
  converged_columns = ifelse(is.na(nconv), NA_integer_, nconv)
) {
  list(
    vectors = as.matrix(vectors),
    backend = backend,
    lambda_max = lambda_max,
    nconv = nconv,
    niter = niter,
    nops = nops,
    mprod = mprod,
    convergence_known = isTRUE(convergence_known),
    converged_columns = as.integer(converged_columns)
  )
}

# Automatic selection preserves the dense threshold policy while keeping an
# explicit RSpectra request literal.
auto_candidate_provider <- function(
  B,
  eig_k,
  lambda_max = NULL,
  verbose = FALSE,
  dense_n = 100L,
  dense_fraction = 0.5
) {
  n <- ncol(B)

  if (
    use_dense_eigensolver(
      n,
      eig_k,
      dense_n = dense_n,
      dense_fraction = dense_fraction
    )
  ) {
    report_progress("Using dense eigenvalue decomposition", verbose = verbose)
    return(dense_candidate_provider(B, eig_k))
  }

  report_progress("Calling rspectra", verbose = verbose)
  rspectra_candidate_provider(
    B = B,
    eig_k = eig_k,
    lambda_max = lambda_max,
    verbose = verbose
  )
}

# RSpectra uses a shifted largest-algebraic formulation and treats nconv as a
# hard convergence gate before shared LTSA postprocessing sees the candidates.
rspectra_candidate_provider <- function(
  B,
  eig_k,
  ...,
  lambda_max = NULL,
  verbose = FALSE,
  shift_eps = 1e-6
) {
  varargs <- list(...)
  n <- ncol(B)

  lambda_probe <- NULL
  if (is.null(lambda_max)) {
    report_progress("Finding largest eigenvalue", verbose = verbose)
    lambda_probe <- rspectra_lambda_max_probe(B, varargs)
    lambda_max <- lambda_probe$value
  } else {
    lambda_max <- validate_backend_lambda_max(
      lambda_max,
      B,
      backend = "RSpectra"
    )
  }
  shift <- lambda_max + shift_margin(lambda_max, shift_eps)
  shifted_operator <- shift_for_smallest(B, shift)

  report_progress("Decomposing shifted matrix", verbose = verbose)
  opts <- rspectra_options(eig_k = eig_k, n = n)
  opts <- merge_named_lists(opts, varargs)
  args <- list(
    A = shifted_operator,
    k = eig_k,
    which = "LA",
    opts = opts
  )
  res <- do.call(RSpectra::eigs_sym, args)
  nconv <- res$nconv %||% NA_integer_
  if (!is.na(nconv) && nconv < eig_k) {
    stop(
      "RSpectra failed to converge enough LTSA candidate vectors: ",
      nconv,
      " / ",
      eig_k,
      call. = FALSE
    )
  }
  scored_candidates <- sort_and_score_candidates(
    B = B,
    vectors = res$vectors,
    eig_k = eig_k,
    shifted_values = res$values,
    lambda_max = lambda_max,
    backend = "RSpectra"
  )
  report_progress(
    "RSpectra converged ",
    ifelse(is.na(nconv), eig_k, nconv),
    " / ",
    eig_k,
    " LTSA candidate vectors; max scaled residual = ",
    signif(max(scored_candidates$residuals$scaled_residuals), 4),
    verbose = verbose
  )

  new_eigen_candidates(
    vectors = scored_candidates$vectors,
    backend = "rspectra",
    lambda_max = lambda_max,
    nconv = nconv,
    niter = res$niter %||% NA_integer_,
    nops = res$nops %||% NA_integer_,
    convergence_known = TRUE,
    converged_columns = ifelse(is.na(nconv), eig_k, nconv)
  )
}

irlba_lambda_max_probe <- function(B) {
  args <- list(
    A = B,
    nv = 1L,
    nu = 0L
  )
  probe <- do.call(irlba::irlba, args)
  lambda_max <- validate_backend_lambda_max(
    probe$d,
    B,
    backend = "irlba"
  )

  list(
    value = lambda_max,
    niter = probe$iter %||% NA_integer_,
    mprod = probe$mprod %||% NA_integer_,
    opts = args
  )
}

svdr_lambda_max_probe <- function(B) {
  args <- list(
    x = B,
    k = 1L
  )
  probe <- do.call(irlba::svdr, args)
  lambda_max <- validate_backend_lambda_max(probe$d, B, backend = "svdr")

  list(
    value = lambda_max,
    mprod = probe$mprod %||% NA_integer_,
    opts = args
  )
}

sort_and_score_candidates <- function(
  B,
  vectors,
  eig_k,
  lambda_max,
  shifted_values = NULL,
  backend
) {
  if (is.null(vectors) || ncol(vectors) < eig_k) {
    stop(
      backend,
      " returned fewer LTSA candidate vectors than requested",
      call. = FALSE
    )
  }

  vectors <- as.matrix(vectors[, seq_len(eig_k), drop = FALSE])
  values <- rayleigh_values(B, vectors)
  ord <- order(values)
  values <- values[ord]
  vectors <- vectors[, ord, drop = FALSE]
  if (!is.null(shifted_values) && length(shifted_values) >= eig_k) {
    shifted_values <- shifted_values[seq_len(eig_k)][ord]
  }
  residuals <- ritz_residuals(B, vectors, values, lambda_max)

  list(
    vectors = vectors,
    values = values,
    shifted_values = shifted_values,
    residuals = residuals
  )
}

# Shifted solves are re-valued against the original LTSA matrix before their
# candidate vectors reach common selection. Final selection uses B.
new_shifted_candidates <- function(
  B,
  vectors,
  eig_k,
  backend,
  lambda_max,
  shifted_values = NULL,
  niter = NA_integer_,
  nops = NA_integer_,
  mprod = NA_integer_,
  verbose = FALSE
) {
  scored_candidates <- sort_and_score_candidates(
    B = B,
    vectors = vectors,
    eig_k = eig_k,
    shifted_values = shifted_values,
    lambda_max = lambda_max,
    backend = backend
  )
  report_progress(
    backend,
    " returned ",
    eig_k,
    " LTSA candidate vectors; max scaled residual = ",
    signif(max(scored_candidates$residuals$scaled_residuals), 4),
    verbose = verbose
  )

  new_eigen_candidates(
    vectors = scored_candidates$vectors,
    backend = backend,
    lambda_max = lambda_max,
    nconv = NA_integer_,
    niter = niter,
    nops = nops,
    mprod = mprod,
    convergence_known = FALSE,
    converged_columns = NA_integer_
  )
}

# irlba candidate provider. It solves the shifted problem by requesting right
# singular vectors, then all final ordering is recomputed by Rayleigh values of
# the original LTSA matrix.
irlba_candidate_provider <- function(
  B,
  eig_k,
  ...,
  lambda_max = NULL,
  verbose = FALSE,
  shift_eps = 1e-6
) {
  varargs <- list(...)

  lambda_probe <- NULL
  if (is.null(lambda_max)) {
    report_progress("Finding largest eigenvalue", verbose = verbose)
    lambda_probe <- irlba_lambda_max_probe(B)
    lambda_max <- lambda_probe$value
  } else {
    lambda_max <- validate_backend_lambda_max(
      lambda_max,
      B,
      backend = "irlba"
    )
  }
  shift <- lambda_max + shift_margin(lambda_max, shift_eps)
  shifted_operator <- shift_for_smallest(B, shift)

  report_progress("Decomposing shifted matrix", verbose = verbose)
  args <- merge_named_lists(
    list(
      A = shifted_operator,
      nv = eig_k,
      nu = 0L
    ),
    varargs
  )
  res <- do.call(irlba::irlba, args)

  new_shifted_candidates(
    B = B,
    vectors = res$v,
    eig_k = eig_k,
    backend = "irlba",
    lambda_max = lambda_max,
    shifted_values = res$d,
    niter = res$iter %||% NA_integer_,
    mprod = res$mprod %||% NA_integer_,
    verbose = verbose
  )
}

# svdr candidate provider. Like irlba, this produces a candidate subspace only;
# the common Rayleigh-Ritz selector decides which nonconstant directions to
# return.
svdr_candidate_provider <- function(
  B,
  eig_k,
  ...,
  lambda_max = NULL,
  verbose = FALSE,
  shift_eps = 1e-6
) {
  varargs <- list(...)

  lambda_probe <- NULL
  if (is.null(lambda_max)) {
    report_progress("Finding largest eigenvalue", verbose = verbose)
    lambda_probe <- svdr_lambda_max_probe(B)
    lambda_max <- lambda_probe$value
  } else {
    lambda_max <- validate_backend_lambda_max(
      lambda_max,
      B,
      backend = "svdr"
    )
  }
  shift <- lambda_max + shift_margin(lambda_max, shift_eps)
  shifted_operator <- shift_for_smallest(B, shift)

  report_progress("Decomposing shifted matrix", verbose = verbose)
  args <- merge_named_lists(
    list(
      x = shifted_operator,
      k = eig_k
    ),
    varargs
  )
  res <- do.call(irlba::svdr, args)

  new_shifted_candidates(
    B = B,
    vectors = res$v,
    eig_k = eig_k,
    backend = "svdr",
    lambda_max = lambda_max,
    shifted_values = res$d,
    mprod = res$mprod %||% NA_integer_,
    verbose = verbose
  )
}

# Dense path used by automatic selection and explicit dense requests. The
# unused controls keep the provider interface uniform across backends.
dense_candidate_provider <- function(
  B,
  eig_k,
  lambda_max = NULL,
  verbose = FALSE,
  backend = "dense_eigen"
) {
  dense <- as.matrix(B)
  eig <- eigen(dense, symmetric = TRUE)
  ord <- order(eig$values)
  values_all <- eig$values[ord]
  vectors_all <- eig$vectors[, ord, drop = FALSE]
  lambda_max <- max(values_all)

  take <- seq_len(eig_k)
  vectors <- vectors_all[, take, drop = FALSE]

  new_eigen_candidates(
    vectors = vectors,
    backend = backend,
    lambda_max = lambda_max,
    nconv = eig_k,
    convergence_known = TRUE,
    converged_columns = eig_k
  )
}
