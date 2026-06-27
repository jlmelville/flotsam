# Null-aware Rayleigh-Ritz selection and fixed-width diagnostics for LTSA
# eigenanalysis.

# Wrapper-level Rayleigh-Ritz extraction. Iterative solvers commonly do their
# own Ritz extraction internally, but this step answers the LTSA-specific
# question: within the returned candidate span, which vectors best diagonalize
# the original LTSA matrix after the known constant null direction is removed?
ltsa_ritz_select <- function(
  B,
  vectors,
  ndim,
  nullvec = ltsa_default_null_vector(nrow(B)),
  lambda_max = NULL,
  drop_tol = 1e-8,
  rank_tol = 1e-10,
  near_zero_tol = ltsa_near_zero_tol(lambda_max),
  zero_tol = ltsa_gap_zero_tol(lambda_max)
) {
  vectors <- as.matrix(vectors)
  n <- nrow(B)
  ndim <- as.integer(ndim)
  if (length(ndim) != 1L || is.na(ndim) || ndim < 1L) {
    stop("ndim must be a positive integer", call. = FALSE)
  }
  if (nrow(B) != ncol(B)) {
    stop("LTSA matrix must be square", call. = FALSE)
  }
  if (nrow(vectors) != n) {
    stop(
      "LTSA candidate vectors must match the matrix dimension",
      call. = FALSE
    )
  }
  if (ncol(vectors) < ndim) {
    stop("Can't find enough LTSA candidate vectors", call. = FALSE)
  }

  nullvec <- ltsa_normalize_null_vector(nullvec, n)
  projected <- vectors - nullvec %*% crossprod(nullvec, vectors)
  original_norms <- sqrt(colSums(vectors * vectors))
  projected_norms <- sqrt(colSums(projected * projected))
  keep <- projected_norms > drop_tol * pmax(1, original_norms)
  projected <- projected[, keep, drop = FALSE]

  dropped_null_columns <- sum(!keep)
  if (ncol(projected) == 0L) {
    rank_after_null <- 0L
  } else {
    qrp <- qr(projected, tol = rank_tol)
    rank_after_null <- qrp$rank
  }
  if (rank_after_null < ndim) {
    stop(
      "LTSA candidate subspace rank after null projection is ",
      rank_after_null,
      ", less than ndim = ",
      ndim,
      call. = FALSE
    )
  }

  Q <- qr.Q(qrp)[, seq_len(rank_after_null), drop = FALSE]
  BQ <- as.matrix(B %*% Q)
  H <- as.matrix(crossprod(Q, BQ))
  H <- 0.5 * (H + t(H))

  eig <- eigen(H, symmetric = TRUE)
  ord <- order(eig$values)
  values_all <- eig$values[ord]
  coef_all <- eig$vectors[, ord, drop = FALSE]
  vectors_all <- Q %*% coef_all
  residuals_all <- ltsa_ritz_residuals(B, vectors_all, values_all, lambda_max)

  if (length(values_all) > ndim) {
    boundary_gap <- values_all[[ndim + 1L]] - values_all[[ndim]]
    global_gap <- boundary_gap / ltsa_residual_scale(lambda_max)
    local_gap <- boundary_gap /
      max(
        abs(values_all[[ndim + 1L]]),
        abs(values_all[[ndim]]),
        zero_tol
      )
  } else {
    boundary_gap <- NA_real_
    global_gap <- NA_real_
    local_gap <- NA_real_
  }

  take <- seq_len(ndim)
  near_zero_nonconstant_count <- sum(abs(values_all) <= near_zero_tol)
  near_zero_thresholds <- ltsa_near_zero_thresholds(lambda_max)
  near_zero_nonconstant_counts <- ltsa_near_zero_counts(
    values_all,
    near_zero_thresholds
  )
  reported_ritz_values <- values_all[seq_len(min(20L, length(values_all)))]

  list(
    vectors = vectors_all[, take, drop = FALSE],
    values = values_all[take],
    all_vectors = vectors_all,
    all_values = values_all,
    absolute_residuals = residuals_all$absolute_residuals[take],
    scaled_residuals = residuals_all$scaled_residuals[take],
    all_absolute_residuals = residuals_all$absolute_residuals,
    all_scaled_residuals = residuals_all$scaled_residuals,
    residual_scale = residuals_all$residual_scale,
    rank_after_null = rank_after_null,
    dropped_null_columns = dropped_null_columns,
    boundary_gap = boundary_gap,
    global_gap = global_gap,
    local_gap = local_gap,
    zero_tol = zero_tol,
    boundary_gap_relative = global_gap,
    near_zero_nonconstant_count = near_zero_nonconstant_count,
    near_zero_tol = near_zero_tol,
    near_zero_thresholds = near_zero_thresholds,
    near_zero_nonconstant_counts = near_zero_nonconstant_counts,
    reported_ritz_values = reported_ritz_values,
    kept_candidate_columns = sum(keep)
  )
}

select_ltsa_embedding_vectors <- function(
  B,
  vectors,
  ndim,
  nullvec = ltsa_default_null_vector(nrow(B)),
  lambda_max = NULL,
  ...
) {
  rr <- ltsa_ritz_select(
    B = B,
    vectors = vectors,
    ndim = ndim,
    nullvec = nullvec,
    lambda_max = lambda_max,
    ...
  )
  rr$vectors
}

ltsa_fixed_exact_dense_backend <- function(backend) {
  backend %in% c("dense_eigen", "eig", "eigen", "fullsvd", "dense_svd")
}

ltsa_fixed_backend_metadata <- function(eig_res) {
  backend <- eig_res$backend %||% "unknown"
  backend <- as.character(backend[[1L]])
  exact_dense <- ltsa_fixed_exact_dense_backend(backend)
  out <- list(
    name = backend,
    convergence_known = isTRUE(eig_res$convergence_known) || exact_dense
  )

  if (identical(backend, "rspectra")) {
    out$nconv <- as.integer(eig_res$nconv %||% NA_integer_)
    out$niter <- as.integer(eig_res$niter %||% NA_integer_)
    out$nops <- as.integer(eig_res$nops %||% NA_integer_)
  } else {
    niter <- as.integer(eig_res$niter %||% NA_integer_)
    mprod <- as.integer(eig_res$mprod %||% NA_integer_)
    if (!is.na(niter)) {
      out$iter <- niter
    }
    if (!is.na(mprod)) {
      out$mprod <- mprod
    }
    if (exact_dense) {
      out$exact_dense <- TRUE
    }
  }

  out
}

ltsa_fixed_backend_failure_messages <- function(eig_res, eig_k) {
  backend <- eig_res$backend %||% "unknown"
  backend <- as.character(backend[[1L]])
  nconv <- as.integer(eig_res$nconv %||% NA_integer_)
  converged_columns <- as.integer(eig_res$converged_columns %||% nconv)

  if (identical(backend, "rspectra") && !is.na(nconv) && nconv < eig_k) {
    return(paste0(
      "RSpectra converged fewer LTSA candidate vectors than requested: ",
      nconv,
      " / ",
      eig_k,
      "."
    ))
  }
  if (
    isTRUE(eig_res$convergence_known) &&
      !is.na(converged_columns) &&
      converged_columns < eig_k
  ) {
    return(paste0(
      "Backend reported fewer converged LTSA candidate vectors than ",
      "requested: ",
      converged_columns,
      " / ",
      eig_k,
      "."
    ))
  }

  character()
}

ltsa_fixed_diagnostic_guidance <- function() {
  paste0(
    "These diagnostics are not completeness certificates; consider ",
    "increasing eig_k or using stricter backend settings if the result ",
    "looks suspicious."
  )
}

ltsa_fixed_ritz_diagnostics <- function(
  eig_res,
  rr,
  eig_k,
  ndim,
  resid_tol,
  gap_tol
) {
  backend <- ltsa_fixed_backend_metadata(eig_res)
  lambda_max <- eig_res$lambda_max %||% NA_real_
  values <- as.numeric(rr$values)
  residuals <- as.numeric(rr$scaled_residuals)
  max_residual <- if (length(residuals) == 0L) {
    Inf
  } else {
    max(residuals)
  }

  invalid_messages <- character()
  if (length(values) < ndim || ncol(rr$vectors) < ndim) {
    invalid_messages <- c(
      invalid_messages,
      "Fewer than ndim usable nonconstant Ritz vectors were selected."
    )
  }
  if (rr$rank_after_null < ndim) {
    invalid_messages <- c(
      invalid_messages,
      paste0(
        "Post-null candidate rank ",
        rr$rank_after_null,
        " is less than ndim = ",
        ndim,
        "."
      )
    )
  }
  if (!all(is.finite(values)) || !all(is.finite(rr$vectors))) {
    invalid_messages <- c(
      invalid_messages,
      "Selected Ritz values or vectors contain non-finite entries."
    )
  }
  if (!is.finite(max_residual) || max_residual > resid_tol) {
    invalid_messages <- c(
      invalid_messages,
      paste0(
        "Selected scaled residuals exceed tolerance: max = ",
        signif(max_residual, 4),
        ", tolerance = ",
        signif(resid_tol, 4),
        "."
      )
    )
  }
  invalid_messages <- c(
    invalid_messages,
    ltsa_fixed_backend_failure_messages(eig_res, eig_k)
  )

  warning_messages <- character()
  if (length(rr$all_values) < ndim + 1L) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "Candidate span contains fewer than ndim + 1 post-null Ritz ",
        "values; no spare boundary direction is available."
      )
    )
  } else if (!is.finite(rr$global_gap)) {
    warning_messages <- c(
      warning_messages,
      "Ritz boundary gap after the selected block is unavailable."
    )
  } else if (rr$global_gap < gap_tol) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "Weak Ritz boundary gap after the selected block: ",
        signif(rr$global_gap, 4),
        " < ",
        signif(gap_tol, 4),
        "."
      )
    )
  }
  if (!isTRUE(backend$convergence_known)) {
    warning_messages <- c(
      warning_messages,
      "Backend does not provide a native convergence certificate."
    )
  }
  resid_ok <- is.finite(max_residual) && max_residual <= resid_tol
  if (
    ltsa_partial_near_zero_block(
      list(
        near_zero_nonconstant_count = rr$near_zero_nonconstant_count,
        acceptance = list(
          rank_ok = rr$rank_after_null >= ndim,
          resid_ok = resid_ok
        )
      ),
      ndim = ndim
    )
  ) {
    warning_messages <- c(
      warning_messages,
      "Only part of a near-zero selected Ritz block is present."
    )
  }

  status <- if (length(invalid_messages) > 0L) {
    "invalid"
  } else if (length(warning_messages) > 0L) {
    "warning"
  } else {
    "ok"
  }
  messages <- c(invalid_messages, warning_messages)
  if (length(messages) > 0L) {
    messages <- c(messages, ltsa_fixed_diagnostic_guidance())
  }

  list(
    method = backend$name,
    eig_k = eig_k,
    values = values,
    ritz_values = as.numeric(rr$all_values),
    residuals = residuals,
    rank = rr$rank_after_null,
    lambda_max = lambda_max,
    status = status,
    messages = messages,
    backend = backend
  )
}

ltsa_fixed_ritz_eig <- function(
  B,
  ndim,
  provider,
  provider_args = list(),
  nullvec = ltsa_default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE
) {
  eig_k <- ltsa_validate_eig_k(eig_k = eig_k, ndim = ndim, n = nrow(B))
  eig_res <- ltsa_call_candidate_provider(
    provider = provider,
    B = B,
    eig_k = eig_k,
    provider_args = provider_args,
    lambda_max = NULL,
    verbose = verbose
  )
  eig_res <- ltsa_as_candidate_result(
    eig_res,
    B = B,
    eig_k = eig_k,
    backend = eig_res$backend %||% "unknown",
    lambda_max = eig_res$lambda_max %||% NULL,
    convergence_known = eig_res$convergence_known %||% FALSE
  )
  lambda_max <- eig_res$lambda_max %||% NA_real_
  rr <- ltsa_ritz_select(
    B = eig_res$matrix,
    vectors = eig_res$vectors,
    ndim = ndim,
    nullvec = nullvec,
    lambda_max = lambda_max
  )
  eigen <- ltsa_fixed_ritz_diagnostics(
    eig_res = eig_res,
    rr = rr,
    eig_k = eig_k,
    ndim = ndim,
    resid_tol = resid_tol,
    gap_tol = gap_tol
  )

  list(
    vectors = rr$vectors,
    values = rr$values,
    eigen = eigen,
    backend = eigen$backend,
    lambda_max = eigen$lambda_max,
    eig_k = eig_k
  )
}

# Internal iterative wrappers retained as fixed-width shims. Adaptive widening
# and rescue policy are rejected at the public argument boundary.
ltsa_rspectra_ritz_eig <- function(
  B,
  ndim,
  ...,
  nullvec = ltsa_default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE
) {
  provider_args <- list(...)
  ltsa_reject_rescue_policy_args(provider_args)

  ltsa_fixed_ritz_eig(
    B = B,
    ndim = ndim,
    provider = ltsa_rspectra_candidate_provider,
    provider_args = provider_args,
    nullvec = nullvec,
    eig_k = eig_k,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    verbose = verbose
  )
}

ltsa_irlba_ritz_eig <- function(
  B,
  ndim,
  ...,
  nullvec = ltsa_default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE
) {
  provider_args <- list(...)
  ltsa_reject_rescue_policy_args(provider_args)

  ltsa_fixed_ritz_eig(
    B = B,
    ndim = ndim,
    provider = ltsa_irlba_candidate_provider,
    provider_args = provider_args,
    nullvec = nullvec,
    eig_k = eig_k,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    verbose = verbose
  )
}

ltsa_svdr_ritz_eig <- function(
  B,
  ndim,
  ...,
  nullvec = ltsa_default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE
) {
  provider_args <- list(...)
  ltsa_reject_rescue_policy_args(provider_args)

  ltsa_fixed_ritz_eig(
    B = B,
    ndim = ndim,
    provider = ltsa_svdr_candidate_provider,
    provider_args = provider_args,
    nullvec = nullvec,
    eig_k = eig_k,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    verbose = verbose
  )
}

# Fixed-width diagnostics still flag the specific case where the selected block
# appears to cut through a near-zero nonconstant cluster.
ltsa_partial_near_zero_block <- function(res, ndim) {
  if (
    is.null(res) ||
      !isTRUE(res$acceptance$rank_ok) ||
      !isTRUE(res$acceptance$resid_ok)
  ) {
    return(FALSE)
  }

  near_zero_count <- res$near_zero_nonconstant_count %||% 0L
  near_zero_count > 0L && near_zero_count < ndim
}
