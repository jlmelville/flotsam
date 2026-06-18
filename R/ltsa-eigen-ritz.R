# Null-aware Rayleigh-Ritz selection, acceptance diagnostics, ambiguity
# warnings, and strict-rescue helpers for LTSA eigenanalysis.

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

ltsa_initial_ritz_candidate_k <- function(ndim, n, initial_extra = 4L) {
  min(n - 1L, 1L + ndim + 1L + initial_extra)
}

ltsa_max_ritz_candidate_k <- function(ndim, n, max_extra = 40L) {
  min(n - 1L, 1L + ndim + 1L + max_extra)
}

ltsa_next_ritz_candidate_k <- function(eig_k, max_k) {
  min(max_k, max(eig_k + 2L, as.integer(ceiling(1.5 * eig_k))))
}

# Acceptance is intentionally split into hard checks and diagnostic checks.
# Post-null rank and scaled residuals are hard quality gates. The boundary gap
# is advisory by default: a weak gap can mean the embedding is part of a larger
# low-energy eigenspace rather than numerically wrong.
accept_ltsa_ritz <- function(
  rr,
  ndim,
  lambda_max,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  require_gap = FALSE
) {
  max_scaled_residual <- max(rr$scaled_residuals)
  gap <- rr$global_gap %||% rr$boundary_gap_relative
  rank_ok <- rr$rank_after_null >= ndim
  resid_ok <- is.finite(max_scaled_residual) && max_scaled_residual <= resid_tol
  gap_available <- is.finite(gap)
  gap_ok <- gap_available && gap >= gap_tol
  gap_status <- if (!gap_available) {
    "unavailable"
  } else if (gap_ok) {
    "ok"
  } else {
    "weak"
  }

  list(
    ok = rank_ok && resid_ok && (!require_gap || gap_ok),
    rank_ok = rank_ok,
    resid_ok = resid_ok,
    gap_ok = gap_ok,
    gap_available = gap_available,
    gap_status = gap_status,
    gap_issue_only = rank_ok && resid_ok && !gap_ok,
    global_gap = gap,
    local_gap = rr$local_gap %||% NA_real_,
    zero_tol = rr$zero_tol %||% NA_real_,
    rel_gap = gap,
    max_scaled_residual = max_scaled_residual,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    lambda_max = lambda_max,
    require_gap = require_gap
  )
}

ltsa_result_energy_key <- function(res, ndim) {
  values <- as.numeric(res$values)
  if (length(values) < ndim) {
    return(list(lambda_ndim = Inf, trace = Inf, near_zero_count = 0L))
  }

  list(
    lambda_ndim = values[[ndim]],
    trace = sum(values[seq_len(ndim)]),
    near_zero_count = res$near_zero_nonconstant_count %||% 0L
  )
}

ltsa_energy_better <- function(candidate, incumbent, ndim) {
  candidate_key <- ltsa_result_energy_key(candidate, ndim)
  incumbent_key <- ltsa_result_energy_key(incumbent, ndim)

  scale <- max(
    abs(candidate_key$lambda_ndim),
    abs(incumbent_key$lambda_ndim),
    1
  )
  tol <- max(100 * .Machine$double.eps * scale, 1e-8 * scale)
  if (candidate_key$lambda_ndim < incumbent_key$lambda_ndim - tol) {
    return(TRUE)
  }
  if (candidate_key$lambda_ndim > incumbent_key$lambda_ndim + tol) {
    return(FALSE)
  }

  trace_scale <- max(abs(candidate_key$trace), abs(incumbent_key$trace), 1)
  trace_tol <- max(100 * .Machine$double.eps * trace_scale, 1e-8 * trace_scale)
  if (candidate_key$trace < incumbent_key$trace - trace_tol) {
    return(TRUE)
  }
  if (candidate_key$trace > incumbent_key$trace + trace_tol) {
    return(FALSE)
  }

  candidate_key$near_zero_count > incumbent_key$near_zero_count
}

ltsa_rescue_candidate <- function(candidate, incumbent, ndim) {
  if (is.null(candidate)) {
    return(incumbent)
  }
  if (
    !isTRUE(candidate$acceptance$rank_ok) ||
      !isTRUE(candidate$acceptance$resid_ok)
  ) {
    return(incumbent)
  }
  if (is.null(incumbent)) {
    return(candidate)
  }
  if (
    !isTRUE(incumbent$acceptance$rank_ok) ||
      !isTRUE(incumbent$acceptance$resid_ok)
  ) {
    return(candidate)
  }
  if (ltsa_energy_better(candidate, incumbent, ndim)) {
    return(candidate)
  }
  incumbent
}

# Strict rescue targets one specific failure signature: residual-good,
# rank-good output where the selected block appears to contain only part of a
# near-zero nonconstant cluster. Rayleigh-Ritz cannot invent a missing
# direction, so the rescue path reruns the backend with stricter settings and
# also tries a combined candidate bank.
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

ltsa_strict_rescue_needed <- function(res, ndim) {
  ltsa_partial_near_zero_block(res, ndim)
}

ltsa_mark_partial_near_zero_block <- function(res, ndim) {
  partial <- ltsa_partial_near_zero_block(res, ndim)
  res$partial_near_zero_block <- partial
  res$acceptance$partial_near_zero_block <- partial
  res
}

ltsa_maybe_warn_partial_near_zero_block <- function(
  res,
  ndim,
  strict_rescue
) {
  res <- ltsa_mark_partial_near_zero_block(res, ndim)
  if (!isTRUE(strict_rescue) && isTRUE(res$partial_near_zero_block)) {
    warning(
      "LTSA selected a residual-good/rank-good embedding but found only ",
      res$near_zero_nonconstant_count,
      " near-zero nonconstant mode",
      if (isTRUE(res$near_zero_nonconstant_count == 1L)) "" else "s",
      " for ndim = ",
      ndim,
      ". This pattern can indicate a missing near-zero LTSA coordinate; ",
      "consider strict_rescue = TRUE or a larger candidate budget.",
      call. = FALSE
    )
  }
  res
}

# Ambiguity warnings are only issued for residual-good, rank-good results. They
# communicate spectral non-uniqueness without turning a plausible embedding into
# an error.
ltsa_spectral_ambiguity_issues <- function(res, ndim) {
  if (
    is.null(res) ||
      !isTRUE(res$acceptance$rank_ok) ||
      !isTRUE(res$acceptance$resid_ok)
  ) {
    return(character())
  }

  issues <- character()
  near_zero_count <- res$near_zero_nonconstant_count %||% 0L
  if (near_zero_count > ndim) {
    issues <- c(
      issues,
      paste0(
        near_zero_count,
        " near-zero nonconstant modes were found for ndim = ",
        ndim
      )
    )
  }

  if (
    isTRUE(res$acceptance$gap_available) &&
      !isTRUE(res$acceptance$gap_ok)
  ) {
    local_gap <- res$acceptance$local_gap %||% NA_real_
    gap_tol <- res$acceptance$gap_tol %||% NA_real_
    local_gap_ok <- is.finite(local_gap) &&
      is.finite(gap_tol) &&
      local_gap >= gap_tol
    if (!local_gap_ok) {
      issues <- c(
        issues,
        paste0(
          "the boundary gap after ndim is weak: global gap = ",
          signif(res$acceptance$global_gap %||% NA_real_, 4),
          ", local gap = ",
          signif(local_gap, 4),
          ", tolerance = ",
          signif(gap_tol, 4)
        )
      )
    }
  }

  issues
}

ltsa_maybe_warn_spectral_ambiguity <- function(res, ndim) {
  issues <- ltsa_spectral_ambiguity_issues(res, ndim)
  if (length(issues) == 0L) {
    return(res)
  }

  warning(
    "LTSA eigenanalysis found an ambiguous low-energy eigenspace; ",
    "embedding may not be unique up to only rotation/sign. Diagnostics: ",
    paste(issues, collapse = "; "),
    ". This can happen when the neighborhood graph is disconnected or ",
    "weakly connected, n_neighbors is too small, or ndim cuts through a ",
    "low-energy eigenspace.",
    call. = FALSE
  )
  res$acceptance$spectral_ambiguity_warning <- TRUE
  res$acceptance$spectral_ambiguity_issues <- issues
  res
}

# Small adapter kept for callers that want an explicit RSpectra provider call
# rather than going through the shared adaptive driver.
ltsa_rs_eig_call <- function(
  B,
  eig_k,
  varargs,
  lambda_max = NULL,
  verbose = FALSE
) {
  ltsa_call_candidate_provider(
    provider = ltsa_rspectra_candidate_provider,
    B = B,
    eig_k = eig_k,
    provider_args = varargs,
    lambda_max = lambda_max,
    verbose = verbose
  )
}

ltsa_strict_rescue_args <- function(
  varargs,
  strict_rescue_tol,
  strict_rescue_maxitr
) {
  strict_args <- varargs
  if (is.null(strict_args$tol)) {
    strict_args$tol <- strict_rescue_tol
  } else {
    strict_args$tol <- min(strict_args$tol, strict_rescue_tol)
  }
  if (is.null(strict_args$maxitr)) {
    strict_args$maxitr <- strict_rescue_maxitr
  } else {
    strict_args$maxitr <- max(strict_args$maxitr, strict_rescue_maxitr)
  }
  strict_args
}

ltsa_attempt_eig_ks <- function(attempts) {
  if (is.null(attempts) || length(attempts) == 0L) {
    return(integer())
  }

  eig_ks <- vapply(
    attempts,
    function(attempt) {
      eig_k <- attempt$eig_k %||% NA_integer_
      if (length(eig_k) != 1L || is.na(eig_k)) {
        return(NA_integer_)
      }
      as.integer(eig_k)
    },
    integer(1)
  )
  eig_ks[!is.na(eig_ks) & eig_ks > 0L]
}

ltsa_strict_rescue_eig_k <- function(selected, ndim, n, strict_rescue_extra) {
  candidates <- c(
    selected$eig_k %||% NA_integer_,
    selected$acceptance$diagnostic_final_eig_k %||% NA_integer_,
    ltsa_attempt_eig_ks(selected$attempts),
    ndim + strict_rescue_extra
  )
  candidates <- as.integer(candidates)
  candidates <- candidates[!is.na(candidates) & candidates > 0L]

  min(n - 1L, max(candidates))
}

ltsa_attempt_opt_integer <- function(opts, name) {
  if (is.null(opts) || !is.list(opts) || is.null(opts[[name]])) {
    return(NA_integer_)
  }

  value <- opts[[name]]
  if (length(value) == 0L || is.na(value[[1L]])) {
    return(NA_integer_)
  }

  out <- suppressWarnings(as.integer(value[[1L]]))
  if (length(out) != 1L || is.na(out)) {
    return(NA_integer_)
  }
  out
}

ltsa_candidate_reference_projection <- function(
  candidate_vectors,
  reference_vectors,
  nullvec = NULL,
  rank_tol = 1e-10
) {
  candidate_vectors <- as.matrix(candidate_vectors)
  reference_vectors <- as.matrix(reference_vectors)
  if (nrow(candidate_vectors) != nrow(reference_vectors)) {
    stop("candidate and reference vector row counts differ", call. = FALSE)
  }
  if (!is.null(nullvec)) {
    nullvec <- ltsa_normalize_null_vector(nullvec, nrow(candidate_vectors))
    candidate_vectors <- candidate_vectors -
      nullvec %*% crossprod(nullvec, candidate_vectors)
    reference_vectors <- reference_vectors -
      nullvec %*% crossprod(nullvec, reference_vectors)
  }
  if (ncol(reference_vectors) == 0L) {
    return(list(
      reference_projection_norms = numeric(),
      reference_projection_min_norm = NA_real_,
      reference_candidate_space_rank = NA_integer_,
      reference_space_rank = NA_integer_,
      reference_overlap_singular_values = numeric(),
      reference_overlap_min_singular_value = NA_real_
    ))
  }

  candidate_qrp <- qr(candidate_vectors, tol = rank_tol)
  candidate_rank <- candidate_qrp$rank
  reference_qrp <- qr(reference_vectors, tol = rank_tol)
  reference_rank <- reference_qrp$rank
  if (candidate_rank == 0L) {
    norms <- rep(0, ncol(reference_vectors))
  } else {
    candidate_Q <- qr.Q(candidate_qrp)[, seq_len(candidate_rank), drop = FALSE]
    projected_norms <- sqrt(colSums(crossprod(candidate_Q, reference_vectors)^2))
    reference_norms <- sqrt(colSums(reference_vectors * reference_vectors))
    norms <- projected_norms / pmax(reference_norms, .Machine$double.eps)
    norms <- pmax(0, pmin(1, norms))
  }
  singular_values <- if (reference_rank == 0L) {
    numeric()
  } else if (candidate_rank == 0L) {
    rep(0, reference_rank)
  } else {
    candidate_Q <- qr.Q(candidate_qrp)[, seq_len(candidate_rank), drop = FALSE]
    reference_Q <- qr.Q(reference_qrp)[, seq_len(reference_rank), drop = FALSE]
    sv <- svd(crossprod(candidate_Q, reference_Q), nu = 0, nv = 0)$d
    sv <- pmax(0, pmin(1, sv))
    if (length(sv) < reference_rank) {
      sv <- c(sv, rep(0, reference_rank - length(sv)))
    }
    sv
  }

  list(
    reference_projection_norms = norms,
    reference_projection_min_norm = if (length(norms) == 0L) {
      NA_real_
    } else {
      min(norms)
    },
    reference_candidate_space_rank = as.integer(candidate_rank),
    reference_space_rank = as.integer(reference_rank),
    reference_overlap_singular_values = singular_values,
    reference_overlap_min_singular_value = if (length(singular_values) == 0L) {
      NA_real_
    } else {
      min(singular_values)
    }
  )
}

ltsa_finalize_attempt_observability <- function(
  res,
  reference_vectors = NULL,
  nullvec = NULL,
  retain_candidate_spaces = FALSE
) {
  if (is.null(res) || is.null(res$attempts) || length(res$attempts) == 0L) {
    return(res)
  }
  reference_vectors <- reference_vectors %||% res$vectors
  if (is.null(reference_vectors)) {
    return(res)
  }

  res$attempts <- lapply(res$attempts, function(attempt) {
    candidate_vectors <- attempt$.candidate_vectors
    if (!is.null(candidate_vectors)) {
      coverage <- tryCatch(
        ltsa_candidate_reference_projection(
          candidate_vectors = candidate_vectors,
          reference_vectors = reference_vectors,
          nullvec = nullvec
        ),
        error = function(e) list(
          reference_projection_norms = NA_real_,
          reference_projection_min_norm = NA_real_,
          reference_candidate_space_rank = NA_integer_
        )
      )
      attempt$reference_projection_norms <-
        coverage$reference_projection_norms
      attempt$reference_projection_min_norm <-
        coverage$reference_projection_min_norm
      attempt$reference_candidate_space_rank <-
        coverage$reference_candidate_space_rank
      attempt$reference_space_rank <-
        coverage$reference_space_rank
      attempt$reference_overlap_singular_values <-
        coverage$reference_overlap_singular_values
      attempt$reference_overlap_min_singular_value <-
        coverage$reference_overlap_min_singular_value
    }
    if (!isTRUE(retain_candidate_spaces)) {
      attempt$.candidate_vectors <- NULL
    }
    attempt
  })

  res
}

ltsa_attempt_summary <- function(
  eig_k,
  eig_res,
  rr,
  acceptance,
  strict_rescue = FALSE,
  bank = FALSE,
  candidate_elapsed = NA_real_,
  ritz_elapsed = NA_real_,
  capture_candidate_space = FALSE
) {
  out <- list(
    eig_k = eig_k,
    backend = eig_res$backend %||% NA_character_,
    nconv = eig_res$nconv,
    niter = eig_res$niter %||% NA_integer_,
    nops = eig_res$nops %||% NA_integer_,
    mprod = eig_res$mprod %||% NA_integer_,
    ncv = ltsa_attempt_opt_integer(eig_res$opts, "ncv"),
    candidate_vector_rows = if (!is.null(eig_res$vectors)) {
      nrow(eig_res$vectors)
    } else {
      NA_integer_
    },
    candidate_vector_cols = if (!is.null(eig_res$vectors)) {
      ncol(eig_res$vectors)
    } else {
      NA_integer_
    },
    returned_columns = eig_res$returned_columns %||% ncol(eig_res$vectors),
    converged_columns = eig_res$converged_columns %||%
      eig_res$nconv %||% NA_integer_,
    convergence_known = eig_res$convergence_known %||% NA,
    candidate_elapsed = candidate_elapsed,
    ritz_elapsed = ritz_elapsed,
    rank_after_null = rr$rank_after_null,
    max_scaled_residual = acceptance$max_scaled_residual,
    global_gap = rr$global_gap,
    local_gap = rr$local_gap,
    zero_tol = rr$zero_tol,
    boundary_gap_relative = rr$boundary_gap_relative,
    near_zero_nonconstant_count = rr$near_zero_nonconstant_count,
    near_zero_nonconstant_counts = rr$near_zero_nonconstant_counts,
    partial_near_zero_block = ltsa_partial_near_zero_block(
      list(
        near_zero_nonconstant_count = rr$near_zero_nonconstant_count,
        acceptance = acceptance
      ),
      ndim = ncol(rr$vectors)
    ),
    first_ritz_values = rr$reported_ritz_values,
    gap_status = acceptance$gap_status,
    accepted = acceptance$ok,
    strict_rescue = strict_rescue,
    bank = bank
  )
  if (isTRUE(capture_candidate_space)) {
    out$.candidate_vectors <- eig_res$vectors
  }
  out
}

ltsa_attempt_error_summary <- function(
  eig_k,
  error,
  eig_res = NULL,
  strict_rescue = FALSE,
  bank = FALSE,
  candidate_elapsed = NA_real_,
  ritz_elapsed = NA_real_,
  capture_candidate_space = FALSE
) {
  out <- list(
    eig_k = eig_k,
    backend = eig_res$backend %||% NA_character_,
    nconv = eig_res$nconv %||% NA_integer_,
    niter = eig_res$niter %||% NA_integer_,
    nops = eig_res$nops %||% NA_integer_,
    mprod = eig_res$mprod %||% NA_integer_,
    ncv = ltsa_attempt_opt_integer(eig_res$opts, "ncv"),
    candidate_vector_rows = if (!is.null(eig_res$vectors)) {
      nrow(eig_res$vectors)
    } else {
      NA_integer_
    },
    candidate_vector_cols = if (!is.null(eig_res$vectors)) {
      ncol(eig_res$vectors)
    } else {
      NA_integer_
    },
    returned_columns = eig_res$returned_columns %||%
      if (!is.null(eig_res$vectors)) ncol(eig_res$vectors) else NA_integer_,
    converged_columns = eig_res$converged_columns %||%
      eig_res$nconv %||% NA_integer_,
    convergence_known = eig_res$convergence_known %||% NA,
    candidate_elapsed = candidate_elapsed,
    ritz_elapsed = ritz_elapsed,
    strict_rescue = strict_rescue,
    bank = bank,
    error = error
  )
  if (isTRUE(capture_candidate_space) && !is.null(eig_res$vectors)) {
    out$.candidate_vectors <- eig_res$vectors
  }
  out
}

ltsa_with_attempt <- function(eig_res, rr, acceptance, attempts, attempt) {
  ltsa_with_ritz(
    eig_res = eig_res,
    rr = rr,
    acceptance = acceptance,
    attempts = c(attempts, list(attempt))
  )
}

ltsa_bank_eig_result <- function(
  strict_eig_res,
  bank_vectors,
  lambda_max,
  bank_backend = NULL
) {
  bank_eig_res <- strict_eig_res
  bank_eig_res$vectors <- bank_vectors
  bank_eig_res$values <- ltsa_rayleigh_values(
    strict_eig_res$matrix,
    bank_vectors
  )
  bank_eig_res$backend <- bank_backend %||%
    paste0(strict_eig_res$backend %||% "candidate", "_bank")
  bank_eig_res$eig_k <- ncol(bank_vectors)
  bank_eig_res$niter <- NA_integer_
  bank_eig_res$nops <- NA_integer_
  bank_eig_res$mprod <- NA_integer_
  bank_eig_res$opts <- NULL
  bank_eig_res$returned_columns <- ncol(bank_vectors)
  bank_eig_res$converged_columns <- ncol(bank_vectors)
  bank_eig_res$nconv <- ncol(bank_vectors)
  bank_eig_res$convergence_known <- TRUE
  residuals <- ltsa_ritz_residuals(
    strict_eig_res$matrix,
    bank_vectors,
    bank_eig_res$values,
    lambda_max
  )
  bank_eig_res$absolute_residuals <- residuals$absolute_residuals
  bank_eig_res$scaled_residuals <- residuals$scaled_residuals
  bank_eig_res$residual_scale <- residuals$residual_scale
  bank_eig_res
}

ltsa_maybe_strict_rescue <- function(
  B,
  selected,
  ndim,
  nullvec,
  lambda_max,
  provider,
  provider_args,
  strict_rescue_arg_mapper,
  strict_rescue_controls,
  resid_tol,
  gap_tol,
  strict_rescue_extra,
  bank_backend = NULL,
  verbose = FALSE
) {
  if (!ltsa_strict_rescue_needed(selected, ndim)) {
    return(selected)
  }

  strict_eig_k <- ltsa_strict_rescue_eig_k(
    selected = selected,
    ndim = ndim,
    n = nrow(B),
    strict_rescue_extra = strict_rescue_extra
  )
  strict_args <- do.call(
    strict_rescue_arg_mapper,
    c(list(provider_args), strict_rescue_controls)
  )

  strict_candidate_elapsed <- NA_real_
  strict_candidate_timing <- system.time({
    strict_eig_res <- tryCatch(
      ltsa_call_candidate_provider(
        provider = provider,
        B = B,
        eig_k = strict_eig_k,
        provider_args = strict_args,
        lambda_max = lambda_max,
        verbose = verbose
      ),
      error = function(e) e
    )
  })
  strict_candidate_elapsed <- unname(strict_candidate_timing[["elapsed"]])
  if (inherits(strict_eig_res, "error")) {
    selected$attempts <- c(
      selected$attempts,
      list(ltsa_attempt_error_summary(
        eig_k = strict_eig_k,
        error = conditionMessage(strict_eig_res),
        strict_rescue = TRUE,
        candidate_elapsed = strict_candidate_elapsed
      ))
    )
    return(selected)
  }

  strict_ritz_elapsed <- NA_real_
  strict_ritz_timing <- system.time({
    strict_rr <- ltsa_ritz_select(
      B = strict_eig_res$matrix,
      vectors = strict_eig_res$vectors,
      ndim = ndim,
      nullvec = nullvec,
      lambda_max = strict_eig_res$lambda_max
    )
  })
  strict_ritz_elapsed <- unname(strict_ritz_timing[["elapsed"]])
  strict_acceptance <- accept_ltsa_ritz(
    rr = strict_rr,
    ndim = ndim,
    lambda_max = strict_eig_res$lambda_max,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    require_gap = FALSE
  )
  strict_attempt <- ltsa_attempt_summary(
    eig_k = strict_eig_k,
    eig_res = strict_eig_res,
    rr = strict_rr,
    acceptance = strict_acceptance,
    strict_rescue = TRUE,
    candidate_elapsed = strict_candidate_elapsed,
    ritz_elapsed = strict_ritz_elapsed,
    capture_candidate_space = TRUE
  )
  strict_result <- ltsa_with_attempt(
    eig_res = strict_eig_res,
    rr = strict_rr,
    acceptance = strict_acceptance,
    attempts = selected$attempts,
    attempt = strict_attempt
  )
  strict_result$acceptance$return_reason <- "strict_rescue"

  bank_vectors <- cbind(selected$candidate_vectors, strict_eig_res$vectors)
  bank_eig_res <- ltsa_bank_eig_result(
    strict_eig_res = strict_eig_res,
    bank_vectors = bank_vectors,
    lambda_max = strict_eig_res$lambda_max,
    bank_backend = bank_backend
  )
  bank_ritz_elapsed <- NA_real_
  bank_ritz_timing <- system.time({
    bank_rr <- ltsa_ritz_select(
      B = bank_eig_res$matrix,
      vectors = bank_vectors,
      ndim = ndim,
      nullvec = nullvec,
      lambda_max = strict_eig_res$lambda_max
    )
  })
  bank_ritz_elapsed <- unname(bank_ritz_timing[["elapsed"]])
  bank_acceptance <- accept_ltsa_ritz(
    rr = bank_rr,
    ndim = ndim,
    lambda_max = strict_eig_res$lambda_max,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    require_gap = FALSE
  )
  bank_attempt <- ltsa_attempt_summary(
    eig_k = ncol(bank_vectors),
    eig_res = bank_eig_res,
    rr = bank_rr,
    acceptance = bank_acceptance,
    strict_rescue = TRUE,
    bank = TRUE,
    ritz_elapsed = bank_ritz_elapsed,
    capture_candidate_space = TRUE
  )
  bank_result <- ltsa_with_attempt(
    eig_res = bank_eig_res,
    rr = bank_rr,
    acceptance = bank_acceptance,
    attempts = strict_result$attempts,
    attempt = bank_attempt
  )
  bank_result$acceptance$return_reason <- "strict_rescue_bank"

  rescued <- ltsa_rescue_candidate(strict_result, selected, ndim)
  rescued <- ltsa_rescue_candidate(bank_result, rescued, ndim)
  if (!identical(rescued, selected)) {
    rescued$acceptance$strict_rescue_used <- TRUE
    rescued$acceptance$strict_rescue_eig_k <- strict_eig_k
  }
  rescued
}

ltsa_with_ritz <- function(eig_res, rr, acceptance, attempts) {
  eig_res$candidate_vectors <- eig_res$vectors
  eig_res$candidate_values <- eig_res$values
  eig_res$vectors <- rr$vectors
  eig_res$values <- rr$values
  eig_res$ritz <- rr
  eig_res$acceptance <- acceptance
  eig_res$attempts <- attempts
  eig_res$absolute_residuals <- rr$absolute_residuals
  eig_res$scaled_residuals <- rr$scaled_residuals
  eig_res$residual_scale <- rr$residual_scale
  eig_res$rank_after_null <- rr$rank_after_null
  eig_res$boundary_gap <- rr$boundary_gap
  eig_res$global_gap <- rr$global_gap
  eig_res$local_gap <- rr$local_gap
  eig_res$zero_tol <- rr$zero_tol
  eig_res$boundary_gap_relative <- rr$boundary_gap_relative
  eig_res$near_zero_nonconstant_count <- rr$near_zero_nonconstant_count
  eig_res$near_zero_thresholds <- rr$near_zero_thresholds
  eig_res$near_zero_nonconstant_counts <- rr$near_zero_nonconstant_counts
  eig_res$reported_ritz_values <- rr$reported_ritz_values
  eig_res <- ltsa_mark_partial_near_zero_block(eig_res, ncol(rr$vectors))
  eig_res
}
