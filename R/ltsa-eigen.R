ltsa_iterative_search_k <- function(ndim, n) {
  base_eig_k <- ndim + 1L
  search_eig_k <- max(ndim + 3L, 2L * base_eig_k)
  min(n - 1L, search_eig_k)
}

symmetrize_ltsa_matrix <- function(B) {
  B <- 0.5 * (B + Matrix::t(B))

  if (methods::is(B, "sparseMatrix")) {
    B <- methods::as(B, "generalMatrix")
    B <- methods::as(B, "CsparseMatrix")
    B <- methods::as(B, "dMatrix")
    Matrix::drop0(B)
  } else {
    as.matrix(B)
  }
}

ltsa_matrix_max_abs <- function(B) {
  if (methods::is(B, "sparseMatrix")) {
    if (length(B@x) == 0L) {
      return(0)
    }
    return(max(abs(B@x)))
  }
  max(abs(B))
}

ltsa_matrix_is_effectively_zero <- function(B) {
  ltsa_matrix_max_abs(B) <= sqrt(.Machine$double.eps)
}

ltsa_default_ncv <- function(n, eig_k) {
  min(n, max(20L, 4L * eig_k + 10L))
}

ltsa_use_dense_eig <- function(n, eig_k, dense_n = 100L, dense_fraction = 0.5) {
  n <= dense_n || eig_k >= dense_fraction * n
}

ltsa_shift_margin <- function(lambda_max, shift_eps = 1e-6) {
  scale <- max(1, abs(lambda_max))
  max(shift_eps * scale, 100 * .Machine$double.eps * scale)
}

ltsa_shift_for_smallest <- function(B, shift) {
  shifted <- shift * Matrix::Diagonal(nrow(B)) - B
  if (methods::is(shifted, "sparseMatrix")) {
    shifted <- methods::as(shifted, "generalMatrix")
    shifted <- methods::as(shifted, "CsparseMatrix")
    shifted <- methods::as(shifted, "dMatrix")
    Matrix::drop0(shifted)
  } else {
    as.matrix(shifted)
  }
}

ltsa_normalize_columns <- function(V) {
  V <- as.matrix(V)
  norms <- sqrt(colSums(V * V))
  sweep(V, 2L, ifelse(norms == 0, 1, norms), "/")
}

ltsa_rayleigh_values <- function(B, vectors) {
  vectors <- as.matrix(vectors)
  BV <- as.matrix(B %*% vectors)
  denom <- colSums(vectors * vectors)
  colSums(vectors * BV) / denom
}

ltsa_ritz_residuals <- function(B, vectors, values, lambda_max) {
  vectors <- ltsa_normalize_columns(vectors)
  BV <- as.matrix(B %*% vectors)
  residual <- BV - sweep(vectors, 2L, values, "*")
  absolute_residual <- sqrt(colSums(residual * residual))
  residual_scale <- ltsa_residual_scale(lambda_max)

  list(
    absolute_residuals = absolute_residual,
    scaled_residuals = absolute_residual / residual_scale,
    residual_scale = residual_scale
  )
}

ltsa_residual_scale <- function(lambda_max) {
  if (
    is.null(lambda_max) ||
      length(lambda_max) < 1L ||
      !is.finite(lambda_max[[1L]])
  ) {
    return(1)
  }

  max(as.numeric(lambda_max[[1L]]), 1)
}

ltsa_default_null_vector <- function(n) {
  rep(1, n)
}

ltsa_normalize_null_vector <- function(nullvec, n) {
  if (length(nullvec) != n || any(!is.finite(nullvec))) {
    stop(
      "LTSA null vector must be finite and match the matrix dimension",
      call. = FALSE
    )
  }

  null_norm <- sqrt(sum(nullvec * nullvec))
  if (!is.finite(null_norm) || null_norm <= 0) {
    stop("LTSA null vector must have positive norm", call. = FALSE)
  }

  nullvec / null_norm
}

ltsa_near_zero_tol <- function(lambda_max, tol = sqrt(.Machine$double.eps)) {
  tol * ltsa_residual_scale(lambda_max)
}

ltsa_gap_zero_tol <- function(lambda_max, tol = sqrt(.Machine$double.eps)) {
  tol * ltsa_residual_scale(lambda_max)
}

ltsa_near_zero_threshold_scales <- function() {
  c(
    "1e-08" = 1e-8,
    "1e-07" = 1e-7,
    "1e-06" = 1e-6,
    "1e-05" = 1e-5
  )
}

ltsa_near_zero_thresholds <- function(lambda_max) {
  scales <- ltsa_near_zero_threshold_scales()
  stats::setNames(scales * ltsa_residual_scale(lambda_max), names(scales))
}

ltsa_near_zero_counts <- function(values, thresholds) {
  values <- as.numeric(values)
  stats::setNames(
    vapply(
      thresholds,
      function(threshold) {
        as.integer(sum(abs(values) <= threshold))
      },
      integer(1)
    ),
    names(thresholds)
  )
}

ltsa_format_ritz_values <- function(values, n = 10L) {
  values <- as.numeric(values)
  if (length(values) == 0L) {
    return("")
  }
  paste(signif(values[seq_len(min(n, length(values)))], 4), collapse = ", ")
}

ltsa_format_named_counts <- function(counts) {
  if (length(counts) == 0L) {
    return("")
  }
  paste(paste0(names(counts), "=", as.integer(counts)), collapse = ", ")
}

ltsa_validate_lambda_max <- function(lambda_max, B) {
  if (length(lambda_max) < 1L || !is.finite(lambda_max[[1L]])) {
    stop(
      "RSpectra largest-eigenvalue probe returned a non-finite value",
      call. = FALSE
    )
  }

  lambda_max <- max(0, as.numeric(lambda_max[[1L]]))
  if (lambda_max <= 0 && !ltsa_matrix_is_effectively_zero(B)) {
    stop(
      "RSpectra largest-eigenvalue probe returned a non-positive value for ",
      "a nonzero LTSA matrix",
      call. = FALSE
    )
  }
  lambda_max
}

ltsa_validate_lambda_probe <- function(probe, B) {
  nconv <- probe$nconv
  if (!is.null(nconv) && length(nconv) > 0L && nconv < 1L) {
    stop("RSpectra largest-eigenvalue probe did not converge", call. = FALSE)
  }

  ltsa_validate_lambda_max(probe$values, B)
}

ltsa_lambda_max_probe <- function(B, varargs) {
  opts <- rs_opt(eig_k = 1L, n = nrow(B))
  opts <- lmerge(opts, varargs)
  opts$retvec <- FALSE

  args <- list(
    A = B,
    k = 1L,
    which = "LA",
    opts = opts
  )
  probe <- do.call(RSpectra::eigs_sym, args)
  lambda_max <- ltsa_validate_lambda_probe(probe, B)

  list(
    value = lambda_max,
    nconv = probe$nconv %||% NA_integer_,
    niter = probe$niter %||% NA_integer_,
    nops = probe$nops %||% NA_integer_,
    which = "LA",
    opts = opts
  )
}

dense_ltsa_eig <- function(B, eig_k, backend = "dense_eigen") {
  dense <- as.matrix(B)
  eig <- eigen(dense, symmetric = TRUE)
  ord <- order(eig$values)
  values_all <- eig$values[ord]
  vectors_all <- eig$vectors[, ord, drop = FALSE]
  lambda_max <- max(values_all)

  take <- seq_len(eig_k)
  values <- values_all[take]
  vectors <- vectors_all[, take, drop = FALSE]
  residuals <- ltsa_ritz_residuals(B, vectors, values, lambda_max)

  list(
    vectors = vectors,
    values = values,
    nconv = eig_k,
    niter = NA_integer_,
    nops = NA_integer_,
    returned_columns = ncol(vectors),
    converged_columns = eig_k,
    absolute_residuals = residuals$absolute_residuals,
    scaled_residuals = residuals$scaled_residuals,
    residual_scale = residuals$residual_scale,
    backend = backend,
    lambda_max = lambda_max,
    lambda_probe = NULL,
    shift = NA_real_,
    shift_eps = NA_real_,
    shift_policy = "dense",
    solve_which = NA_character_,
    eig_k = eig_k,
    matrix = B
  )
}

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

ltsa_strict_rescue_needed <- function(res, ndim) {
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

ltsa_rs_eig_call <- function(
  B,
  eig_k,
  varargs,
  lambda_max = NULL,
  verbose = FALSE
) {
  do.call(
    rs_eig,
    c(
      list(X = B, k = eig_k, lambda_max = lambda_max, verbose = verbose),
      varargs
    )
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

ltsa_attempt_summary <- function(
  eig_k,
  eig_res,
  rr,
  acceptance,
  strict_rescue = FALSE,
  bank = FALSE
) {
  list(
    eig_k = eig_k,
    nconv = eig_res$nconv,
    rank_after_null = rr$rank_after_null,
    max_scaled_residual = acceptance$max_scaled_residual,
    global_gap = rr$global_gap,
    local_gap = rr$local_gap,
    zero_tol = rr$zero_tol,
    boundary_gap_relative = rr$boundary_gap_relative,
    near_zero_nonconstant_count = rr$near_zero_nonconstant_count,
    near_zero_nonconstant_counts = rr$near_zero_nonconstant_counts,
    first_ritz_values = rr$reported_ritz_values,
    gap_status = acceptance$gap_status,
    accepted = acceptance$ok,
    strict_rescue = strict_rescue,
    bank = bank
  )
}

ltsa_with_attempt <- function(eig_res, rr, acceptance, attempts, attempt) {
  ltsa_with_ritz(
    eig_res = eig_res,
    rr = rr,
    acceptance = acceptance,
    attempts = c(attempts, list(attempt))
  )
}

ltsa_bank_eig_result <- function(strict_eig_res, bank_vectors, lambda_max) {
  bank_eig_res <- strict_eig_res
  bank_eig_res$vectors <- bank_vectors
  bank_eig_res$values <- ltsa_rayleigh_values(
    strict_eig_res$matrix,
    bank_vectors
  )
  bank_eig_res$backend <- "rspectra_bank"
  bank_eig_res$eig_k <- ncol(bank_vectors)
  bank_eig_res$returned_columns <- ncol(bank_vectors)
  bank_eig_res$converged_columns <- ncol(bank_vectors)
  bank_eig_res$nconv <- ncol(bank_vectors)
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
  varargs,
  resid_tol,
  gap_tol,
  strict_rescue_tol,
  strict_rescue_maxitr,
  strict_rescue_extra,
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
  strict_args <- ltsa_strict_rescue_args(
    varargs,
    strict_rescue_tol = strict_rescue_tol,
    strict_rescue_maxitr = strict_rescue_maxitr
  )

  strict_eig_res <- tryCatch(
    ltsa_rs_eig_call(
      B = B,
      eig_k = strict_eig_k,
      varargs = strict_args,
      lambda_max = lambda_max,
      verbose = verbose
    ),
    error = function(e) e
  )
  if (inherits(strict_eig_res, "error")) {
    selected$attempts <- c(
      selected$attempts,
      list(list(
        eig_k = strict_eig_k,
        error = conditionMessage(strict_eig_res),
        strict_rescue = TRUE
      ))
    )
    return(selected)
  }

  strict_rr <- ltsa_ritz_select(
    B = strict_eig_res$matrix,
    vectors = strict_eig_res$vectors,
    ndim = ndim,
    nullvec = nullvec,
    lambda_max = strict_eig_res$lambda_max
  )
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
    strict_rescue = TRUE
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
    lambda_max = strict_eig_res$lambda_max
  )
  bank_rr <- ltsa_ritz_select(
    B = bank_eig_res$matrix,
    vectors = bank_vectors,
    ndim = ndim,
    nullvec = nullvec,
    lambda_max = strict_eig_res$lambda_max
  )
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
    bank = TRUE
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
  eig_res
}

ltsa_rspectra_ritz_eig <- function(
  B,
  ndim,
  ...,
  nullvec = ltsa_default_null_vector(nrow(B)),
  initial_extra = 4L,
  max_extra = 40L,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  gap_expansion_steps = 2L,
  strict_rescue = TRUE,
  strict_rescue_tol = 1e-10,
  strict_rescue_maxitr = 5000L,
  strict_rescue_extra = 5L,
  verbose = FALSE
) {
  n <- nrow(B)
  varargs <- list(...)
  eig_k <- ltsa_initial_ritz_candidate_k(ndim, n, initial_extra)
  max_k <- ltsa_max_ritz_candidate_k(ndim, n, max_extra)
  gap_expansion_steps <- as.integer(gap_expansion_steps)
  if (
    length(gap_expansion_steps) != 1L ||
      is.na(gap_expansion_steps) ||
      gap_expansion_steps < 0L
  ) {
    stop("gap_expansion_steps must be a non-negative integer", call. = FALSE)
  }
  strict_rescue <- as.logical(strict_rescue)
  if (length(strict_rescue) != 1L || is.na(strict_rescue)) {
    stop("strict_rescue must be TRUE or FALSE", call. = FALSE)
  }
  strict_rescue_maxitr <- as.integer(strict_rescue_maxitr)
  if (
    length(strict_rescue_maxitr) != 1L ||
      is.na(strict_rescue_maxitr) ||
      strict_rescue_maxitr < 1L
  ) {
    stop("strict_rescue_maxitr must be a positive integer", call. = FALSE)
  }
  strict_rescue_extra <- as.integer(strict_rescue_extra)
  if (
    length(strict_rescue_extra) != 1L ||
      is.na(strict_rescue_extra) ||
      strict_rescue_extra < 1L
  ) {
    stop("strict_rescue_extra must be a positive integer", call. = FALSE)
  }
  if (
    length(strict_rescue_tol) != 1L ||
      is.na(strict_rescue_tol) ||
      !is.finite(strict_rescue_tol) ||
      strict_rescue_tol <= 0
  ) {
    stop("strict_rescue_tol must be positive", call. = FALSE)
  }
  lambda_max <- NULL
  attempts <- list()
  best <- NULL
  first_gap_good <- NULL
  last_error <- NULL
  gap_expansions <- 0L

  repeat {
    eig_res <- tryCatch(
      ltsa_rs_eig_call(
        B = B,
        eig_k = eig_k,
        varargs = varargs,
        lambda_max = lambda_max,
        verbose = verbose
      ),
      error = function(e) e
    )
    if (inherits(eig_res, "error")) {
      last_error <- conditionMessage(eig_res)
      attempts[[length(attempts) + 1L]] <- list(
        eig_k = eig_k,
        error = last_error
      )
    } else {
      lambda_max <- eig_res$lambda_max
      rr <- tryCatch(
        ltsa_ritz_select(
          B = eig_res$matrix,
          vectors = eig_res$vectors,
          ndim = ndim,
          nullvec = nullvec,
          lambda_max = lambda_max
        ),
        error = function(e) e
      )

      if (inherits(rr, "error")) {
        last_error <- conditionMessage(rr)
        attempts[[length(attempts) + 1L]] <- list(
          eig_k = eig_k,
          nconv = eig_res$nconv,
          error = last_error
        )
      } else {
        acceptance <- accept_ltsa_ritz(
          rr = rr,
          ndim = ndim,
          lambda_max = lambda_max,
          resid_tol = resid_tol,
          gap_tol = gap_tol,
          require_gap = FALSE
        )
        attempts[[length(attempts) + 1L]] <- ltsa_attempt_summary(
          eig_k = eig_k,
          eig_res = eig_res,
          rr = rr,
          acceptance = acceptance
        )
        tsmessage(
          "Rayleigh-Ritz LTSA postprocess: post-null rank = ",
          rr$rank_after_null,
          ", max selected scaled residual = ",
          signif(acceptance$max_scaled_residual, 4),
          ", global boundary gap = ",
          signif(rr$global_gap, 4),
          ", local boundary gap = ",
          signif(rr$local_gap, 4),
          ", zero_tol = ",
          signif(rr$zero_tol, 4),
          ", near-zero nonconstant modes = ",
          rr$near_zero_nonconstant_count,
          ", threshold counts = [",
          ltsa_format_named_counts(rr$near_zero_nonconstant_counts),
          "], first Ritz values = [",
          ltsa_format_ritz_values(rr$reported_ritz_values),
          "]"
        )
        best <- ltsa_with_ritz(eig_res, rr, acceptance, attempts)

        if (acceptance$rank_ok && acceptance$resid_ok && acceptance$gap_ok) {
          best$acceptance$return_reason <- "residual_rank_gap_ok"
          if (strict_rescue) {
            best <- ltsa_maybe_strict_rescue(
              B = B,
              selected = best,
              ndim = ndim,
              nullvec = nullvec,
              lambda_max = lambda_max,
              varargs = varargs,
              resid_tol = resid_tol,
              gap_tol = gap_tol,
              strict_rescue_tol = strict_rescue_tol,
              strict_rescue_maxitr = strict_rescue_maxitr,
              strict_rescue_extra = strict_rescue_extra,
              verbose = verbose
            )
          }
          return(best)
        }
        if (acceptance$gap_issue_only) {
          if (is.null(first_gap_good)) {
            first_gap_good <- best
            first_gap_good$acceptance$return_reason <- "weak_first_residual_rank_good"
            first_gap_good$acceptance$first_weak_gap_eig_k <- eig_k
          }
          if (gap_expansions >= gap_expansion_steps || eig_k >= max_k) {
            selected <- first_gap_good %||% best
            selected$attempts <- attempts
            selected$acceptance$return_reason <- "weak_first_residual_rank_good"
            selected$acceptance$gap_expansions <- gap_expansions
            selected$acceptance$diagnostic_final_eig_k <- eig_k
            if (strict_rescue) {
              selected <- ltsa_maybe_strict_rescue(
                B = B,
                selected = selected,
                ndim = ndim,
                nullvec = nullvec,
                lambda_max = lambda_max,
                varargs = varargs,
                resid_tol = resid_tol,
                gap_tol = gap_tol,
                strict_rescue_tol = strict_rescue_tol,
                strict_rescue_maxitr = strict_rescue_maxitr,
                strict_rescue_extra = strict_rescue_extra,
                verbose = verbose
              )
            }
            if (isTRUE(selected$acceptance$strict_rescue_used)) {
              tsmessage(
                "LTSA strict rescue found lower-energy residual-good Ritz ",
                "vectors after an incomplete near-zero selected block; ",
                "returning strict-rescue result with selected values = [",
                ltsa_format_ritz_values(selected$values),
                "]"
              )
            } else {
              tsmessage(
                "LTSA Ritz boundary gap ",
                if (acceptance$gap_available) "remains below tolerance" else
                  "is unavailable",
                " after diagnostic expansion to ",
                eig_k,
                " candidates; returning residual-good, rank-good Ritz vectors."
              )
            }
            return(selected)
          }

          gap_expansions <- gap_expansions + 1L
          attempts[[length(attempts)]]$weak_gap_expansion <- gap_expansions
        }
      }
    }

    if (eig_k >= max_k) {
      break
    }
    eig_k <- ltsa_next_ritz_candidate_k(eig_k, max_k)
  }

  if (is.null(best)) {
    stop(
      "LTSA Rayleigh-Ritz selection failed after requesting up to ",
      max_k,
      " candidate vectors: ",
      last_error %||% "no usable candidate subspace",
      call. = FALSE
    )
  }

  if (!best$acceptance$rank_ok) {
    stop(
      "LTSA Rayleigh-Ritz selection did not find enough post-null rank after ",
      "requesting up to ",
      best$eig_k,
      " candidate vectors",
      call. = FALSE
    )
  }
  if (strict_rescue) {
    best <- ltsa_maybe_strict_rescue(
      B = B,
      selected = best,
      ndim = ndim,
      nullvec = nullvec,
      lambda_max = best$lambda_max %||% lambda_max,
      varargs = varargs,
      resid_tol = resid_tol,
      gap_tol = gap_tol,
      strict_rescue_tol = strict_rescue_tol,
      strict_rescue_maxitr = strict_rescue_maxitr,
      strict_rescue_extra = strict_rescue_extra,
      verbose = verbose
    )
  }
  if (!best$acceptance$resid_ok) {
    warning(
      "LTSA Rayleigh-Ritz residuals remain above tolerance after requesting ",
      best$eig_k,
      " candidate vectors: max scaled residual = ",
      signif(best$acceptance$max_scaled_residual, 4),
      ", tolerance = ",
      signif(resid_tol, 4),
      call. = FALSE
    )
  } else if (!best$acceptance$gap_ok) {
    tsmessage(
      "LTSA Ritz boundary gap ",
      if (best$acceptance$gap_available) "remains below tolerance" else
        "is unavailable",
      " after requesting ",
      best$eig_k,
      " candidates: global gap = ",
      signif(best$acceptance$rel_gap, 4),
      ", local gap = ",
      signif(best$acceptance$local_gap, 4),
      ", tolerance = ",
      signif(gap_tol, 4),
      "; returning residual-good, rank-good Ritz vectors."
    )
  }

  best
}

rs_opt <- function(eig_k = NULL, n = NULL) {
  opts <- list(tol = 1e-6)
  if (!is.null(eig_k) && !is.null(n)) {
    opts$ncv <- ltsa_default_ncv(n, eig_k)
  }
  opts
}

# Avoid shift-invert here; near-zero LTSA/Laplacian eigenvalues can make sparse
# factorizations hang or skip eigenvectors.
# https://github.com/yixuan/spectra/issues/126
rs_eig <-
  function(
    X,
    k = ncol(X) - 1,
    ...,
    lambda_max = NULL,
    verbose = FALSE,
    shift_eps = 1e-6,
    dense_n = 100L,
    dense_fraction = 0.5
  ) {
    varargs <- list(...)
    B <- symmetrize_ltsa_matrix(X)
    n <- ncol(B)
    eig_k <- as.integer(k)
    if (length(eig_k) != 1L || is.na(eig_k) || eig_k < 1L || eig_k >= n) {
      stop(
        "eig_k must be a positive integer less than the matrix dimension",
        call. = FALSE
      )
    }

    if (
      length(varargs) == 0L &&
        ltsa_use_dense_eig(
          n,
          eig_k,
          dense_n = dense_n,
          dense_fraction = dense_fraction
        )
    ) {
      tsmessage("Using dense eigenvalue decomposition")
      return(dense_ltsa_eig(B, eig_k))
    }

    lambda_probe <- NULL
    if (is.null(lambda_max)) {
      tsmessage("Finding largest eigenvalue")
      lambda_probe <- ltsa_lambda_max_probe(B, varargs)
      lambda_max <- lambda_probe$value
    } else {
      lambda_max <- ltsa_validate_lambda_max(lambda_max, B)
    }
    shift <- lambda_max + ltsa_shift_margin(lambda_max, shift_eps)
    X_shift <- ltsa_shift_for_smallest(B, shift)

    tsmessage("Decomposing shifted matrix")
    opts <- rs_opt(eig_k = eig_k, n = n)
    opts <- lmerge(opts, varargs)
    args <-
      list(
        A = X_shift,
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
    if (is.null(res$vectors) || ncol(res$vectors) < eig_k) {
      stop(
        "RSpectra returned fewer LTSA candidate vectors than requested",
        call. = FALSE
      )
    }

    vectors <- as.matrix(res$vectors[, seq_len(eig_k), drop = FALSE])
    values <- ltsa_rayleigh_values(B, vectors)
    ord <- order(values)
    values <- values[ord]
    vectors <- vectors[, ord, drop = FALSE]
    shifted_values <- res$values
    if (!is.null(shifted_values) && length(shifted_values) >= eig_k) {
      shifted_values <- shifted_values[seq_len(eig_k)][ord]
    }
    residuals <- ltsa_ritz_residuals(B, vectors, values, lambda_max)
    tsmessage(
      "RSpectra converged ",
      ifelse(is.na(nconv), eig_k, nconv),
      " / ",
      eig_k,
      " LTSA candidate vectors; max scaled residual = ",
      signif(max(residuals$scaled_residuals), 4)
    )

    list(
      vectors = vectors,
      values = values,
      shifted_values = shifted_values,
      nconv = nconv,
      niter = res$niter %||% NA_integer_,
      nops = res$nops %||% NA_integer_,
      returned_columns = ncol(res$vectors),
      converged_columns = ifelse(is.na(nconv), eig_k, nconv),
      absolute_residuals = residuals$absolute_residuals,
      scaled_residuals = residuals$scaled_residuals,
      residual_scale = residuals$residual_scale,
      backend = "rspectra",
      lambda_max = lambda_max,
      lambda_probe = lambda_probe,
      shift = shift,
      shift_eps = shift_eps,
      shift_policy = "lambda_max_plus_margin",
      solve_which = "LA",
      eig_k = eig_k,
      opts = opts,
      matrix = B
    )
  }

irlba_eig <-
  function(X, k = ncol(X) - 1, ..., lambda_max = NULL, verbose = FALSE) {
    varargs <- list(...)

    if (is.null(lambda_max)) {
      tsmessage("Finding largest eigenvalue")
      args1 <- list(
        A = X,
        nv = 1,
        nu = 0
      )
      lambda_max <- do.call(irlba::irlba, args1)$d
    }
    lm2 <- 2.0 * lambda_max
    X_shift <- shift_lap(X, lm2)

    tsmessage("Decomposing shifted matrix")
    args <- lmerge(
      list(
        A = X_shift,
        nv = k,
        nu = 0
      ),
      varargs
    )
    res <- do.call(irlba::irlba, args)
    res$v
  }

svdr_eig <- function(
  X,
  k = ncol(X) - 1,
  ...,
  lambda_max = NULL,
  verbose = FALSE
) {
  varargs <- list(...)

  if (is.null(lambda_max)) {
    tsmessage("Finding largest eigenvalue")
    args1 <- list(
      x = X,
      k = 1
    )
    lambda_max <- do.call(irlba::svdr, args1)$d
  }
  lm2 <- 2.0 * lambda_max
  X_shift <- shift_lap(X, lm2)

  tsmessage("Decomposing shifted matrix")
  args <- lmerge(
    list(
      x = X_shift,
      k = k
    ),
    varargs
  )
  res <- do.call(irlba::svdr, args)
  res$v
}
