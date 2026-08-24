# Null-aware Rayleigh-Ritz selection and diagnostics for LTSA eigenanalysis.

# Wrapper-level Rayleigh-Ritz extraction. Iterative solvers commonly do their
# own Ritz extraction internally, but this step answers the LTSA-specific
# question: within the returned candidate span, which vectors best diagonalize
# the original LTSA matrix after the known constant null direction is removed?
scaled_ritz_boundary_gap <- function(values, boundary_dim, lambda_max) {
  if (length(values) <= boundary_dim) {
    return(NA_real_)
  }
  boundary_gap <- values[[boundary_dim + 1L]] - values[[boundary_dim]]
  boundary_gap / residual_scale(lambda_max)
}

select_ritz_vectors <- function(
  B,
  vectors,
  ndim,
  null_vector = default_null_vector(nrow(B)),
  lambda_max = NULL,
  drop_tol = 1e-8,
  rank_tol = 1e-10,
  near_zero_tol = near_zero_tolerance(lambda_max)
) {
  vectors <- as.matrix(vectors)
  n <- nrow(B)
  if (
    !is.numeric(ndim) ||
      length(ndim) != 1L ||
      is.na(ndim) ||
      !is.finite(ndim) ||
      ndim < 1L ||
      ndim != floor(ndim)
  ) {
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

  # Remove the known null direction and rank-check the remaining span.
  null_vector <- normalize_null_vector(null_vector, n)
  projected <- vectors - null_vector %*% crossprod(null_vector, vectors)
  original_norms <- sqrt(colSums(vectors * vectors))
  projected_norms <- sqrt(colSums(projected * projected))
  keep <- projected_norms > drop_tol * pmax(1, original_norms)
  projected <- projected[, keep, drop = FALSE]

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

  # Rayleigh-Ritz extraction inside the projected candidate span.
  Q <- qr.Q(qrp)[, seq_len(rank_after_null), drop = FALSE]
  BQ <- as.matrix(B %*% Q)
  H <- as.matrix(crossprod(Q, BQ))
  H <- 0.5 * (H + t(H))

  eig <- eigen(H, symmetric = TRUE)
  ord <- order(eig$values)
  values_all <- eig$values[ord]
  coef_all <- eig$vectors[, ord, drop = FALSE]
  vectors_all <- Q %*% coef_all

  # Residual, gap, and near-zero diagnostics
  residuals_all <- ritz_residuals(B, vectors_all, values_all, lambda_max)

  global_gap <- scaled_ritz_boundary_gap(values_all, ndim, lambda_max)

  take <- seq_len(ndim)
  near_zero_nonconstant_count <- sum(abs(values_all) <= near_zero_tol)

  list(
    vectors = vectors_all[, take, drop = FALSE],
    values = values_all[take],
    ritz_values = values_all,
    residuals = residuals_all$scaled_residuals[take],
    rank_after_null = rank_after_null,
    global_gap = global_gap,
    near_zero_nonconstant_count = near_zero_nonconstant_count
  )
}

prefix_ritz_result <- function(ritz_result, ndim, lambda_max) {
  take <- seq_len(ndim)
  list(
    vectors = ritz_result$vectors[, take, drop = FALSE],
    values = ritz_result$values[take],
    ritz_values = ritz_result$ritz_values,
    residuals = ritz_result$residuals[take],
    rank_after_null = ritz_result$rank_after_null,
    global_gap = scaled_ritz_boundary_gap(
      ritz_result$ritz_values,
      ndim,
      lambda_max
    ),
    near_zero_nonconstant_count = ritz_result$near_zero_nonconstant_count
  )
}

backend_metadata <- function(candidate_result) {
  backend <- candidate_result$backend %||% "unknown"
  backend <- as.character(backend[[1L]])
  exact_dense <- backend %in%
    c(
      "dense_eigen",
      "eig",
      "eigen"
    )
  out <- list(
    name = backend,
    convergence_known = isTRUE(candidate_result$convergence_known) ||
      exact_dense
  )

  if (identical(backend, "rspectra")) {
    out$nconv <- as.integer(candidate_result$nconv %||% NA_integer_)
    out$niter <- as.integer(candidate_result$niter %||% NA_integer_)
    out$nops <- as.integer(candidate_result$nops %||% NA_integer_)
  } else {
    niter <- as.integer(candidate_result$niter %||% NA_integer_)
    mprod <- as.integer(candidate_result$mprod %||% NA_integer_)
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

backend_failure_messages <- function(candidate_result, eig_k) {
  backend <- candidate_result$backend %||% "unknown"
  backend <- as.character(backend[[1L]])
  nconv <- as.integer(candidate_result$nconv %||% NA_integer_)
  converged_columns <- as.integer(
    candidate_result$converged_columns %||% nconv
  )

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
    isTRUE(candidate_result$convergence_known) &&
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

diagnose_ritz <- function(
  candidate_result,
  ritz_result,
  eig_k,
  ndim,
  resid_tol,
  gap_tol,
  dimension_name = "ndim",
  block_name = "requested embedding",
  boundary_name = "selected block"
) {
  backend <- backend_metadata(candidate_result)
  lambda_max <- candidate_result$lambda_max %||% NA_real_
  values <- as.numeric(ritz_result$values)
  residuals <- as.numeric(ritz_result$residuals)
  max_residual <- if (length(residuals) == 0L) {
    Inf
  } else {
    max(residuals)
  }

  invalid_messages <- character()
  if (length(values) < ndim || ncol(ritz_result$vectors) < ndim) {
    invalid_messages <- c(
      invalid_messages,
      paste0(
        "Fewer than ",
        dimension_name,
        " usable nonconstant Ritz vectors were selected."
      )
    )
  }
  if (ritz_result$rank_after_null < ndim) {
    invalid_messages <- c(
      invalid_messages,
      paste0(
        "Post-null candidate rank ",
        ritz_result$rank_after_null,
        " is less than ",
        dimension_name,
        " = ",
        ndim,
        "."
      )
    )
  }
  if (!all(is.finite(values)) || !all(is.finite(ritz_result$vectors))) {
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
    backend_failure_messages(candidate_result, eig_k)
  )

  warning_messages <- character()
  if (length(ritz_result$ritz_values) < ndim + 1L) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "Candidate span contains fewer than ",
        dimension_name,
        " + 1 post-null Ritz ",
        "values; no spare boundary direction is available."
      )
    )
  } else if (!is.finite(ritz_result$global_gap)) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "Ritz boundary gap after the ",
        boundary_name,
        " is unavailable."
      )
    )
  } else if (ritz_result$global_gap < gap_tol) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "Weak Ritz boundary gap after the ",
        boundary_name,
        ": ",
        signif(ritz_result$global_gap, 4),
        " < ",
        signif(gap_tol, 4),
        "."
      )
    )
  }

  # The observed candidate span shows truncation when more near-zero modes
  # are present than the requested embedding dimensions. This is an
  # identifiability warning, not evidence of clustering.
  near_zero_block_truncated <-
    ritz_result$near_zero_nonconstant_count > ndim
  if (near_zero_block_truncated) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "The ",
        block_name,
        " cuts through a larger near-zero eigenspace; ",
        "individual coordinates are not uniquely identifiable. Increase ",
        "eig_k to inspect more of the observed block, or request a larger ",
        "spectral_dim with output = \"result\" to retain more modes from the ",
        "same operator. A stricter tolerance does not resolve a genuine ",
        "repeated eigenspace."
      )
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

  list(
    method = backend$name,
    eig_k = eig_k,
    values = values,
    ritz_values = as.numeric(ritz_result$ritz_values),
    residuals = residuals,
    rank = ritz_result$rank_after_null,
    lambda_max = lambda_max,
    status = status,
    messages = messages,
    backend = backend,
    diagnostics = list(
      near_zero_block_truncated = isTRUE(near_zero_block_truncated)
    )
  )
}

validate_eigen_controls <- function(args, eig_method, output) {
  if (is.null(args)) {
    args <- list()
  }
  if (length(args) == 0L) {
    return(
      list(
        provider_args = list(),
        resid_tol = 1e-5,
        gap_tol = 1e-4
      )
    )
  }

  if (identical(output, "B")) {
    stop(
      "output = \"B\" does not accept arguments in `...`",
      call. = FALSE
    )
  }

  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep.int("", length(args))
  }
  if (any(is.na(arg_names) | !nzchar(arg_names))) {
    stop(
      "Every argument in `...` must be named",
      call. = FALSE
    )
  }
  if (anyDuplicated(arg_names)) {
    duplicate_name <- arg_names[duplicated(arg_names)][[1L]]
    stop(
      "Argument `",
      duplicate_name,
      "` is supplied more than once",
      call. = FALSE
    )
  }

  diagnostic_names <- c("resid_tol", "gap_tol")
  method_names <- list(
    auto = c("dense_n", "dense_fraction"),
    rspectra = c("tol", "maxitr", "ncv", "shift_eps"),
    irlba = c("tol", "maxit", "reorth", "shift_eps"),
    svdr = c("tol", "it", "extra", "shift_eps"),
    eig = character()
  )
  supported_names <- unique(c(
    diagnostic_names,
    method_names[[eig_method]]
  ))
  unsupported_names <- setdiff(arg_names, supported_names)
  if (length(unsupported_names) > 0L) {
    stop(
      "Argument `",
      unsupported_names[[1L]],
      "` is not supported for eig_method = \"",
      eig_method,
      "\"",
      call. = FALSE
    )
  }

  resid_tol <- args$resid_tol %||% 1e-5
  gap_tol <- args$gap_tol %||% 1e-4
  resid_tol <- check_positive_finite(resid_tol, "resid_tol")
  gap_tol <- check_positive_finite(gap_tol, "gap_tol")

  if ("shift_eps" %in% arg_names) {
    args$shift_eps <- check_positive_finite(args$shift_eps, "shift_eps")
  }
  if ("dense_n" %in% arg_names) {
    args$dense_n <- check_whole_number(args$dense_n, "dense_n", min = 0)
  }
  if ("dense_fraction" %in% arg_names) {
    args$dense_fraction <- check_positive_finite(
      args$dense_fraction,
      "dense_fraction"
    )
    if (args$dense_fraction > 1) {
      stop("dense_fraction must be <= 1", call. = FALSE)
    }
  }

  provider_args <- args[!(arg_names %in% diagnostic_names)]

  list(
    provider_args = provider_args,
    resid_tol = resid_tol,
    gap_tol = gap_tol
  )
}

check_positive_finite <- function(x, name) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0
  ) {
    stop(name, " must be a finite positive number", call. = FALSE)
  }
  as.numeric(x)
}

run_eigenanalysis <- function(
  B,
  ndim,
  eig_method,
  eig_k,
  eigen_args,
  null_vector,
  verbose,
  spectral_dim = ndim
) {
  provider <- switch(
    eig_method,
    auto = auto_candidate_provider,
    rspectra = {
      report_progress("Calling rspectra", verbose = verbose)
      rspectra_candidate_provider
    },
    irlba = {
      report_progress("Calling irlba", verbose = verbose)
      irlba_candidate_provider
    },
    svdr = {
      report_progress("Calling irlba svdr", verbose = verbose)
      svdr_candidate_provider
    },
    eig = {
      report_progress("Using full eigenvalue decomposition", verbose = verbose)
      dense_candidate_provider
    }
  )

  provider_args <- eigen_args$provider_args

  solve_with_ritz(
    B = B,
    ndim = ndim,
    provider = provider,
    provider_args = provider_args,
    null_vector = null_vector,
    eig_k = eig_k,
    resid_tol = eigen_args$resid_tol,
    gap_tol = eigen_args$gap_tol,
    verbose = verbose,
    spectral_dim = spectral_dim
  )
}

solve_with_ritz <- function(
  B,
  ndim,
  provider,
  provider_args = list(),
  null_vector = default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE,
  spectral_dim = ndim
) {
  eig_k <- validate_spectral_eig_k(
    eig_k = eig_k,
    ndim = ndim,
    spectral_dim = spectral_dim,
    n = nrow(B)
  )
  B <- symmetrize_operator(B)
  provider_args <- provider_args %||% list()
  candidate_result <- do.call(
    provider,
    c(
      list(
        B = B,
        eig_k = eig_k,
        lambda_max = NULL,
        verbose = verbose
      ),
      provider_args
    )
  )
  lambda_max <- candidate_result$lambda_max %||% NA_real_
  ritz_result <- select_ritz_vectors(
    B = B,
    vectors = candidate_result$vectors,
    ndim = spectral_dim,
    null_vector = null_vector,
    lambda_max = lambda_max
  )
  displayed_ritz_result <- if (spectral_dim > ndim) {
    prefix_ritz_result(
      ritz_result = ritz_result,
      ndim = ndim,
      lambda_max = lambda_max
    )
  } else {
    ritz_result
  }
  eigen <- diagnose_ritz(
    candidate_result = candidate_result,
    ritz_result = displayed_ritz_result,
    eig_k = eig_k,
    ndim = ndim,
    resid_tol = resid_tol,
    gap_tol = gap_tol
  )

  result <- list(
    vectors = ritz_result$vectors,
    values = ritz_result$values,
    eigen = eigen,
    backend = eigen$backend,
    lambda_max = eigen$lambda_max,
    eig_k = eig_k
  )
  if (spectral_dim > ndim) {
    result$spectral_eigen <- diagnose_ritz(
      candidate_result = candidate_result,
      ritz_result = ritz_result,
      eig_k = eig_k,
      ndim = spectral_dim,
      resid_tol = resid_tol,
      gap_tol = gap_tol,
      dimension_name = "spectral_dim",
      block_name = "retained spectral block",
      boundary_name = "retained block"
    )
    result$displayed_boundary_gap <- displayed_ritz_result$global_gap
    result$spectral_boundary_gap <- ritz_result$global_gap
  }
  result
}
