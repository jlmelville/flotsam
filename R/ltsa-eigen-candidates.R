# Backend-neutral candidate result objects and candidate providers for
# RSpectra, irlba, svdr, and dense LTSA diagnostic solves.

# Backend-neutral candidate result contract. Providers may use different
# eigensolver APIs, but the fixed-width LTSA path only needs candidate vectors,
# Rayleigh values against the matrix of interest, convergence metadata when
# available, and residual diagnostics.
ltsa_candidate_result <- function(
  vectors,
  values = NULL,
  shifted_values = NULL,
  backend,
  eig_k,
  matrix,
  lambda_max = NULL,
  lambda_probe = NULL,
  nconv = NA_integer_,
  niter = NA_integer_,
  nops = NA_integer_,
  mprod = NA_integer_,
  opts = NULL,
  convergence_known = FALSE,
  returned_columns = ncol(as.matrix(vectors)),
  converged_columns = ifelse(is.na(nconv), NA_integer_, nconv)
) {
  vectors <- as.matrix(vectors)
  if (is.null(values)) {
    values <- ltsa_rayleigh_values(matrix, vectors)
  }

  list(
    vectors = vectors,
    values = as.numeric(values),
    shifted_values = shifted_values,
    backend = backend,
    eig_k = as.integer(eig_k),
    matrix = matrix,
    lambda_max = lambda_max,
    lambda_probe = lambda_probe,
    nconv = nconv,
    niter = niter,
    nops = nops,
    mprod = mprod,
    opts = opts,
    convergence_known = isTRUE(convergence_known),
    returned_columns = as.integer(returned_columns),
    converged_columns = as.integer(converged_columns)
  )
}

ltsa_as_candidate_result <- function(
  res,
  B,
  eig_k,
  backend = res$backend %||% "unknown",
  lambda_max = NULL,
  convergence_known = FALSE
) {
  if (is.null(res$vectors)) {
    stop(
      "LTSA candidate provider did not return candidate vectors",
      call. = FALSE
    )
  }

  matrix <- res$matrix %||% B
  values <- res$values
  if (is.null(values)) {
    values <- ltsa_rayleigh_values(matrix, res$vectors)
  }

  candidate <- ltsa_candidate_result(
    vectors = res$vectors,
    values = values,
    shifted_values = if ("shifted_values" %in% names(res)) {
      res$shifted_values
    } else {
      NULL
    },
    backend = res$backend %||% backend,
    eig_k = res$eig_k %||% eig_k,
    matrix = matrix,
    lambda_max = res$lambda_max %||% lambda_max,
    lambda_probe = if ("lambda_probe" %in% names(res)) {
      res$lambda_probe
    } else {
      NULL
    },
    nconv = res$nconv %||% NA_integer_,
    niter = res$niter %||% NA_integer_,
    nops = res$nops %||% NA_integer_,
    mprod = res$mprod %||% NA_integer_,
    opts = res$opts %||% NULL,
    convergence_known = res$convergence_known %||% convergence_known,
    returned_columns = res$returned_columns %||% ncol(as.matrix(res$vectors)),
    converged_columns = res$converged_columns %||%
      (res$nconv %||% NA_integer_)
  )
  extra_names <- setdiff(names(res), names(candidate))
  candidate[extra_names] <- res[extra_names]
  candidate
}

ltsa_call_candidate_provider <- function(
  provider,
  B,
  eig_k,
  provider_args = list(),
  lambda_max = NULL,
  verbose = FALSE
) {
  if (!is.function(provider)) {
    stop("LTSA candidate provider must be a function", call. = FALSE)
  }
  if (is.null(provider_args)) {
    provider_args <- list()
  }
  if (!is.list(provider_args)) {
    stop("LTSA candidate provider arguments must be a list", call. = FALSE)
  }

  do.call(
    provider,
    c(
      list(
        B = B,
        eig_k = eig_k,
        lambda_max = lambda_max,
        verbose = verbose
      ),
      provider_args
    )
  )
}

# RSpectra candidate provider. The actual solve is delegated to rs_eig(), which
# uses the shifted largest-algebraic formulation and enforces RSpectra nconv.
ltsa_rspectra_candidate_provider <- function(
  B,
  eig_k,
  ...,
  lambda_max = NULL,
  verbose = FALSE
) {
  res <- do.call(
    rs_eig,
    c(
      list(X = B, k = eig_k, lambda_max = lambda_max, verbose = verbose),
      list(...)
    )
  )

  ltsa_as_candidate_result(
    res,
    B = B,
    eig_k = eig_k,
    backend = res$backend %||% "rspectra",
    lambda_max = lambda_max,
    convergence_known = TRUE
  )
}

ltsa_irlba_lambda_max_probe <- function(B) {
  args <- list(
    A = B,
    nv = 1L,
    nu = 0L
  )
  probe <- do.call(irlba::irlba, args)
  lambda_max <- ltsa_validate_backend_lambda_max(probe$d, B, backend = "irlba")

  list(
    value = lambda_max,
    niter = probe$iter %||% NA_integer_,
    mprod = probe$mprod %||% NA_integer_,
    opts = args
  )
}

ltsa_svdr_lambda_max_probe <- function(B) {
  args <- list(
    x = B,
    k = 1L
  )
  probe <- do.call(irlba::svdr, args)
  lambda_max <- ltsa_validate_backend_lambda_max(probe$d, B, backend = "svdr")

  list(
    value = lambda_max,
    mprod = probe$mprod %||% NA_integer_,
    opts = args
  )
}

ltsa_validate_candidate_eig_k <- function(eig_k, n) {
  eig_k <- as.integer(eig_k)
  if (length(eig_k) != 1L || is.na(eig_k) || eig_k < 1L || eig_k >= n) {
    stop(
      "eig_k must be a positive integer less than the matrix dimension",
      call. = FALSE
    )
  }
  eig_k
}

# Candidate results for shifted solves are always re-valued against the
# original LTSA matrix. The shifted backend values are retained only as backend
# metadata; final selection uses Rayleigh values of B.
ltsa_shifted_candidate_result <- function(
  B,
  vectors,
  eig_k,
  backend,
  lambda_max,
  lambda_probe,
  shift,
  shift_eps,
  shifted_values = NULL,
  niter = NA_integer_,
  nops = NA_integer_,
  mprod = NA_integer_,
  opts = NULL,
  returned_columns = ncol(as.matrix(vectors))
) {
  if (is.null(vectors) || ncol(vectors) < eig_k) {
    stop(
      backend,
      " returned fewer LTSA candidate vectors than requested",
      call. = FALSE
    )
  }

  vectors <- as.matrix(vectors[, seq_len(eig_k), drop = FALSE])
  values <- ltsa_rayleigh_values(B, vectors)
  ord <- order(values)
  values <- values[ord]
  vectors <- vectors[, ord, drop = FALSE]
  if (!is.null(shifted_values) && length(shifted_values) >= eig_k) {
    shifted_values <- shifted_values[seq_len(eig_k)][ord]
  }
  residuals <- ltsa_ritz_residuals(B, vectors, values, lambda_max)
  tsmessage(
    backend,
    " returned ",
    eig_k,
    " LTSA candidate vectors; max scaled residual = ",
    signif(max(residuals$scaled_residuals), 4)
  )

  candidate <- ltsa_candidate_result(
    vectors = vectors,
    values = values,
    shifted_values = shifted_values,
    backend = backend,
    eig_k = eig_k,
    matrix = B,
    lambda_max = lambda_max,
    lambda_probe = lambda_probe,
    nconv = NA_integer_,
    niter = niter,
    nops = nops,
    mprod = mprod,
    opts = opts,
    convergence_known = FALSE,
    returned_columns = returned_columns,
    converged_columns = NA_integer_
  )
  candidate$absolute_residuals <- residuals$absolute_residuals
  candidate$scaled_residuals <- residuals$scaled_residuals
  candidate$residual_scale <- residuals$residual_scale
  candidate$shift <- shift
  candidate$shift_eps <- shift_eps
  candidate$shift_policy <- "lambda_max_plus_margin"
  candidate$solve_which <- "largest_singular"
  candidate
}

# irlba candidate provider. It solves the shifted problem by requesting right
# singular vectors, then all final ordering is recomputed by Rayleigh values of
# the original LTSA matrix.
ltsa_irlba_candidate_provider <- function(
  B,
  eig_k,
  ...,
  lambda_max = NULL,
  verbose = FALSE,
  shift_eps = 1e-6,
  dense_n = 100L,
  dense_fraction = 0.5
) {
  varargs <- list(...)
  B <- symmetrize_ltsa_matrix(B)
  n <- ncol(B)
  eig_k <- ltsa_validate_candidate_eig_k(eig_k, n)

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
    return(ltsa_as_candidate_result(
      dense_ltsa_eig(B, eig_k),
      B = B,
      eig_k = eig_k,
      backend = "dense_eigen",
      convergence_known = TRUE
    ))
  }

  lambda_probe <- NULL
  if (is.null(lambda_max)) {
    tsmessage("Finding largest eigenvalue")
    lambda_probe <- ltsa_irlba_lambda_max_probe(B)
    lambda_max <- lambda_probe$value
  } else {
    lambda_max <- ltsa_validate_backend_lambda_max(
      lambda_max,
      B,
      backend = "irlba"
    )
  }
  shift <- lambda_max + ltsa_shift_margin(lambda_max, shift_eps)
  X_shift <- ltsa_shift_for_smallest(B, shift)

  tsmessage("Decomposing shifted matrix")
  args <- lmerge(
    list(
      A = X_shift,
      nv = eig_k,
      nu = 0L
    ),
    varargs
  )
  res <- do.call(irlba::irlba, args)

  ltsa_shifted_candidate_result(
    B = B,
    vectors = res$v,
    eig_k = eig_k,
    backend = "irlba",
    lambda_max = lambda_max,
    lambda_probe = lambda_probe,
    shift = shift,
    shift_eps = shift_eps,
    shifted_values = res$d,
    niter = res$iter %||% NA_integer_,
    mprod = res$mprod %||% NA_integer_,
    opts = args,
    returned_columns = ncol(res$v)
  )
}

# svdr candidate provider. Like irlba, this produces a candidate subspace only;
# the common Rayleigh-Ritz selector decides which nonconstant directions to
# return.
ltsa_svdr_candidate_provider <- function(
  B,
  eig_k,
  ...,
  lambda_max = NULL,
  verbose = FALSE,
  shift_eps = 1e-6,
  dense_n = 100L,
  dense_fraction = 0.5
) {
  varargs <- list(...)
  B <- symmetrize_ltsa_matrix(B)
  n <- ncol(B)
  eig_k <- ltsa_validate_candidate_eig_k(eig_k, n)

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
    return(ltsa_as_candidate_result(
      dense_ltsa_eig(B, eig_k),
      B = B,
      eig_k = eig_k,
      backend = "dense_eigen",
      convergence_known = TRUE
    ))
  }

  lambda_probe <- NULL
  if (is.null(lambda_max)) {
    tsmessage("Finding largest eigenvalue")
    lambda_probe <- ltsa_svdr_lambda_max_probe(B)
    lambda_max <- lambda_probe$value
  } else {
    lambda_max <- ltsa_validate_backend_lambda_max(
      lambda_max,
      B,
      backend = "svdr"
    )
  }
  shift <- lambda_max + ltsa_shift_margin(lambda_max, shift_eps)
  X_shift <- ltsa_shift_for_smallest(B, shift)

  tsmessage("Decomposing shifted matrix")
  args <- lmerge(
    list(
      x = X_shift,
      k = eig_k
    ),
    varargs
  )
  res <- do.call(irlba::svdr, args)

  ltsa_shifted_candidate_result(
    B = B,
    vectors = res$v,
    eig_k = eig_k,
    backend = "svdr",
    lambda_max = lambda_max,
    lambda_probe = lambda_probe,
    shift = shift,
    shift_eps = shift_eps,
    shifted_values = res$d,
    mprod = res$mprod %||% NA_integer_,
    opts = args,
    returned_columns = ncol(res$v)
  )
}

# Dense diagnostic path used for small matrices or when a requested candidate
# count is large relative to the matrix size.
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
