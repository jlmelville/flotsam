# Eigensolver backends return candidate vectors for the final Rayleigh-Ritz
# selection step, plus compact backend metadata for diagnostics.
new_ltsa_candidates <- function(
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

# RSpectra uses a shifted largest-algebraic formulation and treats nconv as a
# hard convergence gate before shared LTSA postprocessing sees the candidates.
ltsa_rspectra_candidate_provider <- function(
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
  n <- ncol(B)

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
    lambda_max <- ltsa_validate_backend_lambda_max(
      lambda_max,
      B,
      backend = "RSpectra"
    )
  }
  shift <- lambda_max + ltsa_shift_margin(lambda_max, shift_eps)
  X_shift <- ltsa_shift_for_smallest(B, shift)

  tsmessage("Decomposing shifted matrix")
  opts <- ltsa_rspectra_opts(eig_k = eig_k, n = n)
  opts <- lmerge(opts, varargs)
  args <- list(
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
  finished <- ltsa_sort_and_score_candidates(
    B = B,
    vectors = res$vectors,
    eig_k = eig_k,
    shifted_values = res$values,
    lambda_max = lambda_max,
    backend = "RSpectra"
  )
  tsmessage(
    "RSpectra converged ",
    ifelse(is.na(nconv), eig_k, nconv),
    " / ",
    eig_k,
    " LTSA candidate vectors; max scaled residual = ",
    signif(max(finished$residuals$scaled_residuals), 4)
  )

  candidate <- new_ltsa_candidates(
    vectors = finished$vectors,
    values = finished$values,
    shifted_values = finished$shifted_values,
    backend = "rspectra",
    eig_k = eig_k,
    matrix = B,
    lambda_max = lambda_max,
    lambda_probe = lambda_probe,
    nconv = nconv,
    niter = res$niter %||% NA_integer_,
    nops = res$nops %||% NA_integer_,
    opts = opts,
    convergence_known = TRUE,
    returned_columns = ncol(res$vectors),
    converged_columns = ifelse(is.na(nconv), eig_k, nconv)
  )
  candidate$absolute_residuals <- finished$residuals$absolute_residuals
  candidate$scaled_residuals <- finished$residuals$scaled_residuals
  candidate$residual_scale <- finished$residuals$residual_scale
  candidate$shift <- shift
  candidate$shift_eps <- shift_eps
  candidate$shift_policy <- "lambda_max_plus_margin"
  candidate$solve_which <- "LA"
  candidate
}

ltsa_irlba_lambda_max_probe <- function(B) {
  args <- list(
    A = B,
    nv = 1L,
    nu = 0L
  )
  probe <- do.call(irlba::irlba, args)
  lambda_max <- ltsa_validate_backend_lambda_max(
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

ltsa_sort_and_score_candidates <- function(
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
  values <- ltsa_rayleigh_values(B, vectors)
  ord <- order(values)
  values <- values[ord]
  vectors <- vectors[, ord, drop = FALSE]
  if (!is.null(shifted_values) && length(shifted_values) >= eig_k) {
    shifted_values <- shifted_values[seq_len(eig_k)][ord]
  }
  residuals <- ltsa_ritz_residuals(B, vectors, values, lambda_max)

  list(
    vectors = vectors,
    values = values,
    shifted_values = shifted_values,
    residuals = residuals
  )
}

# Shifted solves are re-valued against the original LTSA matrix. The shifted
# backend values are retained only as metadata; final selection uses B.
new_ltsa_shifted_candidates <- function(
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
  finished <- ltsa_sort_and_score_candidates(
    B = B,
    vectors = vectors,
    eig_k = eig_k,
    shifted_values = shifted_values,
    lambda_max = lambda_max,
    backend = backend
  )
  tsmessage(
    backend,
    " returned ",
    eig_k,
    " LTSA candidate vectors; max scaled residual = ",
    signif(max(finished$residuals$scaled_residuals), 4)
  )

  candidate <- new_ltsa_candidates(
    vectors = finished$vectors,
    values = finished$values,
    shifted_values = finished$shifted_values,
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
  candidate$absolute_residuals <- finished$residuals$absolute_residuals
  candidate$scaled_residuals <- finished$residuals$scaled_residuals
  candidate$residual_scale <- finished$residuals$residual_scale
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
  n <- ncol(B)

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

  new_ltsa_shifted_candidates(
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
  n <- ncol(B)

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

  new_ltsa_shifted_candidates(
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

  candidate <- new_ltsa_candidates(
    vectors = vectors,
    values = values,
    backend = backend,
    eig_k = eig_k,
    matrix = B,
    lambda_max = lambda_max,
    nconv = eig_k,
    convergence_known = TRUE,
    returned_columns = ncol(vectors),
    converged_columns = eig_k
  )
  lmerge(
    candidate,
    list(
      absolute_residuals = residuals$absolute_residuals,
      scaled_residuals = residuals$scaled_residuals,
      residual_scale = residuals$residual_scale,
      shift = NA_real_,
      shift_eps = NA_real_,
      shift_policy = "dense",
      solve_which = NA_character_
    )
  )
}

ltsa_eig_candidate_provider <- function(
  B,
  eig_k,
  ...,
  lambda_max = NULL,
  verbose = FALSE
) {
  dense_ltsa_eig(B, eig_k, backend = "eig")
}
