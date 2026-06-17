# Low-level eigensolver calls retained for compatibility with existing internal
# LTSA tests and helper wrappers.

# RSpectra backend call. This avoids shift-invert near zero by estimating the
# largest eigenvalue, forming shift * I - B, and solving for largest algebraic
# eigenvectors of the shifted problem. nconv is treated as a hard convergence
# gate before shared LTSA postprocessing sees the candidate block.
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
    ltsa_irlba_candidate_provider(
      B = X,
      eig_k = k,
      ...,
      lambda_max = lambda_max,
      verbose = verbose
    )$vectors
  }

svdr_eig <- function(
  X,
  k = ncol(X) - 1,
  ...,
  lambda_max = NULL,
  verbose = FALSE
) {
  ltsa_svdr_candidate_provider(
    B = X,
    eig_k = k,
    ...,
    lambda_max = lambda_max,
    verbose = verbose
  )$vectors
}
