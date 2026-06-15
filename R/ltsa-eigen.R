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
  residual_scale <- max(lambda_max, 1)

  list(
    absolute_residuals = absolute_residual,
    scaled_residuals = absolute_residual / residual_scale,
    residual_scale = residual_scale
  )
}

ltsa_validate_lambda_max <- function(lambda_max, B) {
  if (length(lambda_max) < 1L || !is.finite(lambda_max[[1L]])) {
    stop("RSpectra largest-eigenvalue probe returned a non-finite value", call. = FALSE)
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

# RSpectra sometimes fails to return the trivial constant vector, so we can't
# blindly drop the first vector. We also want to make sure we have the
# eigenvectors in the expected order. We fetch an "over-complete" set of
# vectors and then sort them by their Rayleigh quotient on the original matrix.
# A small Rayleigh quotient means the vector is close to a small-eigenvalue
# direction of B.
select_ltsa_embedding_vectors <- function(B, vectors, ndim) {
  if (ncol(vectors) < ndim) {
    stop("Can't find enough vectors", call. = FALSE)
  }

  rayleigh <- ltsa_rayleigh_values(B, vectors)
  vectors <- vectors[, order(rayleigh), drop = FALSE]

  first <- vectors[, 1L]
  centered_norm <- sqrt(sum((first - mean(first))^2))
  first_norm <- sqrt(sum(first^2))
  drop_trivial <- first_norm > 0 && centered_norm <= 1e-4 * first_norm
  start <- if (drop_trivial && ncol(vectors) > ndim) 2L else 1L
  end <- start + ndim - 1L
  if (end > ncol(vectors)) {
    stop("Can't find enough vectors", call. = FALSE)
  }

  vectors[, start:end, drop = FALSE]
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
  function(X,
           k = ncol(X) - 1,
           ...,
           lambda_max = NULL,
           verbose = FALSE,
           shift_eps = 1e-6,
           dense_n = 100L,
           dense_fraction = 0.5) {
    varargs <- list(...)
    B <- symmetrize_ltsa_matrix(X)
    n <- ncol(B)
    eig_k <- as.integer(k)
    if (length(eig_k) != 1L || is.na(eig_k) || eig_k < 1L || eig_k >= n) {
      stop("eig_k must be a positive integer less than the matrix dimension", call. = FALSE)
    }

    if (
      length(varargs) == 0L &&
        ltsa_use_dense_eig(n, eig_k, dense_n = dense_n, dense_fraction = dense_fraction)
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
      stop("RSpectra returned fewer LTSA candidate vectors than requested", call. = FALSE)
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
