# Matrix, scaling, residual, near-zero, and lambda-probe helpers shared by all
# LTSA iterative eigenanalysis backends.

default_eig_k <- function(ndim, n) {
  ndim <- check_whole_number(ndim, "ndim", min = 1)
  n <- check_whole_number(n, "n", min = 2)

  min(n - 1L, max(12L, ndim + 2L))
}

validate_eig_k <- function(eig_k, ndim, n, dimension_name = "ndim") {
  if (is.null(eig_k)) {
    return(default_eig_k(ndim = ndim, n = n))
  }

  ndim <- check_whole_number(ndim, "ndim", min = 1)
  n <- check_whole_number(n, "n", min = 2)

  if (
    !is.numeric(eig_k) ||
      length(eig_k) != 1L ||
      is.na(eig_k) ||
      !is.finite(eig_k) ||
      eig_k != floor(eig_k) ||
      eig_k < ndim + 1L ||
      eig_k >= n
  ) {
    stop(
      "eig_k must satisfy ",
      dimension_name,
      " + 1 <= eig_k < n",
      call. = FALSE
    )
  }
  as.integer(eig_k)
}

validate_spectral_eig_k <- function(eig_k, ndim, spectral_dim, n) {
  dimension_name <- if (spectral_dim > ndim) {
    "spectral_dim"
  } else {
    "ndim"
  }
  eig_k <- validate_eig_k(
    eig_k = eig_k,
    ndim = spectral_dim,
    n = n,
    dimension_name = dimension_name
  )
  if (spectral_dim > ndim && eig_k < spectral_dim + 2L) {
    stop(
      "eig_k must satisfy spectral_dim + 2 <= eig_k < n when ",
      "spectral_dim exceeds ndim",
      call. = FALSE
    )
  }
  eig_k
}

symmetrize_operator <- function(B) {
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

matrix_max_abs <- function(B) {
  if (methods::is(B, "sparseMatrix")) {
    if (length(B@x) == 0L) {
      return(0)
    }
    return(max(abs(B@x)))
  }
  max(abs(B))
}

matrix_is_effectively_zero <- function(B) {
  matrix_max_abs(B) <= sqrt(.Machine$double.eps)
}

default_ncv <- function(n, eig_k) {
  min(n, max(20L, 4L * eig_k + 10L))
}

rspectra_options <- function(eig_k = NULL, n = NULL) {
  opts <- list(tol = 1e-6)
  if (!is.null(eig_k) && !is.null(n)) {
    opts$ncv <- default_ncv(n, eig_k)
  }
  opts
}

use_dense_eigensolver <- function(
  n,
  eig_k,
  dense_n = 100L,
  dense_fraction = 0.5
) {
  n <= dense_n || eig_k >= dense_fraction * n
}

shift_margin <- function(lambda_max, shift_eps = 1e-6) {
  scale <- max(1, abs(lambda_max))
  max(shift_eps * scale, 100 * .Machine$double.eps * scale)
}

shift_for_smallest <- function(B, shift) {
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

normalize_columns <- function(V) {
  V <- as.matrix(V)
  norms <- sqrt(colSums(V * V))
  sweep(V, 2L, ifelse(norms == 0, 1, norms), "/")
}

rayleigh_values <- function(B, vectors) {
  vectors <- as.matrix(vectors)
  BV <- as.matrix(B %*% vectors)
  denom <- colSums(vectors * vectors)
  colSums(vectors * BV) / denom
}

ritz_residuals <- function(B, vectors, values, lambda_max) {
  vectors <- normalize_columns(vectors)
  if (length(values) != ncol(vectors)) {
    stop(
      "internal error: residual value/vector length mismatch",
      call. = FALSE
    )
  }
  BV <- as.matrix(B %*% vectors)
  residual <- BV - sweep(vectors, 2L, values, "*")
  absolute_residual <- sqrt(colSums(residual * residual))
  scale <- residual_scale(lambda_max)

  list(
    absolute_residuals = absolute_residual,
    scaled_residuals = absolute_residual / scale,
    residual_scale = scale
  )
}

residual_scale <- function(lambda_max) {
  if (
    is.null(lambda_max) ||
      length(lambda_max) < 1L ||
      !is.finite(lambda_max[[1L]])
  ) {
    return(1)
  }

  max(as.numeric(lambda_max[[1L]]), 1)
}

default_null_vector <- function(n) {
  rep(1, n)
}

normalize_null_vector <- function(null_vector, n) {
  if (length(null_vector) != n || any(!is.finite(null_vector))) {
    stop(
      "LTSA null vector must be finite and match the matrix dimension",
      call. = FALSE
    )
  }

  null_norm <- sqrt(sum(null_vector * null_vector))
  if (!is.finite(null_norm) || null_norm <= 0) {
    stop("LTSA null vector must have positive norm", call. = FALSE)
  }

  null_vector / null_norm
}

near_zero_tolerance <- function(lambda_max, tol = sqrt(.Machine$double.eps)) {
  tol * residual_scale(lambda_max)
}

validate_backend_lambda_max <- function(lambda_max, B, backend) {
  if (length(lambda_max) < 1L || !is.finite(lambda_max[[1L]])) {
    stop(
      backend,
      " largest-eigenvalue probe returned a non-finite value",
      call. = FALSE
    )
  }

  lambda_max <- max(0, as.numeric(lambda_max[[1L]]))
  if (lambda_max <= 0 && !matrix_is_effectively_zero(B)) {
    stop(
      backend,
      " largest-eigenvalue probe returned a non-positive value for ",
      "a nonzero LTSA matrix",
      call. = FALSE
    )
  }
  lambda_max
}

validate_lambda_probe <- function(probe, B) {
  nconv <- probe$nconv
  if (!is.null(nconv) && length(nconv) > 0L && nconv < 1L) {
    stop("RSpectra largest-eigenvalue probe did not converge", call. = FALSE)
  }

  validate_backend_lambda_max(probe$values, B, backend = "RSpectra")
}

rspectra_lambda_max_probe <- function(B, varargs) {
  opts <- rspectra_options(eig_k = 1L, n = nrow(B))
  opts <- merge_named_lists(opts, varargs)
  opts$retvec <- FALSE

  args <- list(
    A = B,
    k = 1L,
    which = "LA",
    opts = opts
  )
  probe <- do.call(RSpectra::eigs_sym, args)
  lambda_max <- validate_lambda_probe(probe, B)

  list(
    value = lambda_max,
    nconv = probe$nconv %||% NA_integer_,
    niter = probe$niter %||% NA_integer_,
    nops = probe$nops %||% NA_integer_,
    which = "LA",
    opts = opts
  )
}
