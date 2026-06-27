# Matrix, scaling, residual, near-zero, and lambda-probe helpers shared by all
# LTSA iterative eigenanalysis backends.

ltsa_default_eig_k <- function(ndim, n) {
  ndim <- check_whole_number(ndim, "ndim", min = 1)
  n <- check_whole_number(n, "n", min = 2)

  min(n - 1L, max(12L, ndim + 2L))
}

ltsa_validate_eig_k <- function(eig_k, ndim, n) {
  if (is.null(eig_k)) {
    return(ltsa_default_eig_k(ndim = ndim, n = n))
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
    stop("eig_k must satisfy ndim + 1 <= eig_k < n", call. = FALSE)
  }
  as.integer(eig_k)
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

ltsa_rspectra_opts <- function(eig_k = NULL, n = NULL) {
  opts <- list(tol = 1e-6)
  if (!is.null(eig_k) && !is.null(n)) {
    opts$ncv <- ltsa_default_ncv(n, eig_k)
  }
  opts
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

ltsa_validate_lambda_max <- function(lambda_max, B) {
  ltsa_validate_backend_lambda_max(lambda_max, B, backend = "RSpectra")
}

ltsa_validate_backend_lambda_max <- function(lambda_max, B, backend) {
  if (length(lambda_max) < 1L || !is.finite(lambda_max[[1L]])) {
    stop(
      backend,
      " largest-eigenvalue probe returned a non-finite value",
      call. = FALSE
    )
  }

  lambda_max <- max(0, as.numeric(lambda_max[[1L]]))
  if (lambda_max <= 0 && !ltsa_matrix_is_effectively_zero(B)) {
    stop(
      backend,
      " largest-eigenvalue probe returned a non-positive value for ",
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
  opts <- ltsa_rspectra_opts(eig_k = 1L, n = nrow(B))
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
