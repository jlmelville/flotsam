ltsa_iterative_search_k <- function(ndim, n) {
  base_k <- ndim + 1L
  search_k <- max(ndim + 3L, 2L * base_k)
  if (search_k >= 0.5 * n) {
    return(base_k)
  }
  min(n - 1L, search_k)
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

  rayleigh <- colSums(vectors * (B %*% vectors)) / colSums(vectors * vectors)
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

rs_opt <- function() {
  list(
    tol = 1e-6
  )
}

# Avoid shift-invert here; near-zero LTSA/Laplacian eigenvalues can make sparse
# factorizations hang or skip eigenvectors.
# https://github.com/yixuan/spectra/issues/126
rs_eig <-
  function(X, k = ncol(X) - 1, ..., lambda_max = NULL, verbose = FALSE) {
    varargs <- list(...)
    if (is.null(lambda_max)) {
      tsmessage("Finding largest eigenvalue")
      args1 <- list(
        A = X,
        k = 1,
        opts = list(retvec = FALSE)
      )
      lambda_max <- do.call(RSpectra::eigs_sym, args1)$values
    }
    lm2 <- 2.0 * lambda_max
    X_shift <- shift_lap(X, lm2)

    tsmessage("Decomposing shifted matrix")
    args <-
      list(
        A = X_shift,
        k = k,
        which = "LM",
        opts = rs_opt()
      )
    args$opts <- lmerge(args$opts, varargs)
    do.call(RSpectra::eigs_sym, args)
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
