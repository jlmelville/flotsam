#' Local Tangent Space Alignment
#'
#' Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha, 2004)
#' for dimensionality reduction.
#'
#' @param X The input data matrix or data frame with one observation per row. If
#'   a data frame is supplied, non-numeric columns are ignored. At least one
#'   numeric column is required.
#' @param n_neighbors  The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation.
#' @param ndim The dimension of the space to embed into.
#' @param nn_method Method for finding nearest neighbors. Can be one of:
#'   * `"nnd"` Approximate nearest neighbors by Nearest Neighbor Descent.
#'   * `"exact"` Exact nearest neighbors by exhaustively comparing all items.
#'   Slow for large datasets.
#' @param eig_method How to carry out the eigendecomposition. Possible values are:
#'    * `"rspectra"` Use [RSpectra::eigs_sym()].
#'    * `"irlba"` Use [irlba::irlba()].
#'    * `"svdr"` Use [irlba::svdr()].
#'    * `"fullsvd"` Use the [base::svd()] function. This is only feasible for small
#'    datasets and should be used for diagnostic purposes only.
#'    * `"eig"` or `"eigen"` Use the [base::eigen()] function. This is only
#'    feasible for small datasets and should be used for diagnostic purposes
#'    only.
#' @param include_self Should an item be part of its own neighborhood? This has
#'   a minor effect on most results, but work by Zhang and co-workers (2017)
#'   suggests that this is in effect the main difference between LTSA and the
#'   Hessian Locally Linear Embedding (HLLE) method, so setting this to `FALSE`
#'   may allow emulating the HLLE method.
#' @param normalize If `TRUE` calculate the eigendecomposition on a normalized
#'   version of the Laplacian. This may be slightly easier to converge while
#'   giving similar results to the un-normalized case. It may also have suitable
#'   better properties if clustering is to be carried out on the eigenvectors.
#' @param ret_B If `TRUE`, return the matrix instead of the eigenvectors. This
#'   is mainly useful for diagnostic purposes if eigendecomposition is failing.
#' @param n_threads Number of threads to use. Applies only to the nearest
#'   neighbor calculation.
#' @param verbose If `TRUE` log information about progress to the console.
#' @param ... Extra arguments to be passed to the eigendecomposition method
#' specified by `eig_method`. For `"rspectra"`, arguments are passed to the
#' `opts` list. Suitable parameters include:
#'
#'   * `ncv` Number of Lanzcos vectors to use.
#'   * `tol` Tolerance.
#'   * `maxitr` Maximum number of iterations.
#'
#' For `"irlba"` suitable arguments are:
#'
#'   * `work` Working subspace dimension size.
#'   * `tol` Tolerance.
#'   * `maxit` Maximum number of iterations.
#'
#' For `"svdr"` suitable arguments are:
#'
#'   * `extra` Number of extra vectors to use.
#'   * `tol` Tolerance.
#'   * `it`  Maximum number of iterations.
#'
#' For more details see the documentation for [RSpectra::eigs_sym()],
#' [irlba::irlba()] and [irlba::svdr()] functions, respectively. Don't pass
#' other arguments unless you know what you are doing, as it may cause the
#' `ltsa` to fail.
#'
#' @references
#' Zhang, Z., & Zha, H. (2004).
#' Principal manifolds and nonlinear dimensionality reduction via tangent space alignment.
#' *SIAM journal on scientific computing*, *26*(1), 313-338.
#' <https://doi.org/10.1137/S1064827502419154>
#'
#' Zhang, S., Ma, Z., & Tan, H. (2017).
#' On the Equivalence of HLLE and LTSA.
#' *IEEE transactions on cybernetics*, *48*(2), 742-753.
#' <https://doi.org/10.1109/TCYB.2017.2655338>
#'
#' @examples
# Create a "swiss roll": 2D rectangle rolled up in 3D
#' n <- 1000
#' max_z <- 10
#'
#' phi <- stats::runif(n, min = 1.5 * pi, max = 4.5 * pi)
#' x <- phi * cos(phi)
#' y <- phi * sin(phi)
#' z <- stats::runif(n, max = max_z)
#' swiss_roll <- data.frame(x, y, z)
#'
#' # unroll it
#' swiss_ltsa <- ltsa(swiss_roll)
#' plot(swiss_ltsa, col = phi)
#'
#' # compare with PCA
#' swiss_pca <- stats::prcomp(swiss_roll, rank. = 2, scale = FALSE, retx = TRUE)$x
#' plot(swiss_pca, col = phi)
#' @export
ltsa <-
  function(X,
           n_neighbors = 15,
           ndim = 2,
           nn_method = "nnd",
           eig_method = "rspectra",
           include_self = TRUE,
           normalize = FALSE,
           ret_B = FALSE,
           n_threads = 0,
           verbose = FALSE,
           ...) {
    X <- x2m(X)
    args <- validate_ltsa_args(
      X = X,
      n_neighbors = n_neighbors,
      ndim = ndim,
      nn_method = nn_method,
      eig_method = eig_method,
      include_self = include_self,
      normalize = normalize,
      ret_B = ret_B,
      n_threads = n_threads,
      verbose = verbose
    )
    X <- args$X
    k <- args$n_neighbors
    ndim <- args$ndim
    nn_method <- args$nn_method
    eig_method <- args$eig_method
    include_self <- args$include_self
    normalize <- args$normalize
    ret_B <- args$ret_B
    n_threads <- args$n_threads
    verbose <- args$verbose

    nn_fun <- switch(nn_method,
      exact = rnndescent::brute_force_knn,
      nnd = rnndescent::nnd_knn
    )
    tsmessage("Finding nearest neighbors with method '", nn_method, "'")
    nn_args <- list(
      data = X,
      k = ifelse(include_self, k, k + 1),
      n_threads = n_threads,
      verbose = FALSE
    )
    nn <- do.call(nn_fun, nn_args)
    nn$dist <- NULL
    mode(nn$idx) <- "integer"

    n <- nrow(X)

    tsmessage("Getting neighborhoods")
    assembly <- assemble_ltsa_B(
      X = X,
      nn_idx = nn$idx,
      ndim = ndim,
      include_self = include_self,
      verbose = verbose
    )
    B <- assembly$B
    rank_deficient_count <- assembly$rank_deficient_count
    min_local_rank <- assembly$min_local_rank

    if (rank_deficient_count > 0) {
      warning(
        rank_deficient_count,
        " local neighborhoods had numerical rank below ndim = ",
        ndim,
        "; lower-dimensional local bases were used. Minimum local rank was ",
        min_local_rank,
        ".",
        call. = FALSE
      )
    }

    if (ret_B) {
      return(B)
    }

    Dinvs <- NULL
    if (normalize) {
      tsmessage("Forming shifted Lsym")
      sres <- norm_and_shift_L(B)
      Dinvs <- sres$Dinvs
      B <- sres$Lshift
    }

    tsmessage("Performing eigenanalysis")

    out <- tryCatch(
      {
        if (normalize) {
          switch(eig_method,
            irlba = {
              eig_args <- lmerge(
                list(
                  A = B,
                  nv = ndim + 1,
                  nu = 0
                ),
                list(...)
              )
              tsmessage("Calling irlba")
              res <-
                do.call(irlba::irlba, eig_args)$v[, 2:(ndim + 1)]
            },
            svdr = {
              eig_args <- lmerge(
                list(
                  x = B,
                  k = ndim + 1
                ),
                list(...)
              )
              tsmessage("Calling irlba svdr")
              res <-
                do.call(irlba::svdr, eig_args)$v[, 2:(ndim + 1)]
            },
            rspectra = {
              tsmessage("Calling rspectra")

              eig_args <- list(
                A = B,
                k = ndim + 1,
                which = "LM",
                opts = rs_opt()
              )
              eig_args$opts <- lmerge(eig_args$opts, list(...))

              res <-
                do.call(RSpectra::eigs_sym, eig_args)$vectors
              if (ncol(res) != ndim + 1) {
                stop("Can't find enough vectors")
              }
              res <- res[, 2:(ndim + 1)]
            },
            fullsvd = {
              tsmessage("Using full SVD")
              res <-
                svd(as.matrix(B), nv = ndim + 1, nu = 0)$v[, 2:(ndim + 1)]
            },
            eig = {
              tsmessage("Using full eigenvalue decomposition")
              res <- eigen(as.matrix(B))$vectors[, 2:(ndim + 1)]
            }
          )
          Dinvs * res
        } else {
          switch(eig_method,
            rspectra = {
              tsmessage("Calling rspectra")
              k_search <- ltsa_iterative_search_k(ndim, ncol(B))
              res <-
                rs_eig(B, k = k_search, ..., verbose = verbose)$vectors
              if (ncol(res) < ndim) {
                stop("Can't find enough vectors")
              }
              select_ltsa_embedding_vectors(B, res, ndim)
            },
            irlba = {
              k_search <- ltsa_iterative_search_k(ndim, ncol(B))
              res <- irlba_eig(B, k = k_search, ...)
              select_ltsa_embedding_vectors(B, res, ndim)
            },
            svdr = {
              k_search <- ltsa_iterative_search_k(ndim, ncol(B))
              res <- svdr_eig(B, k = k_search, ...)
              select_ltsa_embedding_vectors(B, res, ndim)
            },
            fullsvd = {
              tsmessage("Using full SVD")
              res <- svd(as.matrix(B))$v
              nvec <- ncol(res)
              res <- res[, rev((nvec - ndim):(nvec - 1))]
            },
            eig = {
              tsmessage("Using full eigenvalue decomposition")
              res <- eigen(as.matrix(B), symmetric = TRUE)$vectors
              nvec <- ncol(res)
              res <- res[, rev((nvec - ndim):(nvec - 1))]
            }
          )
        }
      },
      error = function(e) {
        stop("Eigenanalysis failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    tsmessage("Finished")
    out
  }

assemble_ltsa_B <- function(X,
                            nn_idx,
                            ndim,
                            include_self,
                            verbose = FALSE,
                            preserve_pattern = FALSE) {
  if (preserve_pattern) {
    return(assemble_ltsa_B_compact(
      X = X,
      nn_idx = nn_idx,
      ndim = ndim,
      include_self = include_self,
      verbose = verbose,
      preserve_pattern = TRUE
    ))
  }

  assemble_ltsa_B_append(
    X = X,
    nn_idx = nn_idx,
    ndim = ndim,
    include_self = include_self,
    verbose = verbose
  )
}

assemble_ltsa_B_append <- function(X,
                                   nn_idx,
                                   ndim,
                                   include_self,
                                   verbose = FALSE) {
  n <- nrow(X)
  weight_idx <- ltsa_weight_neighborhoods(nn_idx, include_self)
  k <- ncol(weight_idx)
  builder <- ltsa_triplet_builder_create(t(weight_idx), k)

  prev_time <- Sys.time()
  rank_deficient_count <- 0L
  min_local_rank <- ndim
  for (i in seq_len(n)) {
    if (verbose) {
      curr_time <- Sys.time()
      if (difftime(curr_time, prev_time, units = "secs") > 60) {
        tsmessage("processed ", i, " / ", n)
        prev_time <- curr_time
      }
    }

    nni <- weight_idx[i, ]
    local <- ltsa_local_weights(X, nni, ndim)
    if (local$rank < ndim) {
      rank_deficient_count <- rank_deficient_count + 1L
      min_local_rank <- min(min_local_rank, local$rank)
    }
    ltsa_triplet_builder_append(builder, nni, as.vector(local$Wi))
  }

  tsmessage("Assembling sparse matrix")
  components <- ltsa_triplet_builder_finalize(builder)
  B <- ltsa_components_to_dgCMatrix(components, n)

  list(
    B = B,
    rank_deficient_count = rank_deficient_count,
    min_local_rank = min_local_rank
  )
}

assemble_ltsa_B_compact <- function(X,
                                    nn_idx,
                                    ndim,
                                    include_self,
                                    verbose = FALSE,
                                    preserve_pattern = FALSE) {
  n <- nrow(X)
  weight_idx <- ltsa_weight_neighborhoods(nn_idx, include_self)
  k <- ncol(weight_idx)
  weights <- matrix(0, nrow = k * k, ncol = n)

  prev_time <- Sys.time()
  rank_deficient_count <- 0L
  min_local_rank <- ndim
  for (i in seq_len(n)) {
    if (verbose) {
      curr_time <- Sys.time()
      if (difftime(curr_time, prev_time, units = "secs") > 60) {
        tsmessage("processed ", i, " / ", n)
        prev_time <- curr_time
      }
    }

    local <- ltsa_local_weights(X, weight_idx[i, ], ndim)
    if (local$rank < ndim) {
      rank_deficient_count <- rank_deficient_count + 1L
      min_local_rank <- min(min_local_rank, local$rank)
    }
    weights[, i] <- as.vector(local$Wi)
  }

  tsmessage("Assembling sparse matrix")
  components <- ltsa_triplet_assembly_components(
    t(nn_idx),
    ncol(nn_idx),
    t(weight_idx),
    weights,
    k,
    preserve_pattern
  )

  B <- ltsa_components_to_dgCMatrix(components, n)

  list(
    B = B,
    rank_deficient_count = rank_deficient_count,
    min_local_rank = min_local_rank
  )
}

ltsa_components_to_dgCMatrix <- function(components, n) {
  methods::new(
    "dgCMatrix",
    i = components$i,
    p = components$p,
    x = components$x,
    Dim = as.integer(c(n, n))
  )
}

ltsa_weight_neighborhoods <- function(nn_idx, include_self) {
  if (include_self) {
    nn_idx
  } else {
    nn_idx[, -1L, drop = FALSE]
  }
}

ltsa_local_weights <- function(X, nni, ndim) {
  Xi <- (scale(X[nni, , drop = FALSE], center = TRUE, scale = FALSE))

  local_basis <- svecs(Xi, ndim)
  Vi <- local_basis$vectors
  k <- length(nni)
  Gi <- cbind(1 / sqrt(k), Vi)
  Wi <- -tcrossprod(Gi)
  diag(Wi) <- diag(Wi) + 1.0

  list(Wi = Wi, rank = local_basis$rank)
}

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

# get indices of the diagonal from the sparse matrix m
diag_spm <- function(m) {
  is <- m@i
  ps <- m@p
  n <- nrow(m)
  sapply(1:n, function(i, is, ps) {
    begin <- ps[i] + 1
    end <- ps[i + 1]
    begin + which(is[begin:end] == i - 1)
  }, is, ps) - 1
}

lsym_norm <- function(M, D) {
  M@x <- spm_times_scalar(M@p, M@x, D)
  D * M
}

# vI - m
shift_lap <- function(m, v = 2.0) {
  x <- m@x

  x <- -x
  ds <- diag_spm(m)
  x[ds] <- x[ds] + v

  m@x <- x
  m
}

norm_and_shift_L <- function(L) {
  Dinvs <- sqrt(1 / diag(L))
  list(Lshift = shift_lap(lsym_norm(L, Dinvs), 2.0), Dinvs = Dinvs)
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
  function(X,
           k = ncol(X) - 1,
           ...,
           lambda_max = NULL,
           verbose = FALSE) {
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
  function(X,
           k = ncol(X) - 1,
           ...,
           lambda_max = NULL,
           verbose = FALSE) {
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
    args <- lmerge(list(
      A = X_shift,
      nv = k,
      nu = 0
    ), varargs)
    res <- do.call(irlba::irlba, args)
    res$v
  }

svdr_eig <- function(X,
                     k = ncol(X) - 1,
                     ...,
                     lambda_max = NULL,
                     verbose = FALSE) {
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
  args <- lmerge(list(
    x = X_shift,
    k = k
  ), varargs)
  res <- do.call(irlba::svdr, args)
  res$v
}

svecs <- function(X, ndim = 2) {
  max_rank <- min(dim(X))
  ndim <- min(ndim, max_rank)

  if (ndim < 0.5 * max_rank) {
    res <- tryCatch(
      suppressWarnings(irlba::irlba(X, nu = ndim, nv = 0)),
      error = function(e) NULL
    )
    if (is.null(res)) {
      res <- svd(X, nu = ndim, nv = 0)
    }
  } else {
    res <- svd(X, nu = ndim, nv = 0)
  }
  rank <- numerical_rank(res$d, dim(X))
  keep <- seq_len(min(ndim, length(res$d), ncol(res$u)))
  keep <- keep[res$d[keep] > rank$tol]

  list(
    vectors = res$u[, keep, drop = FALSE],
    rank = rank$rank,
    tol = rank$tol
  )
}

numerical_rank <- function(d, dims) {
  if (length(d) == 0 || max(d) == 0) {
    return(list(rank = 0L, tol = 0))
  }
  # standard conservative tolerance for numerical rank
  tol <- max(dims) * max(d) * .Machine$double.eps
  list(rank = sum(d > tol), tol = tol)
}
