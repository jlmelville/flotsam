#' Local Tangent Space Alignment
#'
#' Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha, 2004)
#' for dimensionality reduction.
#'
#' @param X The input data matrix or dataframe with one observation per row.
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
#'    * `"svdr` Use [irlba::svdr()].
#'    * `"fullsvd"` Use the [base::svd()] function. This is only feasible for small
#'    datasets and should be used for diagnostic purposes only.
#'    * `"eigen"` Use the [base::eigen()] function. This is only feasible for
#'    small datasets and should be used for diagnostic purposes only.
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
    k <- n_neighbors
    X <- x2m(X)

    nn_fun <- switch(nn_method,
      exact = rnndescent::brute_force_knn,
      nnd = rnndescent::nnd_knn,
      stop("Unknown nn method '", nn_method, "'")
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
    indim <- ncol(X)

    tsmessage("Allocating space for sparse matrix")
    B <- create_sparse(nn$idx, verbose = verbose)

    Bx <- B@x

    tsmessage("Getting neighborhoods")
    # This is the slow bit
    prev_time <- Sys.time()
    for (i in 1:n) {
      # N
      if (verbose) {
        curr_time <- Sys.time()
        if (difftime(curr_time, prev_time, units = "secs") > 60) {
          tsmessage("processed ", i, " / ", n)
          prev_time <- curr_time
        }
      }

      nni <- nn$idx[i, ]
      if (!include_self) {
        nni <- nni[-1]
      }
      # centered neighborhood k X M
      Xi <- (scale(X[nni, ], center = TRUE, scale = FALSE))

      # get the right singular vectors
      Vi <- svecs(Xi, ndim)

      # binding the 1/sqrt_k column handles the centering part during the crossprod
      Gi <- cbind(1 / sqrt(k), Vi)
      # Because Vs are orthonormal, Moore-Penrose inverse is V.V'
      # I - V.V'
      Wi <- -tcrossprod(Gi)
      diag(Wi) <- diag(Wi) + 1.0

      spi <- find_in_spsq(B, nni)
      Bx[spi] <- Bx[spi] + Wi
    }
    B@x <- Bx

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

    res <- NULL
    out <- tryCatch(
      {
        if (normalize) {
          switch(tolower(eig_method),
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
                opts = rs_opt(ncol(B)),
                sigma = rs_sigma_eps() + 2.0
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
            },
            stop(
              "Unknown method: '",
              eig_method,
              "' for normalized Laplacian"
            )
          )
          Dinvs * res
        } else {
          switch(tolower(eig_method),
            rspectra = {
              tsmessage("Calling rspectra")
              opt <- lmerge(rs_opt(ncol(B)), list(...))
              res <-
                rs_eig(B, k = ndim + 1, list(opt = opt), verbose = verbose)$vectors
              if (ncol(res) != ndim + 1) {
                stop("Can't find enough vectors")
              }
              res[, 2:(ndim + 1)]
            },
            irlba = {
              res <- irlba_eig(B, k = ndim + 1, ...)
              res[, 2:(ndim + 1)]
            },
            svdr = {
              res <- svdr_eig(B, k = ndim + 1, ...)
              res[, 2:(ndim + 1)]
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
            },
            stop(
              "Unknown method '",
              eig_method,
              "' with un-normalized Laplacian"
            )
          )
        }
      },
      error = function(cond) {
        tsmessage(cond)
        NULL
      }
    )
    if (is.null(out)) {
      tsmessage("Eigenanalysis failed")
    }
    tsmessage("Finished")
    out
  }

create_sparse <- function(nn_idx, verbose = FALSE) {
  n <- nrow(nn_idx)
  ij <- nbrhood_triplets(t(nn_idx), ncol(nn_idx))
  m <- methods::new(
    "dgTMatrix",
    i = ij$i,
    j = ij$j,
    x = rep(0, length(ij$i)),
    Dim = as.integer(c(n, n))
  )
  m <- m + Matrix::t(m)
  methods::as(m, "CsparseMatrix")
}

find_in_spsq <- function(m, is) {
  sparse_idxs(m@i, m@p, is)
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

rs_opt <- function(n) {
  list(
    tol = 1e-6,
    initvec = jitter(rep(1, n))
  )
}

rs_sigma_eps <- function() {
  0.0001
}

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
        sigma = rs_sigma_eps() + lm2,
        opts = rs_opt(ncol(X_shift))
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
    irlba::irlba(X, nu = ndim, nv = 0)$u
  } else {
    svd(X, nu = ndim, nv = 0)$u
  }
}
