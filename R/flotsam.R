#' Local Tangent Space Alignment
#'
#' Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha, 2004)
#' for dimensionality reduction.
#'
#' @param X The input data matrix or data frame with one observation per row. If
#'   a data frame is supplied, non-numeric columns are ignored. At least one
#'   numeric column is required.
#' @param n_neighbors  The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation. If `NULL`, the
#'   default is `15` when neighbors are computed, or inferred when a precomputed
#'   graph is supplied as `nn_method`.
#' @param ndim The dimension of the space to embed into.
#' @param nn_method Method for finding nearest neighbors, or a precomputed
#'   nearest-neighbor graph. Can be one of:
#'   * `"nnd"` Approximate nearest neighbors by Nearest Neighbor Descent.
#'   * `"exact"` Exact nearest neighbors by exhaustively comparing all items.
#'   Slow for large datasets.
#'   * A precomputed nearest-neighbor index matrix.
#'   * A nearest-neighbor result object with an `idx` matrix.
#'
#'   Precomputed index matrices should match the raw neighbor-search output used
#'   by `rnndescent`: one row per observation and 1-based indices into `X`. When
#'   `include_self = TRUE`, `ncol(nn_method)` must equal `n_neighbors` and each
#'   row must contain its own row index. When `include_self = FALSE`,
#'   `ncol(nn_method)` must equal `n_neighbors + 1` and the first column must
#'   contain the row index, because that column is dropped before LTSA assembly.
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
#' @param n_assembly_threads Number of threads to use when constructing the LTSA
#'   alignment matrix `B` after nearest neighbors are computed. The default
#'   `1` preserves the serial assembly path. Values greater than `1` opt into
#'   parallel construction of `B`, which can be faster but may increase peak
#'   memory use.
#' @param copy_max_mib Maximum size, in MiB, of the optional row-major dense
#'   copy of `X` used to speed up high-dimensional local Gram assembly. If this
#'   cap is exceeded, no row-major copy is made. The default is 256 MiB. Set to
#'   `0` to disable this copy. This code path is only used when the number of
#'   numeric columns in `X` is greater than the number of neighbors. Outside of
#'   synthetic datasets, you are quite likely to hit this code path.
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
  function(
    X,
    n_neighbors = NULL,
    ndim = 2,
    nn_method = "nnd",
    eig_method = "rspectra",
    include_self = TRUE,
    normalize = FALSE,
    ret_B = FALSE,
    n_threads = 0,
    n_assembly_threads = 1,
    copy_max_mib = 256,
    verbose = FALSE,
    ...
  ) {
    X <- x2m(X)
    neighbor_args <- normalize_ltsa_neighbor_args(
      nn_method = nn_method
    )
    nn_method <- neighbor_args$nn_method
    nn_idx <- neighbor_args$nn_idx

    args <- validate_ltsa_args(
      X = X,
      n_neighbors = n_neighbors,
      ndim = ndim,
      nn_method = nn_method,
      nn_idx = nn_idx,
      eig_method = eig_method,
      include_self = include_self,
      normalize = normalize,
      ret_B = ret_B,
      n_threads = n_threads,
      n_assembly_threads = n_assembly_threads,
      copy_max_mib = copy_max_mib,
      verbose = verbose
    )
    X <- args$X
    k <- args$n_neighbors
    nn_idx <- args$nn_idx
    ndim <- args$ndim
    nn_method <- args$nn_method
    eig_method <- args$eig_method
    include_self <- args$include_self
    normalize <- args$normalize
    ret_B <- args$ret_B
    n_threads <- args$n_threads
    n_assembly_threads <- args$n_assembly_threads
    copy_max_mib <- args$copy_max_mib
    verbose <- args$verbose

    neighbors <- prepare_ltsa_neighbors(
      X = X,
      n_neighbors = k,
      nn_method = nn_method,
      nn_idx = nn_idx,
      include_self = include_self,
      n_threads = n_threads,
      verbose = verbose
    )
    nn_idx <- neighbors$nn_idx

    tsmessage("Getting neighborhoods")
    assembly <- assemble_ltsa_B(
      X = X,
      nn_idx = nn_idx,
      ndim = ndim,
      include_self = include_self,
      n_assembly_threads = n_assembly_threads,
      copy_max_mib = copy_max_mib,
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
    normalized_nullvec <- NULL
    if (normalize) {
      if (eig_method %in% c("rspectra", "irlba", "svdr")) {
        tsmessage("Forming normalized Lsym")
        nres <- norm_lsym_L(B)
        Dinvs <- nres$Dinvs
        normalized_nullvec <- nres$nullvec
        B <- nres$Lsym
      } else {
        tsmessage("Forming shifted Lsym")
        sres <- norm_and_shift_L(B)
        Dinvs <- sres$Dinvs
        B <- sres$Lshift
      }
    }

    tsmessage("Performing eigenanalysis")

    out <- tryCatch(
      {
        if (normalize) {
          switch(
            eig_method,
            irlba = {
              tsmessage("Calling irlba")
              eig_res <- ltsa_irlba_ritz_eig(
                B,
                ndim = ndim,
                ...,
                nullvec = normalized_nullvec,
                verbose = verbose
              )
              res <- eig_res$vectors
            },
            svdr = {
              tsmessage("Calling irlba svdr")
              eig_res <- ltsa_svdr_ritz_eig(
                B,
                ndim = ndim,
                ...,
                nullvec = normalized_nullvec,
                verbose = verbose
              )
              res <- eig_res$vectors
            },
            rspectra = {
              tsmessage("Calling rspectra")
              eig_res <- ltsa_rspectra_ritz_eig(
                B,
                ndim = ndim,
                ...,
                nullvec = normalized_nullvec,
                verbose = verbose
              )
              res <- eig_res$vectors
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
          switch(
            eig_method,
            rspectra = {
              tsmessage("Calling rspectra")
              eig_res <- ltsa_rspectra_ritz_eig(
                B,
                ndim = ndim,
                ...,
                verbose = verbose
              )
              eig_res$vectors
            },
            irlba = {
              tsmessage("Calling irlba")
              eig_res <- ltsa_irlba_ritz_eig(
                B,
                ndim = ndim,
                ...,
                verbose = verbose
              )
              eig_res$vectors
            },
            svdr = {
              tsmessage("Calling irlba svdr")
              eig_res <- ltsa_svdr_ritz_eig(
                B,
                ndim = ndim,
                ...,
                verbose = verbose
              )
              eig_res$vectors
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
