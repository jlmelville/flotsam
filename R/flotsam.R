#' Local Tangent Space Alignment
#'
#' Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha, 2004)
#' for dimensionality reduction.
#'
#' @details
#' LTSA forms a sparse alignment matrix `B` and then returns the low-energy
#' nonconstant eigenvectors of that matrix. For iterative final eigenanalysis,
#' `flotsam` does more than return the raw columns from the eigensolver. It
#' asks the selected backend for a candidate block, projects out the known
#' constant null vector, performs a small Rayleigh-Ritz extraction on the
#' original LTSA matrix, sorts by the resulting Ritz values, and returns the
#' first `ndim` nonconstant Ritz vectors.
#'
#' This postprocessing is intended to make clustered low-energy spectra less
#' fragile. The iterative methods over-request candidate vectors by default,
#' check post-null rank, compute scaled residuals against `B`, inspect the
#' boundary gap after the selected block, and report near-zero nonconstant mode
#' counts when `verbose = TRUE`. A weak boundary gap does not make the result
#' fail, but it can trigger a small number of diagnostic expansion attempts and
#' an ambiguity warning because the requested `ndim` may cut through a larger
#' low-energy eigenspace.
#'
#' The `"rspectra"` path first estimates the largest eigenvalue and solves a
#' shifted largest-algebraic problem instead of using shift-invert near zero.
#' RSpectra's `nconv` metadata is treated as a hard backend convergence check.
#' The `"irlba"` and `"svdr"` paths use the same Rayleigh-Ritz selection and
#' diagnostics, but they rely on post-hoc residual checks because those
#' backends do not expose RSpectra-style convergence counts.
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
#' @param eig_method How to carry out the final eigendecomposition. Possible
#'   values are:
#'    * `"rspectra"` Use [RSpectra::eigs_sym()] to generate candidate vectors,
#'      then apply null-aware Rayleigh-Ritz postprocessing. This is the default
#'      and can report hard backend convergence failures through RSpectra's
#'      `nconv` metadata.
#'    * `"irlba"` Use [irlba::irlba()] to generate candidate vectors, then apply
#'      the same null-aware Rayleigh-Ritz postprocessing. This backend does not
#'      report RSpectra-style convergence counts, so convergence issues are
#'      detected with post-hoc residual diagnostics.
#'    * `"svdr"` Use [irlba::svdr()] to generate candidate vectors, then apply
#'      the same null-aware Rayleigh-Ritz postprocessing. This backend does not
#'      report RSpectra-style convergence counts, so convergence issues are
#'      detected with post-hoc residual diagnostics.
#'    * `"fullsvd"` Use the [base::svd()] function. This is only feasible for
#'      small datasets and should be used for diagnostic purposes only.
#'    * `"eig"` or `"eigen"` Use the [base::eigen()] function. This is only
#'      feasible for small datasets and should be used for diagnostic purposes
#'      only. Dense `"eig"` is the better diagnostic reference when algebraic
#'      eigenvalue ordering matters.
#' @param include_self Should an item be part of its own neighborhood? This has
#'   a minor effect on most results, but work by Zhang and co-workers (2017)
#'   suggests that this is in effect the main difference between LTSA and the
#'   Hessian Locally Linear Embedding (HLLE) method, so setting this to `FALSE`
#'   may allow emulating the HLLE method.
#' @param normalize If `TRUE`, calculate the final decomposition on a normalized
#'   version of the LTSA alignment matrix. For iterative methods (`"rspectra"`,
#'   `"irlba"`, and `"svdr"`), the normalized operator uses the same null-aware
#'   Rayleigh-Ritz selection as the unnormalized path, then transforms the
#'   selected vectors back to output coordinates. This may be slightly easier to
#'   converge while giving similar results to the unnormalized case. It may also
#'   have better properties if clustering is to be carried out on the
#'   eigenvectors.
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
#' `opts` list used by [RSpectra::eigs_sym()]. Suitable parameters include:
#'
#'   * `ncv` Number of Lanczos vectors to use.
#'   * `tol` Tolerance.
#'   * `maxitr` Maximum number of iterations.
#'
#' For `"irlba"`, arguments are passed to [irlba::irlba()]. Suitable arguments
#' include:
#'
#'   * `work` Working subspace dimension size.
#'   * `tol` Tolerance.
#'   * `maxit` Maximum number of iterations.
#'   * `reorth` Whether to use full reorthogonalization.
#'
#' For `"svdr"`, arguments are passed to [irlba::svdr()]. Suitable arguments
#' include:
#'
#'   * `extra` Number of extra vectors to use.
#'   * `tol` Tolerance.
#'   * `it`  Maximum number of iterations.
#'
#' The iterative methods all use shared Rayleigh-Ritz postprocessing, residual
#' diagnostics, and adaptive candidate expansion. However, only RSpectra reports
#' explicit convergence counts; `"irlba"` and `"svdr"` rely on post-hoc residual
#' diagnostics instead.
#'
#' Advanced LTSA eigenanalysis controls are also accepted here. `initial_extra`
#' and `max_extra` control how many candidate vectors beyond the constant null
#' vector, `ndim` embedding vectors, and one boundary-gap vector may be
#' requested. `gap_expansion_steps` controls how many extra attempts are made
#' when rank and residual diagnostics are good but the boundary gap is weak.
#' The default is one extra attempt, confirming a weak-gap result once before
#' returning it.
#' `resid_tol` and `gap_tol` set the scaled residual and global boundary-gap
#' diagnostic tolerances. `strict_rescue`, `strict_rescue_tol`,
#' `strict_rescue_maxitr`, `strict_rescue_maxit`, `strict_rescue_it`, and
#' `strict_rescue_extra` control the stricter rerun used when the selected
#' block appears to contain only part of a near-zero nonconstant eigenspace.
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
