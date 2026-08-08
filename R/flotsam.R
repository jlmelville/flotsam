#' Local Tangent Space Alignment
#'
#' Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha, 2004)
#' for dimensionality reduction.
#'
#' @details
#' `ltsa()` builds an LTSA alignment matrix from local neighborhoods and
#' returns the lowest nonconstant embedding vectors. The default
#' `output = "embedding"` returns only the embedding matrix.
#'
#' Use `output = "result"` to inspect compact diagnostics when convergence or
#' eigenspace ambiguity matters. The result contains the embedding,
#' eigenanalysis diagnostics, assembly diagnostics, and, when
#' `include_B = TRUE`, the assembled unnormalized matrix `B`. The diagnostic
#' `status` and `messages` fields summarize the requested solve; inspect those
#' condition-specific messages when the status is not `"ok"`.
#'
#' Use `output = "B"` to return the assembled unnormalized LTSA matrix and skip
#' final eigenanalysis.
#'
#' When `output = "result"`, `eigen$method` records the requested method and
#' `eigen$backend$name` records the backend actually used. Iterative methods
#' may use the dense fallback described under `Eigensolver controls`.
#'
#' @section Normalized LTSA:
#' `normalize = TRUE` applies symmetric Jacobi scaling to the assembled LTSA
#' alignment matrix `B`. Let `D = diag(B)`. The eigensolve is performed on
#' `D^(-1/2) B D^(-1/2)`, and the selected vectors are mapped back with
#' `D^(-1/2)`. Equivalently, this solves the generalized eigenproblem
#' `B v = lambda D v`.
#'
#' This keeps the same LTSA alignment energy but changes the pointwise mass
#' used by the final eigenproblem, so it can produce a different embedding from
#' the standard formulation. The diagonal `diag(B)` measures pointwise
#' participation in the LTSA alignment energy. Reverse-neighborhood
#' participation often contributes strongly, so the scaling can reduce the
#' influence of points that appear in many neighborhoods (i.e. hubness) and can
#' improve eigensolver convergence.
#'
#' However, other effects, such as boundary behavior, curvature, and local
#' scale, also contribute to `diag(B)` and it's hard to know a priori which
#' will dominate. Additionally, points with small diagonal mass are
#' relatively amplified by the `D^(-1/2)` scaling so the scaling may not always
#' be beneficial. Empirically, in testing across a variety of datasets, when
#' the assumption of a single-smooth manifold is violated, the normalized
#' embeddings often show more local structure and converge faster, compared to
#' the standard formulation.
#'
#' @section Precomputed neighbor input:
#' `nn_method` may be a precomputed 1-based neighbor index matrix or an object
#' with an `idx` matrix. With `include_self = TRUE`, each row must include its
#' own row index and have `n_neighbors` columns. With `include_self = FALSE`,
#' the first column must contain the row index and is dropped before LTSA
#' assembly, so the matrix must have `n_neighbors + 1` columns.
#'
#' @section Performance and memory controls:
#' `copy_max_mib` limits an optional row-major dense copy of `X` used during
#' high-dimensional local Gram assembly. Increase it only when that copy is
#' useful and enough memory is available; set it to `0` to disable the copy.
#'
#' @param X The input data matrix or data frame with one observation per row. If
#'   a data frame is supplied, non-numeric columns are ignored. At least one
#'   numeric column is required.
#' @param n_neighbors  The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation. If `NULL`, the
#'   default is `15` when neighbors are computed, or inferred when a precomputed
#'   graph is supplied as `nn_method`. It must be at least `ndim + 2` so each
#'   local projector has a residual direction beyond the constant and tangent
#'   subspaces.
#' @param ndim The dimension of the space to embed into.
#' @param nn_method Method for finding nearest neighbors, or a precomputed
#'   nearest-neighbor graph. Can be one of:
#'   * `"nnd"` Approximate nearest neighbors by Nearest Neighbor Descent.
#'   * `"exact"` Exact nearest neighbors by exhaustively comparing all items.
#'   Slow for large datasets.
#'   * A precomputed nearest-neighbor index matrix.
#'   * A nearest-neighbor result object with an `idx` matrix.
#' @param eig_method How to carry out the final eigendecomposition. Possible
#'   values are:
#'    * `"rspectra"` Use [RSpectra::eigs_sym()]. This is the default and can
#'      report hard backend convergence failures through RSpectra's `nconv`
#'      metadata.
#'    * `"irlba"` Use [irlba::irlba()] and post-hoc residual diagnostics.
#'    * `"svdr"` Use [irlba::svdr()] and post-hoc residual diagnostics.
#'    * `"eig"` or `"eigen"` Use the [base::eigen()] function. This is only
#'      feasible for small datasets and should be used for diagnostic purposes
#'      only. Dense `"eig"` is the better diagnostic reference when algebraic
#'      eigenvalue ordering matters.
#' @param eig_k Number of candidate vectors requested from the final
#'   eigensolver. If `NULL`, the default is
#'   `min(n - 1L, max(12L, ndim + 2L))`, where `n` is the number of
#'   observations. Must satisfy `ndim + 1 <= eig_k < n`. Larger values give
#'   the Rayleigh-Ritz postprocessing a wider candidate span.
#' @param output What to return:
#'   * `"embedding"` Return the embedding matrix. This is the default.
#'   * `"result"` Return a list containing the embedding, compact eigenanalysis
#'     diagnostics, assembly diagnostics, and optionally `B`.
#'   * `"B"` Return the assembled unnormalized LTSA matrix and skip final
#'     eigenanalysis.
#' @param include_B If `TRUE` and `output = "result"`, include the assembled
#'   unnormalized LTSA matrix `B` in the result object. Ignored for other output
#'   modes.
#' @param include_self Should an item be part of its own neighborhood? This has
#'   a minor effect on most results, but work by Zhang and co-workers (2017)
#'   suggests that this is in effect the main difference between LTSA and the
#'   Hessian Locally Linear Embedding (HLLE) method, so setting this to `FALSE`
#'   may allow emulating the HLLE method.
#' @param normalize If `TRUE`, use the normalized LTSA formulation described in
#'   the "Normalized LTSA" section. The default is `FALSE`, which uses the
#'   standard LTSA formulation.
#' @param n_threads Number of threads to use. Applies only to the nearest
#'   neighbor calculation.
#' @param n_assembly_threads Number of threads to use when constructing the LTSA
#'   alignment matrix `B` after nearest neighbors are computed. The default
#'   `1` preserves the serial assembly path. Values greater than `1` opt into
#'   parallel construction of `B`, which can be faster but may increase peak
#'   memory use.
#' @param copy_max_mib Maximum size, in MiB, of the optional row-major dense
#'   copy of `X` used during high-dimensional local Gram assembly. Set to `0`
#'   to disable this copy.
#' @param verbose If `TRUE` log information about progress to the console.
#' @param ... Optional eigensolver and diagnostic controls. See the
#'   `Eigensolver controls` section.
#'
#' @section Eigensolver controls:
#' Backend controls:
#'
#' * `"rspectra"`: `tol`, `maxitr`, and `ncv`. See
#'   \link[RSpectra:eigs_sym]{RSpectra::eigs_sym()} for their meanings. The
#'   package default for `tol` is `1e-6`.
#' * `"irlba"`: `tol`, `maxit`, and `reorth`.
#'   See \link[irlba:irlba]{irlba::irlba()} for their meanings.
#' * `"svdr"`: `tol`, `it`, and `extra`. See
#'   \link[irlba:svdr]{irlba::svdr()} for their meanings.
#'
#' `resid_tol` and `gap_tol` set diagnostic thresholds for candidate residuals
#' and the Ritz boundary gap. They default to `1e-5` and `1e-4`, respectively,
#' and do not affect backend selection.
#'
#' `dense_n` and `dense_fraction` control the automatic dense fallback. Dense
#' eigenanalysis is selected when `n <= dense_n` or
#' `eig_k >= dense_fraction * n`; their defaults are `100` and `0.5`.
#' `shift_eps` controls the positive shift margin used by iterative methods and
#' defaults to `1e-6`. Backend controls and `shift_eps` use the requested
#' iterative backend instead of the automatic dense fallback. The requested
#' method remains in `eigen$method`, while `eigen$backend$name` identifies the
#' backend actually used.
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
#' # The default return is the embedding matrix.
#' small_iris <- iris[1:75, 1:4]
#' iris_ltsa <- ltsa(
#'   small_iris,
#'   n_neighbors = 12,
#'   nn_method = "exact"
#' )
#'
#' # Request compact eigenanalysis and assembly diagnostics.
#' iris_result <- ltsa(
#'   small_iris,
#'   n_neighbors = 12,
#'   nn_method = "exact",
#'   eig_k = 8,
#'   output = "result"
#' )
#' iris_result$eigen$status
#' iris_result$eigen$messages
#'
#' # Return the assembled unnormalized LTSA matrix.
#' iris_B <- ltsa(
#'   small_iris,
#'   n_neighbors = 12,
#'   nn_method = "exact",
#'   output = "B"
#' )
#' dim(iris_B)
#' @export
ltsa <-
  function(
    X,
    n_neighbors = NULL,
    ndim = 2,
    nn_method = "nnd",
    eig_method = "rspectra",
    eig_k = NULL,
    output = c("embedding", "result", "B"),
    include_B = FALSE,
    include_self = TRUE,
    normalize = FALSE,
    n_threads = 1,
    n_assembly_threads = 1,
    copy_max_mib = 256,
    verbose = FALSE,
    ...
  ) {
    output <- match.arg(output)
    X <- x2m(X)
    neighbor_args <- normalize_ltsa_neighbor_args(
      nn_method = nn_method
    )
    nn_method <- neighbor_args$nn_method
    nn_idx <- neighbor_args$nn_idx

    validated <- validate_ltsa_args(
      X = X,
      n_neighbors = n_neighbors,
      ndim = ndim,
      nn_method = nn_method,
      nn_idx = nn_idx,
      eig_method = eig_method,
      eig_k = eig_k,
      output = output,
      include_B = include_B,
      include_self = include_self,
      normalize = normalize,
      n_threads = n_threads,
      n_assembly_threads = n_assembly_threads,
      copy_max_mib = copy_max_mib,
      verbose = verbose
    )

    eigen_args <- validate_eigen_controls(
      args = list(...),
      eig_method = validated$eig_method,
      output = validated$output
    )

    neighbors <- prepare_ltsa_neighbors(
      X = X,
      n_neighbors = validated$n_neighbors,
      nn_method = validated$nn_method,
      nn_idx = validated$nn_idx,
      include_self = validated$include_self,
      n_threads = validated$n_threads,
      verbose = validated$verbose
    )

    tsmessage("Assembling LTSA matrix", verbose = validated$verbose)
    assembly <- assemble_ltsa_B(
      X = X,
      nn_idx = neighbors$nn_idx,
      ndim = validated$ndim,
      include_self = validated$include_self,
      n_assembly_threads = validated$n_assembly_threads,
      copy_max_mib = validated$copy_max_mib,
      verbose = validated$verbose
    )
    B <- assembly$B

    if (assembly$rank_deficient_count > 0) {
      warning(
        assembly$rank_deficient_count,
        " local neighborhoods had numerical rank below ndim = ",
        validated$ndim,
        "; lower-dimensional local bases were used. Minimum local rank was ",
        assembly$min_local_rank,
        ".",
        call. = FALSE
      )
    }

    if (identical(validated$output, "B")) {
      return(B)
    }

    B_operator <- B
    nullvec <- ltsa_default_null_vector(nrow(B_operator))
    if (validated$normalize) {
      tsmessage("Forming normalized Bsym", verbose = validated$verbose)
      normalized <- ltsa_normalize_sparse_operator(B_operator)
      mass <- normalized$mass
      Dinvs <- normalized$Dinvs
      nullvec <- normalized$nullvec
      B_operator <- normalized$Lsym
    }

    tsmessage("Performing eigenanalysis", verbose = validated$verbose)

    eig_res <- tryCatch(
      {
        ltsa_run_eigenanalysis(
          B = B_operator,
          ndim = validated$ndim,
          eig_method = validated$eig_method,
          eig_k = validated$eig_k,
          eigen_args = eigen_args,
          nullvec = nullvec,
          verbose = validated$verbose
        )
      },
      error = function(e) {
        stop("Eigenanalysis failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    symmetric_embedding <- eig_res$vectors
    embedding <- symmetric_embedding
    if (validated$normalize) {
      embedding <- Dinvs * embedding
    }

    if (identical(validated$output, "embedding")) {
      signal_eigen_status(eig_res$eigen)
      signal_effective_component_status(assembly$diagnostics)
    }

    tsmessage("Finished", verbose = validated$verbose)
    if (identical(validated$output, "embedding")) {
      return(embedding)
    }

    assembly_diagnostics <- assembly$diagnostics %||% list()
    if (!validated$normalize) {
      assembly_diagnostics$component_embedding_overlap <-
        ltsa_component_embedding_overlap(
          eig_res$vectors,
          assembly_diagnostics$component_membership,
          assembly_diagnostics$component_sizes
        )
    }

    eigen_result <- list(
      method = validated$eig_method,
      normalized = isTRUE(validated$normalize),
      eig_k = eig_res$eigen$eig_k,
      values = eig_res$eigen$values,
      ritz_values = eig_res$eigen$ritz_values,
      residuals = eig_res$eigen$residuals,
      rank = eig_res$eigen$rank,
      lambda_max = eig_res$eigen$lambda_max,
      status = eig_res$eigen$status,
      messages = eig_res$eigen$messages,
      backend = eig_res$eigen$backend,
      diagnostics = eig_res$eigen$diagnostics
    )
    if (validated$normalize) {
      eigen_result$normalized_details <- ltsa_normalized_details(
        B = B,
        mass = mass,
        symmetric_embedding = symmetric_embedding,
        embedding = embedding,
        values = eig_res$eigen$values,
        lambda_max = eig_res$eigen$lambda_max,
        reverse_occurrence = assembly$reverse_occurrence,
        component_membership = assembly_diagnostics$component_membership,
        component_sizes = assembly_diagnostics$component_sizes
      )
    }

    result <- list(
      embedding = embedding,
      eigen = eigen_result,
      assembly = lmerge(
        list(
          n_neighbors = as.integer(validated$n_neighbors),
          include_self = isTRUE(validated$include_self),
          neighbor_source = neighbors$source,
          neighbor_elapsed = as.numeric(neighbors$elapsed),
          rank_deficient_count = assembly$rank_deficient_count,
          min_local_rank = assembly$min_local_rank
        ),
        assembly_diagnostics
      )
    )
    if (validated$include_B) {
      result$B <- B
    }
    result
  }

signal_eigen_status <- function(eigen) {
  if (identical(eigen$status, "ok")) {
    return(invisible(NULL))
  }

  message <- paste(eigen$messages, collapse = " ")
  if (identical(eigen$status, "invalid")) {
    stop(
      "LTSA eigenanalysis status is invalid: ",
      message,
      call. = FALSE
    )
  }

  if (identical(eigen$status, "warning")) {
    warning(
      "LTSA eigenanalysis status is warning: ",
      message,
      call. = FALSE
    )
  }

  invisible(NULL)
}

signal_effective_component_status <- function(assembly) {
  if (assembly$component_count <= 1L) {
    return(invisible(NULL))
  }

  displayed_sizes <- head(assembly$component_sizes, 8L)
  size_summary <- paste(displayed_sizes, collapse = ", ")
  if (length(assembly$component_sizes) > length(displayed_sizes)) {
    size_summary <- paste0(size_summary, ", ...")
  }
  warning(
    "The effective-neighborhood co-membership graph has ",
    assembly$component_count,
    " disconnected components (sizes: ",
    size_summary,
    "). Inspect output = \"result\" component diagnostics and increase ",
    "n_neighbors or revise neighbor construction so effective neighborhoods ",
    "overlap.",
    call. = FALSE
  )
  invisible(NULL)
}
