#' Local Tangent Space Alignment
#'
#' Apply the Local Tangent Space Alignment (LTSA) method (Zhang and Zha, 2004)
#' for dimensionality reduction.
#'
#' @details
#' `ltsa()` builds a symmetric positive-semidefinite LTSA alignment operator
#' `B` from local neighborhoods and returns its lowest nonconstant embedding
#' vectors. The constant vector is a known null direction of `B`, although the
#' nullspace can contain other directions.
#'
#' With `output = "result"`, `eigen$status` and `eigen$messages` report solve
#' problems. If `eigen$diagnostics$near_zero_block_truncated` is `TRUE`, the
#' requested dimensions cut through an observed near-zero eigenspace, so
#' individual coordinates may not be identifiable. With `output = "result"`,
#' set `spectral_dim` above `ndim` to retain more modes from the same operator
#' for subspace-based inspection. Increasing `ndim` instead changes the local
#' tangent dimension and rebuilds the operator. A positive
#' `assembly$rank_deficient_count` and
#' its `assembly$min_local_rank` suggest reconsidering the neighborhoods or
#' input. If `assembly$component_count` exceeds one, use
#' `assembly$component_sizes` and `assembly$component_membership` to reconnect
#' the effective-neighborhood graph or analyze its components separately.
#'
#' @section Normalized LTSA:
#' Let `D = diag(B)`. With `normalize = TRUE`, the eigensolve is performed on
#' `D^(-1/2) B D^(-1/2)` and the selected vectors are mapped back with
#' `D^(-1/2)`. This solves the generalized eigenproblem `B v = lambda D v`.
#' With `output = "B"`, this normalized operator is returned directly and
#' eigenanalysis is skipped. The optional `B` in a detailed result remains the
#' raw, unnormalized alignment matrix.
#'
#' Normalized LTSA is a distinct estimator with `D`-weighted orthogonality and
#' centering constraints. It is not a graph-cut estimator or merely a
#' conditioning route to ordinary LTSA.
#'
#' @section Precomputed neighbor input:
#' `nn_method` may be a precomputed 1-based neighbor index matrix or an object
#' with an `idx` matrix. With `include_self = TRUE`, each row must include its
#' own row index and have `n_neighbors` columns. With `include_self = FALSE`,
#' the first column must contain the row index and is dropped before LTSA
#' assembly, so the matrix must have `n_neighbors + 1` columns.
#'
#' @section Assembly resources:
#' Serial assembly generally uses less temporary storage. Parallel assembly
#' trades additional storage for speed and may compound threaded-BLAS
#' oversubscription.
#'
#' @param X The input data matrix or data frame with one observation per row. If
#'   a data frame is supplied, non-numeric columns are ignored. At least one
#'   numeric column is required.
#' @param n_neighbors The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation. If `NULL`, the
#'   computed-neighbor default is the smaller of `15` and the maximum permitted
#'   by the number of observations and `include_self`; for a precomputed graph
#'   supplied as `nn_method`, it is inferred from the graph. It must be at least
#'   `ndim + 2` so each local projector has a residual direction beyond the
#'   constant and tangent subspaces. Data too small to meet that minimum produce
#'   an error.
#' @param ndim The dimension of the space to embed into.
#' @param spectral_dim Number of nonconstant modes to retain from the fixed LTSA
#'   operator. It defaults to `ndim`. A value greater than `ndim` requires
#'   `output = "result"`; the displayed `embedding` remains `ndim`-dimensional
#'   and the larger block is returned as `spectral_embedding`. An expanded
#'   block can contain at most `n - 3` modes so the candidate span can include
#'   the known null direction and a retained-boundary mode.
#' @param nn_method Method for finding nearest neighbors, or a precomputed
#'   nearest-neighbor graph. Can be one of:
#'   * `"nnd"` Approximate nearest neighbors by Nearest Neighbor Descent.
#'   * `"exact"` Exact nearest neighbors by exhaustively comparing all items.
#'   Slow for large datasets.
#'   * A precomputed nearest-neighbor index matrix.
#'   * A nearest-neighbor result object with an `idx` matrix.
#' @param eig_method How to carry out the final eigendecomposition. Possible
#'   values are:
#'    * `"auto"` Use dense [base::eigen()] when `n <= dense_n` or
#'      `eig_k >= dense_fraction * n`, and otherwise use
#'      [RSpectra::eigs_sym()]. This is the default.
#'    * `"rspectra"` Always use [RSpectra::eigs_sym()]. This method can
#'      report hard backend convergence failures through RSpectra's `nconv`
#'      metadata.
#'    * `"irlba"` Always use [irlba::irlba()] and post-hoc residual
#'      diagnostics.
#'    * `"svdr"` Always use [irlba::svdr()] and post-hoc residual diagnostics.
#'    * `"eig"` or `"eigen"` Use the [base::eigen()] function. This is only
#'      feasible for small datasets and should be used for diagnostic purposes
#'      only. Dense `"eig"` is the better diagnostic reference when algebraic
#'      eigenvalue ordering matters.
#' @param eig_k Number of candidate vectors requested from the final
#'   eigensolver. If `NULL`, the default is
#'   `min(n - 1L, max(12L, spectral_dim + 2L))`, where `n` is the number of
#'   observations. With the default `spectral_dim = ndim`, it must satisfy
#'   `ndim + 1 <= eig_k < n`; the minimum can report that no spare boundary
#'   direction is available. When `spectral_dim > ndim`, it must instead
#'   satisfy `spectral_dim + 2 <= eig_k < n`, leaving room for both the known
#'   null direction and a boundary mode after the retained block. Larger values
#'   give the Rayleigh-Ritz postprocessing a wider candidate span. Dense
#'   eigenanalysis computes the full eigensystem, then retains the lowest
#'   `eig_k` candidate vectors for that postprocessing.
#' @param output What to return:
#'   * `"embedding"` Return the embedding matrix. This is the default.
#'   * `"result"` Return a list containing the embedding, compact eigenanalysis
#'     diagnostics, assembly diagnostics, and optionally `B`. When
#'     `spectral_dim > ndim`, it also contains `spectral_embedding` and retained
#'     block diagnostics in `eigen$spectral`.
#'   * `"B"` Skip final eigenanalysis and return the raw alignment matrix when
#'     `normalize = FALSE`, or the normalized operator described above when
#'     `normalize = TRUE`.
#' @param include_B If `TRUE` and `output = "result"`, include the assembled
#'   unnormalized LTSA matrix `B` in the result object. Ignored for other output
#'   modes.
#' @param include_self Should an item be part of its own neighborhood? This has
#'   a minor effect on most results, but work by Zhang and co-workers (2017)
#'   suggests that this is in effect the main difference between LTSA and the
#'   Hessian Locally Linear Embedding (HLLE) method, so setting this to `FALSE`
#'   may allow emulating the HLLE method.
#' @param normalize If `TRUE`, use the normalized LTSA formulation described in
#'   the "Normalized LTSA" section, including for `output = "B"`. The default
#'   is `FALSE`, which uses the standard LTSA formulation.
#' @param n_threads Nonnegative number of threads for nearest-neighbor
#'   calculation. For reproducibility, set `n_threads = 1` (the default). This
#'   does not control LTSA matrix assembly.
#' @param n_assembly_threads Number of threads to use when constructing the LTSA
#'   alignment matrix `B` after nearest neighbors are computed. The default
#'   `1` requests serial assembly. Values greater than `1` request parallel
#'   assembly, which trades additional temporary storage for speed. See the
#'   [threading note](https://jlmelville.github.io/flotsam/articles/numerical-diagnostics.html#threading)
#'   about avoiding oversubscription when BLAS is also multithreaded.
#' @param verbose If `TRUE` log information about progress to the console.
#' @param ... Optional eigensolver and diagnostic controls. See the
#'   `Eigensolver controls` section.
#'
#' @section Eigensolver controls:
#' Automatic-selection controls:
#'
#' * `"auto"`: `dense_n` and `dense_fraction`. Dense eigenanalysis is selected
#'   when `n <= dense_n` or `eig_k >= dense_fraction * n`; their defaults are
#'   `100` and `0.5`.
#'
#' Backend controls for explicit iterative methods:
#'
#' * `"rspectra"`: `tol`, `maxitr`, and `ncv`. See
#'   \link[RSpectra:eigs_sym]{RSpectra::eigs_sym()} for their meanings. The
#'   package default for `tol` is `1e-6`, deliberately looser than RSpectra's
#'   own `1e-10` default; specify `tol` for a like-for-like comparison.
#' * `"irlba"`: `tol`, `maxit`, and `reorth`.
#'   See \link[irlba:irlba]{irlba::irlba()} for their meanings.
#' * `"svdr"`: `tol`, `it`, and `extra`. See
#'   \link[irlba:svdr]{irlba::svdr()} for their meanings.
#'
#' These are the backend controls exposed by `ltsa()`; other arguments accepted
#' by the underlying packages are not automatically available. Values for the
#' listed controls are passed to the selected backend unchanged.
#'
#' `resid_tol` and `gap_tol` set diagnostic thresholds for candidate residuals
#' and the Ritz boundary gap. They default to `1e-5` and `1e-4`, respectively.
#' They are accepted for every eigensolver mode and do not affect backend
#' selection.
#' Dense `"eig"`/`"eigen"` accepts only these diagnostic controls.
#'
#' `shift_eps` controls the positive shift margin used by iterative methods and
#' defaults to `1e-6`. It is accepted only with an explicit iterative method.
#' Backend controls are accepted only with their explicit backend. The
#' requested policy or canonical method is stored in `eigen$method`
#' (`"eigen"` becomes `"eig"`), while `eigen$backend$name` identifies the
#' backend actually used (`"dense_eigen"`, `"rspectra"`, `"irlba"`, or
#' `"svdr"`). `output = "B"` does not accept eigenanalysis controls in `...`.
#'
#' @return With `output = "embedding"`, an `n` by `ndim` embedding matrix. With
#'   `output = "B"`, the raw LTSA alignment matrix when `normalize = FALSE`, or
#'   the normalized operator supplied to eigenanalysis when `normalize = TRUE`.
#'   With `output = "result"`, a list containing `embedding`, compact `eigen`
#'   solve information, and `assembly` neighbor, rank-deficiency, component,
#'   route, and thread information. When `spectral_dim > ndim`,
#'   `spectral_embedding` contains the retained `n` by `spectral_dim` block,
#'   `eigen` continues to describe the displayed prefix, and `eigen$spectral`
#'   describes the retained block. Both scopes report a scaled boundary gap.
#'   The raw, unnormalized `B` is also included when `include_B = TRUE`.
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
#' iris_embedding <- ltsa(
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
#' iris_result$eigen$method
#' iris_result$eigen$backend$name
#'
#' # Retain extra modes without changing the local tangent dimension.
#' iris_block <- ltsa(
#'   small_iris,
#'   n_neighbors = 12,
#'   nn_method = "exact",
#'   eig_method = "eig",
#'   eig_k = 6,
#'   output = "result",
#'   spectral_dim = 4
#' )
#' dim(iris_block$spectral_embedding)
#' iris_block$eigen$diagnostics$scaled_boundary_gap
#' iris_block$eigen$spectral$diagnostics$scaled_boundary_gap
#'
#' # Return the raw LTSA alignment matrix.
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
    eig_method = "auto",
    eig_k = NULL,
    output = c("embedding", "result", "B"),
    include_B = FALSE,
    include_self = TRUE,
    normalize = FALSE,
    n_threads = 1,
    n_assembly_threads = 1,
    verbose = FALSE,
    ...,
    spectral_dim = ndim
  ) {
    output <- match.arg(output)
    X <- prepare_input_matrix(X)
    neighbor_args <- resolve_neighbor_input(
      nn_method = nn_method
    )
    nn_method <- neighbor_args$nn_method
    nn_idx <- neighbor_args$nn_idx

    validated <- validate_arguments(
      X = X,
      n_neighbors = n_neighbors,
      ndim = ndim,
      spectral_dim = spectral_dim,
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
      verbose = verbose
    )

    eigen_args <- validate_eigen_controls(
      args = list(...),
      eig_method = validated$eig_method,
      output = validated$output
    )

    neighbors <- prepare_neighbors(
      X = X,
      n_neighbors = validated$n_neighbors,
      nn_method = validated$nn_method,
      nn_idx = validated$nn_idx,
      include_self = validated$include_self,
      n_threads = validated$n_threads,
      verbose = validated$verbose
    )

    report_progress("Assembling LTSA matrix", verbose = validated$verbose)
    assembly <- assemble_alignment_matrix(
      X = X,
      nn_idx = neighbors$nn_idx,
      ndim = validated$ndim,
      include_self = validated$include_self,
      n_assembly_threads = validated$n_assembly_threads,
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

    operator <- B
    null_vector <- default_null_vector(nrow(operator))
    if (validated$normalize) {
      report_progress("Forming normalized Bsym", verbose = validated$verbose)
      normalization <- normalize_sparse_operator(operator)
      inv_sqrt_diagonal <- normalization$inv_sqrt_diagonal
      null_vector <- normalization$null_vector
      operator <- normalization$normalized_operator
    }

    if (identical(validated$output, "B")) {
      return(operator)
    }

    report_progress("Performing eigenanalysis", verbose = validated$verbose)

    eigenanalysis <- tryCatch(
      {
        run_eigenanalysis(
          B = operator,
          ndim = validated$ndim,
          spectral_dim = validated$spectral_dim,
          eig_method = validated$eig_method,
          eig_k = validated$eig_k,
          eigen_args = eigen_args,
          null_vector = null_vector,
          verbose = validated$verbose
        )
      },
      error = function(e) {
        stop("Eigenanalysis failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    spectral_dimnames <- list(
      rownames(X),
      paste0("LTSA", seq_len(validated$spectral_dim))
    )
    spectral_embedding <- eigenanalysis$vectors
    if (validated$normalize) {
      spectral_embedding <- inv_sqrt_diagonal * spectral_embedding
    }
    dimnames(spectral_embedding) <- spectral_dimnames
    embedding <- spectral_embedding[,
      seq_len(validated$ndim),
      drop = FALSE
    ]

    if (identical(validated$output, "embedding")) {
      signal_eigen_status(eigenanalysis$eigen)
      signal_effective_component_status(assembly$diagnostics)
    }

    report_progress("Finished", verbose = validated$verbose)
    if (identical(validated$output, "embedding")) {
      return(embedding)
    }

    assembly_diagnostics <- assembly$diagnostics %||% list()

    eigen_summary <- list(
      method = validated$eig_method,
      normalized = isTRUE(validated$normalize),
      eig_k = eigenanalysis$eigen$eig_k,
      values = eigenanalysis$eigen$values,
      ritz_values = eigenanalysis$eigen$ritz_values,
      residuals = eigenanalysis$eigen$residuals,
      rank = eigenanalysis$eigen$rank,
      lambda_max = eigenanalysis$eigen$lambda_max,
      status = eigenanalysis$eigen$status,
      messages = eigenanalysis$eigen$messages,
      backend = eigenanalysis$eigen$backend,
      diagnostics = eigenanalysis$eigen$diagnostics
    )
    result <- list(embedding = embedding)
    if (validated$spectral_dim > validated$ndim) {
      eigen_summary$diagnostics$scaled_boundary_gap <-
        eigenanalysis$displayed_boundary_gap
      spectral_eigen <- eigenanalysis$spectral_eigen
      spectral_eigen$diagnostics$scaled_boundary_gap <-
        eigenanalysis$spectral_boundary_gap
      eigen_summary$spectral <- list(
        dimension = validated$spectral_dim,
        values = spectral_eigen$values,
        residuals = spectral_eigen$residuals,
        rank = spectral_eigen$rank,
        status = spectral_eigen$status,
        messages = spectral_eigen$messages,
        diagnostics = spectral_eigen$diagnostics
      )
      result$spectral_embedding <- spectral_embedding
    }
    result$eigen <- eigen_summary
    result$assembly <- merge_named_lists(
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
