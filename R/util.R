`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

format_timestamp <- function() {
  format(Sys.time(), "%T")
}

# Message with a time stamp, controlled explicitly by verbose or force.
report_progress <-
  function(
    ...,
    verbose = FALSE,
    domain = NULL,
    appendLF = TRUE,
    force = FALSE,
    time_stamp = TRUE
  ) {
    if (force || isTRUE(verbose)) {
      msg <- ""
      if (time_stamp) {
        msg <- paste0(format_timestamp(), " ")
      }
      message(msg, ..., domain = domain, appendLF = appendLF)
      utils::flush.console()
    }
  }

# convert data frame to matrix using numeric columns
prepare_input_matrix <- function(X) {
  if (is.data.frame(X)) {
    numeric_cols <- vapply(X, is.numeric, logical(1))
    if (!any(numeric_cols)) {
      stop("X must contain at least one numeric column", call. = FALSE)
    }
    m <- as.matrix(X[, numeric_cols, drop = FALSE])
  } else if (is.matrix(X)) {
    m <- X
  } else {
    stop("X must be a matrix or data frame", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("X must contain numeric values", call. = FALSE)
  }
  storage.mode(m) <- "double"
  m
}

validate_ltsa_args <- function(
  X,
  n_neighbors,
  ndim,
  nn_method,
  nn_idx,
  eig_method,
  eig_k,
  output,
  include_B,
  include_self,
  normalize,
  n_threads,
  n_assembly_threads,
  verbose
) {
  if (!all(is.finite(X))) {
    stop("X must contain only finite numeric values", call. = FALSE)
  }
  if (nrow(X) < 2) {
    stop("X must contain at least two observations", call. = FALSE)
  }
  if (ncol(X) < 1) {
    stop("X must contain at least one column", call. = FALSE)
  }

  ndim <- check_whole_number(ndim, "ndim", min = 1)
  if (!is.null(n_neighbors)) {
    n_neighbors <- check_whole_number(n_neighbors, "n_neighbors", min = 1)
  }
  n_threads <- check_whole_number(n_threads, "n_threads", min = 0)
  n_assembly_threads <- check_whole_number(
    n_assembly_threads,
    "n_assembly_threads",
    min = 1
  )
  include_self <- check_scalar_logical(include_self, "include_self")
  normalize <- check_scalar_logical(normalize, "normalize")
  include_B <- check_scalar_logical(include_B, "include_B")
  verbose <- check_scalar_logical(verbose, "verbose")

  nn_method <- check_choice(
    nn_method,
    c("exact", "nnd"),
    "nn_method"
  )
  eig_method <- check_choice(
    eig_method,
    c("auto", "rspectra", "irlba", "svdr", "eig", "eigen"),
    "eig_method"
  )
  if (eig_method == "eigen") {
    eig_method <- "eig"
  }

  if (ndim >= nrow(X)) {
    stop("ndim must be less than the number of observations", call. = FALSE)
  }
  eig_k <- ltsa_validate_eig_k(eig_k, ndim = ndim, n = nrow(X))

  if (is.null(nn_idx)) {
    if (is.null(n_neighbors)) {
      n_neighbors <- 15L
    }
  } else {
    nn <- validate_precomputed_neighbors(
      nn_idx = nn_idx,
      n_obs = nrow(X),
      n_neighbors = n_neighbors,
      ndim = ndim,
      include_self = include_self
    )
    nn_idx <- nn$nn_idx
    n_neighbors <- nn$n_neighbors
  }

  validate_neighbor_count(
    n_neighbors = n_neighbors,
    ndim = ndim,
    n_obs = nrow(X),
    include_self = include_self
  )

  list(
    X = X,
    n_neighbors = n_neighbors,
    nn_idx = nn_idx,
    ndim = ndim,
    nn_method = nn_method,
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
}

validate_neighbor_count <- function(
  n_neighbors,
  ndim,
  n_obs,
  include_self
) {
  if (n_neighbors < ndim + 2) {
    stop(
      "n_neighbors must be at least ndim + 2 to leave at least one local ",
      "residual direction beyond the constant and tangent subspaces",
      call. = FALSE
    )
  }

  max_neighbors <- if (include_self) n_obs else n_obs - 1L
  if (n_neighbors > max_neighbors) {
    stop(
      "n_neighbors is too large for the number of observations",
      call. = FALSE
    )
  }

  invisible(n_neighbors)
}

validate_precomputed_neighbors <- function(
  nn_idx,
  n_obs,
  n_neighbors,
  ndim,
  include_self
) {
  if (!is.matrix(nn_idx)) {
    stop(
      "Precomputed nearest-neighbor graph must be a numeric matrix",
      call. = FALSE
    )
  }
  if (!is.numeric(nn_idx)) {
    stop(
      "Precomputed nearest-neighbor graph must be a numeric matrix",
      call. = FALSE
    )
  }
  if (nrow(nn_idx) != n_obs) {
    stop(
      "Precomputed nearest-neighbor graph must have one row per observation in X",
      call. = FALSE
    )
  }
  if (ncol(nn_idx) < 1L) {
    stop(
      "Precomputed nearest-neighbor graph must have at least one column",
      call. = FALSE
    )
  }
  if (anyNA(nn_idx) || !all(is.finite(nn_idx))) {
    stop(
      "Precomputed nearest-neighbor graph must contain only finite whole-number indices",
      call. = FALSE
    )
  }
  if (any(nn_idx != floor(nn_idx))) {
    stop(
      "Precomputed nearest-neighbor graph must contain only whole-number indices",
      call. = FALSE
    )
  }
  if (any(nn_idx < 1 | nn_idx > n_obs)) {
    stop(
      "Precomputed nearest-neighbor graph indices must be between 1 and nrow(X)",
      call. = FALSE
    )
  }

  inferred_n_neighbors <- if (include_self) ncol(nn_idx) else ncol(nn_idx) - 1L
  if (is.null(n_neighbors)) {
    n_neighbors <- as.integer(inferred_n_neighbors)
  }

  expected_ncol <- if (include_self) n_neighbors else n_neighbors + 1L
  if (ncol(nn_idx) != expected_ncol) {
    stop(
      "ncol(nn_method) must match n_neighbors and include_self",
      call. = FALSE
    )
  }

  duplicate_row <- which(vapply(
    seq_len(n_obs),
    function(i) anyDuplicated(nn_idx[i, ]) > 0L,
    logical(1)
  ))
  if (length(duplicate_row) > 0L) {
    stop(
      "Precomputed nearest-neighbor graph row ",
      duplicate_row[1L],
      " contains duplicate indices",
      call. = FALSE
    )
  }

  if (include_self) {
    has_self <- vapply(
      seq_len(n_obs),
      function(i) any(nn_idx[i, ] == i),
      logical(1)
    )
    if (!all(has_self)) {
      stop(
        "Each precomputed nearest-neighbor graph row must contain its own row ",
        "index when include_self is TRUE",
        call. = FALSE
      )
    }
  } else if (!all(nn_idx[, 1L] == seq_len(n_obs))) {
    stop(
      "The first column of the precomputed nearest-neighbor graph must ",
      "contain row self indices when include_self is FALSE",
      call. = FALSE
    )
  }

  storage.mode(nn_idx) <- "integer"
  list(nn_idx = nn_idx, n_neighbors = as.integer(n_neighbors))
}

check_whole_number <- function(x, name, min = 0) {
  if (
    !is.numeric(x) ||
      length(x) != 1 ||
      is.na(x) ||
      !is.finite(x) ||
      x < min ||
      x != floor(x) ||
      x > .Machine$integer.max
  ) {
    stop(name, " must be a whole number >= ", min, call. = FALSE)
  }
  as.integer(x)
}

check_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  x
}

check_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop(
      name,
      " must be one of: ",
      paste(choices, collapse = ", "),
      call. = FALSE
    )
  }
  x <- tolower(x)
  if (!x %in% choices) {
    stop(
      name,
      " must be one of: ",
      paste(choices, collapse = ", "),
      call. = FALSE
    )
  }
  x
}

# Add the named overrides to the base list.
# Use to override default values in base with user-supplied values.
merge_named_lists <- function(base, overrides) {
  for (name in names(overrides)) {
    base[[name]] <- overrides[[name]]
  }
  base
}
