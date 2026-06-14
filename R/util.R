stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
# appears only if called from an environment where a logical verbose = TRUE
# OR force = TRUE
tsmessage <-
  function(...,
           domain = NULL,
           appendLF = TRUE,
           force = FALSE,
           time_stamp = TRUE) {
    verbose <- get0("verbose", envir = sys.parent())

    if (force || (!is.null(verbose) && verbose)) {
      msg <- ""
      if (time_stamp) {
        msg <- paste0(stime(), " ")
      }
      message(msg, ..., domain = domain, appendLF = appendLF)
      utils::flush.console()
    }
  }

# convert data frame to matrix using numeric columns
x2m <- function(X) {
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

validate_ltsa_args <- function(X,
                               n_neighbors,
                               ndim,
                               nn_method,
                               eig_method,
                               include_self,
                               normalize,
                               ret_B,
                               n_threads,
                               verbose) {
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
  n_neighbors <- check_whole_number(n_neighbors, "n_neighbors", min = 1)
  n_threads <- check_whole_number(n_threads, "n_threads", min = 0)

  include_self <- check_scalar_logical(include_self, "include_self")
  normalize <- check_scalar_logical(normalize, "normalize")
  ret_B <- check_scalar_logical(ret_B, "ret_B")
  verbose <- check_scalar_logical(verbose, "verbose")

  nn_method <- check_choice(nn_method, c("exact", "nnd"), "nn_method")
  eig_method <- check_choice(
    eig_method,
    c("rspectra", "irlba", "svdr", "fullsvd", "eig", "eigen"),
    "eig_method"
  )
  if (eig_method == "eigen") {
    eig_method <- "eig"
  }

  if (ndim >= nrow(X)) {
    stop("ndim must be less than the number of observations", call. = FALSE)
  }
  if (n_neighbors <= ndim) {
    stop("n_neighbors must be greater than ndim", call. = FALSE)
  }

  max_neighbors <- if (include_self) nrow(X) else nrow(X) - 1
  if (n_neighbors > max_neighbors) {
    stop(
      "n_neighbors is too large for the number of observations",
      call. = FALSE
    )
  }

  list(
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
}

check_whole_number <- function(x, name, min = 0) {
  if (!is.numeric(x) ||
    length(x) != 1 ||
    is.na(x) ||
    !is.finite(x) ||
    x < min ||
    x != floor(x)) {
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
    stop(name, " must be one of: ", paste(choices, collapse = ", "), call. = FALSE)
  }
  x <- tolower(x)
  if (!x %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "), call. = FALSE)
  }
  x
}

# Add the (named) values in l2 to l1.
# Use to override default values in l1 with user-supplied values in l2
lmerge <- function(l1, l2) {
  for (name in names(l2)) {
    l1[[name]] <- l2[[name]]
  }
  l1
}
