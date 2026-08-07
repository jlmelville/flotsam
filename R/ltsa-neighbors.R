normalize_ltsa_neighbor_args <- function(nn_method) {
  if (is.matrix(nn_method)) {
    return(list(nn_method = "nnd", nn_idx = nn_method))
  }

  if (is.list(nn_method) && is.matrix(nn_method$idx)) {
    return(list(nn_method = "nnd", nn_idx = nn_method$idx))
  }

  list(nn_method = nn_method, nn_idx = NULL)
}

prepare_ltsa_neighbors <- function(
  X,
  n_neighbors,
  nn_method,
  nn_idx,
  include_self,
  n_threads,
  verbose = FALSE
) {
  if (!is.null(nn_idx)) {
    tsmessage(
      "Using precomputed nearest-neighbor graph with k = ",
      n_neighbors,
      verbose = verbose
    )
    if (n_threads != 1L) {
      tsmessage(
        "Ignoring n_threads = ",
        n_threads,
        " because precomputed nearest-neighbor graph was supplied",
        verbose = verbose
      )
    }
    return(list(
      nn_idx = nn_idx,
      n_neighbors = n_neighbors,
      source = "precomputed",
      elapsed = NA_real_
    ))
  }

  nn_fun <- switch(
    nn_method,
    exact = rnndescent::brute_force_knn,
    nnd = rnndescent::nnd_knn
  )
  tsmessage(
    "Finding nearest neighbors with method '",
    nn_method,
    "' using n_threads = ",
    n_threads,
    verbose = verbose
  )
  nn_args <- list(
    data = X,
    k = as.integer(min(as.double(n_neighbors) + 1, nrow(X))),
    n_threads = n_threads,
    verbose = FALSE
  )
  elapsed <- system.time({
    nn <- do.call(nn_fun, nn_args)
  })[["elapsed"]]
  nn_idx <- canonicalize_ltsa_computed_neighbors(
    nn_idx = nn$idx,
    nn_dist = nn$dist,
    n_neighbors = n_neighbors,
    include_self = include_self
  )

  list(
    nn_idx = nn_idx,
    n_neighbors = n_neighbors,
    source = nn_method,
    elapsed = elapsed
  )
}

canonicalize_ltsa_computed_neighbors <- function(
  nn_idx,
  nn_dist,
  n_neighbors,
  include_self
) {
  n_obs <- nrow(nn_idx)
  nonself_count <- if (include_self) n_neighbors - 1L else n_neighbors
  row_width <- if (include_self) n_neighbors else n_neighbors + 1L
  canonical_idx <- matrix(NA_integer_, nrow = n_obs, ncol = row_width)

  for (row in seq_len(n_obs)) {
    if (
      anyNA(nn_idx[row, ]) ||
        any(nn_idx[row, ] < 1L | nn_idx[row, ] > n_obs)
    ) {
      stop(
        "Computed nearest-neighbor row ",
        row,
        " contains a missing or out-of-range index",
        call. = FALSE
      )
    }
    candidate_order <- order(nn_dist[row, ], nn_idx[row, ], na.last = TRUE)
    nonself_idx <- nn_idx[row, candidate_order]
    nonself_idx <- unique(nonself_idx[nonself_idx != row])

    if (length(nonself_idx) < nonself_count) {
      stop(
        "Computed nearest-neighbor row ",
        row,
        " cannot provide ",
        nonself_count,
        " unique nonself indices",
        call. = FALSE
      )
    }

    canonical_idx[row, ] <- c(row, head(nonself_idx, nonself_count))
  }

  storage.mode(canonical_idx) <- "integer"
  canonical_idx
}
