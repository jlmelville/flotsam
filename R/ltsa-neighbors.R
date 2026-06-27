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
    tsmessage("Using precomputed nearest-neighbor graph with k = ", n_neighbors)
    if (n_threads > 0L) {
      tsmessage(
        "Ignoring n_threads = ",
        n_threads,
        " because precomputed nearest-neighbor graph was supplied"
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
    n_threads
  )
  nn_args <- list(
    data = X,
    k = ifelse(include_self, n_neighbors, n_neighbors + 1L),
    n_threads = n_threads,
    verbose = FALSE
  )
  elapsed <- system.time({
    nn <- do.call(nn_fun, nn_args)
  })[["elapsed"]]
  nn$dist <- NULL
  mode(nn$idx) <- "integer"

  list(
    nn_idx = nn$idx,
    n_neighbors = n_neighbors,
    source = nn_method,
    elapsed = elapsed
  )
}
