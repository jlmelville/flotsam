assemble_ltsa_B <- function(
  X,
  nn_idx,
  ndim,
  include_self,
  n_assembly_threads = 1L,
  copy_max_mib = 256,
  assembly_max_mib = 4096,
  verbose = FALSE
) {
  n_assembly_threads <- check_whole_number(
    n_assembly_threads,
    "n_assembly_threads",
    min = 1
  )
  copy_max_mib <- check_nonnegative_number(
    copy_max_mib,
    "copy_max_mib"
  )
  assembly_max_mib <- check_nonnegative_number(
    assembly_max_mib,
    "assembly_max_mib"
  )
  row_major_copy_max_bytes <- ltsa_mib_to_bytes(
    copy_max_mib,
    "copy_max_mib"
  )

  n <- nrow(X)
  k <- ncol(nn_idx) - as.integer(!include_self)
  memory <- ltsa_assembly_memory_preflight(
    n = n,
    k = k,
    p = ncol(X),
    ndim = ndim,
    requested_assembly_threads = n_assembly_threads,
    include_self = include_self,
    row_major_copy_max_bytes = row_major_copy_max_bytes,
    assembly_max_mib = assembly_max_mib
  )
  weight_idx <- ltsa_effective_weight_idx(nn_idx, include_self)

  if (verbose) {
    tsmessage(
      "Computing local weights and assembling sparse matrix",
      verbose = verbose
    )
  }
  value_nnt <- t(weight_idx)
  components <- if (memory$selected_route == "serial") {
    ltsa_assemble_local_weights(X, value_nnt, k, ndim, row_major_copy_max_bytes)
  } else {
    ltsa_assemble_local_weights_parallel(
      X,
      value_nnt,
      k,
      ndim,
      n_assembly_threads,
      row_major_copy_max_bytes
    )
  }
  if (
    memory$requested_route == "parallel" &&
      memory$selected_route == "serial"
  ) {
    components$requested_assembly_threads <- n_assembly_threads
    components$parallel_fallback_reason <-
      "parallel_estimate_exceeds_assembly_cap"
  }
  components$memory <- memory
  components <- lmerge(
    lmerge(
      components,
      ltsa_local_rank_diagnostics(
        components$local_ranks,
        ndim,
        max_rank = min(ncol(X), k)
      )
    ),
    ltsa_effective_components(weight_idx, n)
  )
  log_ltsa_assembly_diagnostics(components, verbose)
  B <- ltsa_components_to_dgCMatrix(components, n)

  list(
    B = B,
    rank_deficient_count = components$rank_deficient_count,
    min_local_rank = components$min_local_rank,
    reverse_occurrence = tabulate(weight_idx, nbins = n),
    diagnostics = ltsa_assembly_diagnostics(components)
  )
}

ltsa_effective_weight_idx <- function(nn_idx, include_self) {
  if (include_self) {
    nn_idx
  } else {
    nn_idx[, -1L, drop = FALSE]
  }
}

ltsa_assembly_diagnostics <- function(components) {
  fields <- c(
    "assembly_route",
    "local_solver_route",
    "local_rank_histogram",
    "rank_deficient_neighborhood_indices",
    "component_count",
    "component_sizes",
    "component_membership",
    "requested_assembly_threads",
    "effective_assembly_threads",
    "raw_entries_estimate",
    "raw_bytes_estimate",
    "memory",
    "row_major_used",
    "row_major_fallback_reason",
    "parallel_fallback_reason"
  )
  components[intersect(fields, names(components))]
}

ltsa_mib_to_bytes <- function(mib, name) {
  bytes <- mib * 1024 * 1024
  if (!is.finite(bytes)) {
    stop(name, " is too large", call. = FALSE)
  }
  bytes
}

ltsa_assembly_memory_estimates <- function(
  n,
  k,
  p,
  ndim,
  requested_assembly_threads,
  include_self,
  row_major_copy_max_bytes
) {
  estimates <- ltsa_assembly_memory_estimates_cpp(
    n,
    k,
    p,
    ndim,
    requested_assembly_threads,
    include_self,
    row_major_copy_max_bytes
  )
  for (route in c("serial", "parallel")) {
    estimates[[route]]$modeled_peak_mib_bound <-
      estimates[[route]]$modeled_peak_bytes_bound / 1024^2
    estimates[[route]]$final_sparse_output_mib_bound <-
      estimates[[route]]$components_bytes$final_sparse_output_bound / 1024^2
  }
  estimates
}

ltsa_assembly_memory_preflight <- function(
  n,
  k,
  p,
  ndim,
  requested_assembly_threads,
  include_self,
  row_major_copy_max_bytes,
  assembly_max_mib
) {
  cap_bytes <- ltsa_mib_to_bytes(assembly_max_mib, "assembly_max_mib")
  estimates <- ltsa_assembly_memory_estimates(
    n = n,
    k = k,
    p = p,
    ndim = ndim,
    requested_assembly_threads = requested_assembly_threads,
    include_self = include_self,
    row_major_copy_max_bytes = row_major_copy_max_bytes
  )
  requested_route <- if (requested_assembly_threads > 1L) {
    "parallel"
  } else {
    "serial"
  }
  selected_route <- requested_route

  if (requested_route == "serial") {
    if (estimates$serial$modeled_peak_bytes_bound > cap_bytes) {
      ltsa_stop_assembly_memory_cap(
        requested_route,
        estimates,
        assembly_max_mib
      )
    }
  } else if (estimates$parallel$modeled_peak_bytes_bound > cap_bytes) {
    if (estimates$serial$modeled_peak_bytes_bound <= cap_bytes) {
      selected_route <- "serial"
      warning(
        "Parallel LTSA assembly modeled peak bound (",
        format(estimates$parallel$modeled_peak_mib_bound, digits = 5),
        " MiB) exceeds assembly_max_mib = ",
        format(assembly_max_mib, digits = 5),
        "; using serial assembly with modeled peak bound ",
        format(estimates$serial$modeled_peak_mib_bound, digits = 5),
        " MiB.",
        call. = FALSE
      )
    } else {
      ltsa_stop_assembly_memory_cap(
        requested_route,
        estimates,
        assembly_max_mib
      )
    }
  }

  list(
    estimate_kind = "modeled_storage_bound",
    cap_mib = assembly_max_mib,
    cap_bytes = cap_bytes,
    requested_route = requested_route,
    selected_route = selected_route,
    serial = estimates$serial,
    parallel = if (requested_route == "parallel") estimates$parallel else NULL,
    sizeof_bytes = estimates$sizeof_bytes,
    excludes = c(
      "pre-existing X and neighbor input",
      "allocator and STL capacity overhead beyond modeled objects",
      "thread stacks and LAPACK/BLAS internal allocations",
      "R object headers and diagnostic character-string contents",
      "temporary workspaces internal to R primitives such as unique and match",
      "later eigensolver allocations"
    )
  )
}

ltsa_stop_assembly_memory_cap <- function(
  requested_route,
  estimates,
  assembly_max_mib
) {
  requested <- estimates[[requested_route]]
  serial_text <- if (requested_route == "parallel") {
    paste0(
      "; serial modeled peak bound is ",
      format(estimates$serial$modeled_peak_mib_bound, digits = 5),
      " MiB"
    )
  } else {
    ""
  }
  components <- unlist(requested$components_bytes, use.names = TRUE)
  components <- sort(components[components > 0], decreasing = TRUE)
  component_text <- paste0(
    names(utils::head(components, 4L)),
    "=",
    format(utils::head(components, 4L) / 1024^2, digits = 5),
    " MiB",
    collapse = ", "
  )
  stop(
    "LTSA assembly memory preflight failed before assembly staging: requested ",
    requested_route,
    " modeled peak bound is ",
    format(requested$modeled_peak_mib_bound, digits = 5),
    " MiB, above assembly_max_mib = ",
    format(assembly_max_mib, digits = 5),
    serial_text,
    ". Largest modeled components: ",
    component_text,
    ". Increase assembly_max_mib only when sufficient memory is available, ",
    "or reduce the data or neighborhood size.",
    call. = FALSE
  )
}

ltsa_local_rank_diagnostics <- function(local_ranks, ndim, max_rank) {
  local_ranks <- as.integer(local_ranks)
  if (any(local_ranks < 0L | local_ranks > max_rank)) {
    stop("Invalid local ranks returned by assembly", call. = FALSE)
  }

  histogram <- tabulate(local_ranks + 1L, nbins = max_rank + 1L)
  names(histogram) <- as.character(0:max_rank)
  list(
    local_rank_histogram = histogram,
    rank_deficient_neighborhood_indices = which(local_ranks < ndim)
  )
}

ltsa_effective_components <- function(weight_idx, n) {
  parent <- seq_len(n)
  tree_size <- rep.int(1L, n)

  find_root <- function(node) {
    root <- node
    while (parent[[root]] != root) {
      root <- parent[[root]]
    }
    while (parent[[node]] != node) {
      next_node <- parent[[node]]
      parent[[node]] <<- root
      node <- next_node
    }
    root
  }

  union_nodes <- function(left, right) {
    left_root <- find_root(left)
    right_root <- find_root(right)
    if (left_root == right_root) {
      return(invisible(NULL))
    }
    if (tree_size[[left_root]] < tree_size[[right_root]]) {
      tmp <- left_root
      left_root <- right_root
      right_root <- tmp
    }
    parent[[right_root]] <<- left_root
    tree_size[[left_root]] <<-
      tree_size[[left_root]] + tree_size[[right_root]]
    invisible(NULL)
  }

  for (row in seq_len(nrow(weight_idx))) {
    representative <- weight_idx[[row, 1L]]
    if (ncol(weight_idx) > 1L) {
      for (col in 2:ncol(weight_idx)) {
        union_nodes(representative, weight_idx[[row, col]])
      }
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1L))
  component_roots <- unique(roots)
  membership <- match(roots, component_roots)
  sizes <- tabulate(membership, nbins = length(component_roots))
  list(
    component_count = as.integer(length(component_roots)),
    component_sizes = as.integer(sizes),
    component_membership = as.integer(membership)
  )
}

ltsa_component_embedding_overlap <- function(
  embedding,
  membership,
  sizes,
  point_weights = NULL
) {
  embedding <- as.matrix(embedding)
  if (is.null(point_weights)) {
    point_weights <- rep.int(1, nrow(embedding))
  }
  decomposition <- qr(embedding)
  embedding_rank <- decomposition$rank
  component_contrast_rank <- max(0L, length(sizes) - 1L)
  if (embedding_rank == 0L || component_contrast_rank == 0L) {
    return(list(
      principal_angle_cosines = numeric(),
      projection_energy = 0,
      embedding_rank = as.integer(embedding_rank),
      component_contrast_rank = as.integer(component_contrast_rank)
    ))
  }

  embedding_basis <- qr.Q(decomposition)[,
    seq_len(embedding_rank),
    drop = FALSE
  ]
  component_weights <- as.numeric(rowsum(
    point_weights,
    membership,
    reorder = FALSE
  ))
  component_projection <- rowsum(
    sqrt(point_weights) * embedding_basis,
    membership,
    reorder = FALSE
  ) /
    sqrt(component_weights)
  global_direction <- sqrt(component_weights / sum(component_weights))
  component_projection <- component_projection -
    tcrossprod(
      global_direction,
      as.numeric(crossprod(global_direction, component_projection))
    )
  squared_cosines <- pmax(
    eigen(
      crossprod(component_projection),
      symmetric = TRUE,
      only.values = TRUE
    )$values,
    0
  )
  n_angles <- min(component_contrast_rank, embedding_rank)
  cosines <- sqrt(squared_cosines[seq_len(n_angles)])

  list(
    principal_angle_cosines = pmin(as.numeric(cosines), 1),
    projection_energy = as.numeric(sum(squared_cosines)),
    embedding_rank = as.integer(embedding_rank),
    component_contrast_rank = as.integer(component_contrast_rank)
  )
}

log_ltsa_assembly_diagnostics <- function(components, verbose) {
  if (!verbose) {
    return(invisible(NULL))
  }

  diagnostics <- ltsa_assembly_diagnostics(components)
  tsmessage(
    "LTSA assembly route: ",
    diagnostics$assembly_route,
    "; assembly workers requested/active: ",
    diagnostics$requested_assembly_threads,
    "/",
    diagnostics$effective_assembly_threads,
    verbose = verbose
  )
  tsmessage(
    "LTSA row-major Gram used: ",
    diagnostics$row_major_used,
    verbose = verbose
  )
  if (ltsa_log_fallback_reason(diagnostics$row_major_fallback_reason)) {
    tsmessage(
      "LTSA row-major fallback reason: ",
      diagnostics$row_major_fallback_reason,
      verbose = verbose
    )
  }
  if (ltsa_log_fallback_reason(diagnostics$parallel_fallback_reason)) {
    tsmessage(
      "LTSA parallel fallback reason: ",
      diagnostics$parallel_fallback_reason,
      verbose = verbose
    )
  }
  invisible(NULL)
}

ltsa_log_fallback_reason <- function(reason) {
  if (is.null(reason) || length(reason) != 1L || is.na(reason)) {
    return(FALSE)
  }
  nzchar(reason) && !reason %in% c("not_requested", "not_applicable_svd_route")
}

ltsa_components_to_dgCMatrix <- function(components, n) {
  n <- check_whole_number(n, "n", min = 0)
  if (n >= .Machine$integer.max) {
    stop("n is too large for a dgCMatrix", call. = FALSE)
  }
  if (!is.integer(components$i) || !is.integer(components$p)) {
    stop("Invalid sparse component index type", call. = FALSE)
  }
  if (!is.numeric(components$x)) {
    stop("Invalid sparse component value type", call. = FALSE)
  }
  if (length(components$p) != n + 1L) {
    stop("Invalid sparse column pointer length", call. = FALSE)
  }
  if (length(components$i) != length(components$x)) {
    stop("Invalid sparse index/value lengths", call. = FALSE)
  }
  if (length(components$i) > .Machine$integer.max) {
    stop("Too many non-zero slots for a dgCMatrix", call. = FALSE)
  }
  if (any(components$p < 0L) || any(diff(components$p) < 0L)) {
    stop("Invalid sparse column pointers", call. = FALSE)
  }
  if (
    components$p[[1L]] != 0L || components$p[[n + 1L]] != length(components$i)
  ) {
    stop("Invalid sparse column pointer bounds", call. = FALSE)
  }
  if (
    length(components$i) > 0L &&
      (any(components$i < 0L) || any(components$i >= n))
  ) {
    stop("Invalid sparse row indices", call. = FALSE)
  }

  methods::new(
    "dgCMatrix",
    i = components$i,
    p = components$p,
    x = components$x,
    Dim = as.integer(c(n, n))
  )
}
