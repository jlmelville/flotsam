assemble_ltsa_B <- function(
  X,
  nn_idx,
  ndim,
  include_self,
  n_assembly_threads = 1L,
  verbose = FALSE
) {
  n_assembly_threads <- check_whole_number(
    n_assembly_threads,
    "n_assembly_threads",
    min = 1
  )
  n <- nrow(X)
  k <- ncol(nn_idx) - as.integer(!include_self)
  weight_idx <- ltsa_effective_weight_idx(nn_idx, include_self)

  if (verbose) {
    tsmessage(
      "Computing local weights and assembling sparse matrix",
      verbose = verbose
    )
  }
  transposed_neighbor_indices <- t(weight_idx)
  components <- if (n_assembly_threads == 1L) {
    ltsa_assemble_local_weights(X, transposed_neighbor_indices, k, ndim)
  } else {
    ltsa_assemble_local_weights_parallel(
      X,
      transposed_neighbor_indices,
      k,
      ndim,
      n_assembly_threads
    )
  }
  components <- lmerge(components, ltsa_effective_components(weight_idx, n))
  log_ltsa_assembly_diagnostics(components, verbose)
  B <- ltsa_components_to_dgCMatrix(components, n)

  list(
    B = B,
    rank_deficient_count = components$rank_deficient_count,
    min_local_rank = components$min_local_rank,
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
    "component_count",
    "component_sizes",
    "component_membership",
    "requested_assembly_threads",
    "effective_assembly_threads"
  )
  components[intersect(fields, names(components))]
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
  invisible(NULL)
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
