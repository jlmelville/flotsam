assemble_ltsa_B <- function(
  X,
  nn_idx,
  ndim,
  include_self,
  n_assembly_threads = 1L,
  copy_max_mib = 256,
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
  row_major_copy_max_bytes <- copy_max_mib * 1024 * 1024

  n <- nrow(X)
  weight_idx <- ltsa_weight_neighborhoods(nn_idx, include_self)
  k <- ncol(weight_idx)

  if (verbose) {
    tsmessage("Computing local weights and assembling sparse matrix")
  }
  value_nnt <- t(weight_idx)
  components <- if (n_assembly_threads <= 1L) {
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
  log_ltsa_assembly_diagnostics(components, verbose)
  B <- ltsa_components_to_dgCMatrix(components, n)

  list(
    B = B,
    rank_deficient_count = components$rank_deficient_count,
    min_local_rank = components$min_local_rank,
    diagnostics = ltsa_assembly_diagnostics(components)
  )
}

ltsa_assembly_diagnostics <- function(components) {
  fields <- c(
    "assembly_route",
    "requested_assembly_threads",
    "effective_assembly_threads",
    "raw_entries_estimate",
    "raw_bytes_estimate",
    "duplicate_fallback_count",
    "row_major_used",
    "row_major_fallback_reason",
    "parallel_fallback_reason"
  )
  components[intersect(fields, names(components))]
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
    diagnostics$effective_assembly_threads
  )
  tsmessage(
    "LTSA raw staged entries estimate: ",
    format_ltsa_count(diagnostics$raw_entries_estimate),
    "; raw staging memory estimate: ",
    format_ltsa_bytes(diagnostics$raw_bytes_estimate)
  )
  tsmessage(
    "LTSA duplicate-neighborhood fallback count: ",
    diagnostics$duplicate_fallback_count,
    "; row-major Gram used: ",
    diagnostics$row_major_used
  )
  if (ltsa_log_fallback_reason(diagnostics$row_major_fallback_reason)) {
    tsmessage(
      "LTSA row-major fallback reason: ",
      diagnostics$row_major_fallback_reason
    )
  }
  if (ltsa_log_fallback_reason(diagnostics$parallel_fallback_reason)) {
    tsmessage(
      "LTSA parallel fallback reason: ",
      diagnostics$parallel_fallback_reason
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

format_ltsa_count <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    return(as.character(x))
  }
  format(
    round(x),
    scientific = FALSE,
    big.mark = ltsa_thousands_sep(),
    decimal.mark = ltsa_decimal_mark(),
    trim = TRUE
  )
}

format_ltsa_bytes <- function(bytes) {
  if (!is.numeric(bytes) || length(bytes) != 1L || is.na(bytes)) {
    return(as.character(bytes))
  }
  if (!is.finite(bytes) || bytes < 0) {
    return(as.character(bytes))
  }

  units <- c("bytes", "KiB", "MiB", "GiB", "TiB")
  unit <- 1L
  value <- bytes
  while (value >= 1024 && unit < length(units)) {
    value <- value / 1024
    unit <- unit + 1L
  }

  if (unit == 1L) {
    paste0(format_ltsa_count(value), " ", units[[unit]])
  } else {
    paste0(
      format(
        round(value, 1L),
        nsmall = 1L,
        scientific = FALSE,
        big.mark = ltsa_thousands_sep(),
        decimal.mark = ltsa_decimal_mark(),
        trim = TRUE
      ),
      " ",
      units[[unit]]
    )
  }
}

ltsa_thousands_sep <- function() {
  sep <- Sys.localeconv()[["thousands_sep"]]
  if (!is.character(sep) || length(sep) != 1L || is.na(sep) || !nzchar(sep)) {
    sep <- if (identical(ltsa_decimal_mark(), ",")) "." else ","
  }
  if (identical(sep, ltsa_decimal_mark())) {
    sep <- if (identical(sep, ",")) "." else ","
  }
  sep
}

ltsa_decimal_mark <- function() {
  mark <- getOption("OutDec", ".")
  if (
    !is.character(mark) || length(mark) != 1L || is.na(mark) || !nzchar(mark)
  ) {
    mark <- Sys.localeconv()[["decimal_point"]]
  }
  if (
    !is.character(mark) || length(mark) != 1L || is.na(mark) || !nzchar(mark)
  ) {
    mark <- "."
  }
  mark
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

ltsa_local_weights <- function(X, nni, ndim) {
  ltsa_local_weights_cpp(X, nni, ndim)
}
