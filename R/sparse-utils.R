ltsa_normalize_sparse_operator <- function(L) {
  diagonal <- diag(L)
  if (any(!is.finite(diagonal) | diagonal <= 0)) {
    stop(
      "Cannot normalize the LTSA matrix because its diagonal contains ",
      "non-positive or non-finite entries. Increase n_neighbors or set ",
      "normalize = FALSE.",
      call. = FALSE
    )
  }

  Dinvs <- sqrt(1 / diagonal)
  L_scaled <- L
  L_scaled@x <- spm_times_scalar(L_scaled@p, L_scaled@x, Dinvs)
  list(
    Lsym = Dinvs * L_scaled,
    Dinvs = Dinvs,
    nullvec = 1 / Dinvs
  )
}
