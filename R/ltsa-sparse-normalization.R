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

  inv_sqrt_diagonal <- sqrt(1 / diagonal)
  L_scaled <- L
  L_scaled@x <- scale_csc_columns(
    L_scaled@p,
    L_scaled@x,
    inv_sqrt_diagonal
  )
  list(
    normalized_operator = inv_sqrt_diagonal * L_scaled,
    inv_sqrt_diagonal = inv_sqrt_diagonal,
    null_vector = 1 / inv_sqrt_diagonal
  )
}
