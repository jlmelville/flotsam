normalize_sparse_operator <- function(B) {
  diagonal <- diag(B)
  if (any(!is.finite(diagonal) | diagonal <= 0)) {
    stop(
      "Cannot normalize the LTSA matrix because its diagonal contains ",
      "non-positive or non-finite entries. Increase n_neighbors or set ",
      "normalize = FALSE.",
      call. = FALSE
    )
  }

  inv_sqrt_diagonal <- sqrt(1 / diagonal)
  scaled <- B
  scaled@x <- scale_csc_columns(
    scaled@p,
    scaled@x,
    inv_sqrt_diagonal
  )
  list(
    normalized_operator = inv_sqrt_diagonal * scaled,
    inv_sqrt_diagonal = inv_sqrt_diagonal,
    null_vector = 1 / inv_sqrt_diagonal
  )
}
