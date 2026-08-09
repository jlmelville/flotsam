expect_sparse_equivalent <- function(candidate, reference, tolerance = 1e-12) {
  candidate <- Matrix::drop0(candidate)
  reference <- Matrix::drop0(reference)
  expect_s4_class(candidate, "dgCMatrix")
  expect_s4_class(reference, "dgCMatrix")
  expect_identical(candidate@Dim, reference@Dim)
  expect_equal(
    as.matrix(candidate),
    as.matrix(reference),
    tolerance = tolerance
  )
}

make_precomputed_neighbor_indices <- function(n, n_neighbors, include_self) {
  offsets <- if (include_self) {
    seq.int(0L, n_neighbors - 1L)
  } else {
    seq.int(0L, n_neighbors)
  }
  t(vapply(
    seq_len(n),
    function(i) as.integer((i - 1L + offsets) %% n + 1L),
    integer(length(offsets))
  ))
}
