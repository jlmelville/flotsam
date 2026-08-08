orthonormalize_test_basis <- function(V) {
  qrV <- qr(as.matrix(V))
  qr.Q(qrV)[, seq_len(qrV$rank), drop = FALSE]
}

expect_same_subspace <- function(actual, expected, tolerance = 1e-10) {
  q_actual <- orthonormalize_test_basis(actual)
  q_expected <- orthonormalize_test_basis(expected)
  projector_delta <- tcrossprod(q_actual) - tcrossprod(q_expected)
  projection_distance <- sqrt(sum(projector_delta * projector_delta))

  expect_equal(nrow(q_actual), nrow(q_expected))
  expect_equal(ncol(q_actual), ncol(q_expected))
  expect_lt(projection_distance, tolerance)
}
