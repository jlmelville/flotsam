test_that("sparse normalization implements Jacobi scaling", {
  # Protects the algebra behind the public normalize = TRUE path.
  L <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 3L, 2L, 3L),
    x = c(4, 2, 2, 9, 3, 3, 16),
    dims = c(3L, 3L),
    repr = "C"
  )

  normalized <- flotsam:::ltsa_normalize_sparse_operator(L)
  diagonal <- diag(L)
  Dinvs <- sqrt(1 / diagonal)
  reference <- Matrix::Diagonal(x = Dinvs) %*%
    L %*%
    Matrix::Diagonal(x = Dinvs)

  expect_equal(as.matrix(normalized$normalized_operator), as.matrix(reference))
  expect_equal(normalized$inv_sqrt_diagonal, Dinvs)
  expect_equal(normalized$null_vector, sqrt(diagonal))
})

test_that("normalized public embeddings satisfy the generalized problem", {
  result <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = TRUE,
    normalize = TRUE,
    eig_method = "eig",
    eig_k = 4L,
    output = "result",
    include_B = TRUE
  )

  B <- result$B
  V <- result$embedding
  D <- diag(B)
  values <- result$eigen$values

  expect_equal(
    unname(as.matrix(B %*% V)),
    unname(D * sweep(V, 2L, values, "*")),
    tolerance = 1e-10
  )
  expected_weighted_gram <- diag(ncol(V))
  dimnames(expected_weighted_gram) <- list(colnames(V), colnames(V))
  expect_equal(
    crossprod(V, D * V),
    expected_weighted_gram,
    tolerance = 1e-10
  )
  expect_equal(as.numeric(crossprod(V, D)), rep(0, ncol(V)), tolerance = 1e-10)
})

test_that("normalized embeddings match an independent dense generalized solve", {
  result <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = TRUE,
    normalize = TRUE,
    eig_method = "eig",
    eig_k = 4L,
    output = "result",
    include_B = TRUE
  )

  B <- result$B
  D <- diag(B)
  Dinvs <- sqrt(1 / D)
  symmetric_operator <- Dinvs * sweep(as.matrix(B), 2L, Dinvs, "*")
  dense <- eigen(symmetric_operator, symmetric = TRUE)
  order <- order(dense$values)
  null_direction <- sqrt(D / sum(D))
  null_index <- which.max(abs(crossprod(null_direction, dense$vectors)))
  nonconstant <- order[order != null_index]
  reference <- Dinvs * dense$vectors[, nonconstant[seq_len(2L)], drop = FALSE]

  expect_equal(
    result$eigen$values,
    dense$values[nonconstant[seq_len(2L)]],
    tolerance = 1e-10
  )
  expect_same_subspace(result$embedding, reference, tolerance = 1e-7)
})

test_that("output B returns the operator selected by normalize", {
  args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = TRUE,
    output = "B"
  )

  raw_B <- do.call(ltsa, c(args, list(normalize = FALSE)))
  normalized_B <- do.call(ltsa, c(args, list(normalize = TRUE)))
  inv_sqrt_diagonal <- Matrix::Diagonal(x = 1 / sqrt(diag(raw_B)))
  reference <- inv_sqrt_diagonal %*% raw_B %*% inv_sqrt_diagonal

  expect_equal(
    as.matrix(normalized_B),
    as.matrix(reference),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(as.matrix(normalized_B), as.matrix(raw_B))))

  result <- do.call(
    ltsa,
    c(
      args[names(args) != "output"],
      list(
        normalize = TRUE,
        eig_method = "eig",
        eig_k = 4L,
        output = "result",
        include_B = TRUE
      )
    )
  )
  expect_equal(
    as.matrix(result$B),
    as.matrix(raw_B),
    tolerance = 0
  )
})

test_that("normalized output B still skips eigenanalysis", {
  local_mocked_bindings(
    ltsa_run_eigenanalysis = function(...) {
      stop("eigenanalysis reached", call. = FALSE)
    },
    .package = "flotsam"
  )

  expect_no_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8L,
      normalize = TRUE,
      output = "B"
    )
  )
})
