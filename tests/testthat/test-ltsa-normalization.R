test_that("sparse normalization implements Jacobi scaling", {
  # Protects the algebra behind the public normalize = TRUE path.
  B <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 3L, 2L, 3L),
    x = c(4, 2, 2, 9, 3, 3, 16),
    dims = c(3L, 3L),
    repr = "C"
  )

  normalization <- flotsam:::normalize_sparse_operator(B)
  diagonal <- diag(B)
  inv_sqrt_diagonal <- sqrt(1 / diagonal)
  reference <- Matrix::Diagonal(x = inv_sqrt_diagonal) %*%
    B %*%
    Matrix::Diagonal(x = inv_sqrt_diagonal)

  expect_equal(
    as.matrix(normalization$normalized_operator),
    as.matrix(reference)
  )
  expect_equal(normalization$inv_sqrt_diagonal, inv_sqrt_diagonal)
  expect_equal(normalization$null_vector, sqrt(diagonal))
})

test_that("normalized public path rejects an isolated effective observation", {
  # fmt: skip
  X <- matrix(
    c(
      0, 0, 0, 0,
      1, 0, 0, 1,
      0, 1, 0, 1,
      0, 0, 1, 1,
      1, 1, 0, 2,
      1, 0, 1, 2
    ),
    nrow = 6L,
    byrow = TRUE
  )
  # Observation 6 appears only in the self column that include_self = FALSE
  # removes, leaving a zero diagonal entry in the assembled operator.
  # fmt: skip
  nn_idx <- matrix(
    c(
      1L, 2L, 3L, 4L, 5L,
      2L, 1L, 3L, 4L, 5L,
      3L, 1L, 2L, 4L, 5L,
      4L, 1L, 2L, 3L, 5L,
      5L, 1L, 2L, 3L, 4L,
      6L, 1L, 2L, 3L, 4L
    ),
    nrow = 6L,
    byrow = TRUE
  )
  args <- list(
    X = X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = FALSE,
    output = "B"
  )

  raw_B <- do.call(ltsa, args)
  expect_identical(diag(raw_B)[[6L]], 0)
  expect_error(
    do.call(ltsa, c(args, list(normalize = TRUE))),
    "Cannot normalize the LTSA matrix because its diagonal contains non-positive"
  )
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
  inv_sqrt_diagonal <- sqrt(1 / D)
  symmetric_operator <- inv_sqrt_diagonal *
    sweep(as.matrix(B), 2L, inv_sqrt_diagonal, "*")
  dense <- eigen(symmetric_operator, symmetric = TRUE)
  order <- order(dense$values)
  null_direction <- sqrt(D / sum(D))
  null_index <- which.max(abs(crossprod(null_direction, dense$vectors)))
  nonconstant <- order[order != null_index]
  reference <- inv_sqrt_diagonal *
    dense$vectors[, nonconstant[seq_len(2L)], drop = FALSE]

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
    run_eigenanalysis = function(...) {
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
