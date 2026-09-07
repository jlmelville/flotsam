all_neighbor_fixture <- function(n) {
  matrix(rep.int(seq_len(n), n), nrow = n, byrow = TRUE)
}

assemble_fixed_ltsa <- function(
  X,
  nn_idx,
  ndim,
  n_assembly_threads = 1L
) {
  flotsam:::assemble_alignment_matrix(
    X = X,
    nn_idx = nn_idx,
    ndim = ndim,
    include_self = TRUE,
    n_assembly_threads = n_assembly_threads
  )
}

ltsa_fixed_B <- function(
  X,
  nn_idx,
  ndim,
  normalize = FALSE,
  n_assembly_threads = 1L
) {
  ltsa(
    X = X,
    nn_method = nn_idx,
    ndim = ndim,
    include_self = TRUE,
    normalize = normalize,
    n_assembly_threads = n_assembly_threads,
    output = "B"
  )
}

normalize_dense_ltsa <- function(B) {
  inv_sqrt_diagonal <- 1 / sqrt(diag(B))
  B * tcrossprod(inv_sqrt_diagonal)
}

expect_residual_projector <- function(W, tangent_rank, tolerance = 1e-11) {
  n <- nrow(W)
  J <- matrix(1 / n, n, n)
  tangent <- diag(n) - J - W

  expect_equal(W, t(W), tolerance = tolerance)
  expect_equal(W %*% W, W, tolerance = tolerance)
  expect_equal(drop(W %*% rep(1, n)), rep(0, n), tolerance = tolerance)
  expect_equal(tangent, t(tangent), tolerance = tolerance)
  expect_equal(tangent %*% tangent, tangent, tolerance = tolerance)
  expect_equal(
    drop(tangent %*% rep(1, n)),
    rep(0, n),
    tolerance = tolerance
  )
  expect_equal(sum(svd(tangent, nu = 0L, nv = 0L)$d > 1e-8), tangent_rank)
}

test_that("decimal-constant neighborhoods have the analytical LTSA operators", {
  n <- 6L
  nn_idx <- all_neighbor_fixture(n)
  J <- matrix(1 / n, n, n)
  W <- diag(n) - J
  ordinary <- n * W
  expected <- list(
    `FALSE` = ordinary,
    `TRUE` = normalize_dense_ltsa(ordinary)
  )

  for (p in c(1L, 7L)) {
    X <- matrix(0.1, nrow = n, ncol = p)
    for (normalize in c(FALSE, TRUE)) {
      for (n_assembly_threads in c(1L, 2L)) {
        expect_warning(
          B <- ltsa_fixed_B(
            X,
            nn_idx,
            ndim = 1L,
            normalize = normalize,
            n_assembly_threads = n_assembly_threads
          ),
          "numerical rank below ndim"
        )
        assembled <- assemble_fixed_ltsa(
          X,
          nn_idx,
          ndim = 1L,
          n_assembly_threads = n_assembly_threads
        )

        expect_equal(
          as.matrix(B),
          expected[[as.character(normalize)]],
          tolerance = 1e-13
        )
        if (!normalize) {
          expect_equal(as.matrix(assembled$B), as.matrix(B), tolerance = 0)
        }
        expect_identical(assembled$rank_deficient_count, n)
        expect_identical(assembled$min_local_rank, 0L)
        if (!normalize) {
          expect_residual_projector(
            as.matrix(B) / n,
            tangent_rank = 0L,
            tolerance = 1e-14
          )
        }
      }
    }
  }
})

test_that("retained tangent bases produce residual projectors", {
  n <- 6L
  nn_idx <- all_neighbor_fixture(n)
  contrast <- c(-2, -1, 0, 0, 1, 2)
  signal <- 0.1 + 2^-30 * contrast
  feature_multipliers <- c(1, -2, 0.5, 4, -0.25, 2, -4)
  fixtures <- list(
    svd = matrix(signal, ncol = 1L),
    gram = outer(signal, feature_multipliers)
  )

  tangent <- tcrossprod(contrast) / drop(crossprod(contrast))
  expected_W <- diag(n) - matrix(1 / n, n, n) - tangent

  for (X in fixtures) {
    for (n_assembly_threads in c(1L, 2L)) {
      ordinary <- ltsa_fixed_B(
        X,
        nn_idx,
        ndim = 1L,
        n_assembly_threads = n_assembly_threads
      )
      normalized <- ltsa_fixed_B(
        X,
        nn_idx,
        ndim = 1L,
        normalize = TRUE,
        n_assembly_threads = n_assembly_threads
      )

      W <- as.matrix(ordinary) / n
      expect_equal(W, expected_W, tolerance = 1e-11)
      expect_residual_projector(W, tangent_rank = 1L)
      expect_equal(
        as.matrix(normalized),
        normalize_dense_ltsa(n * expected_W),
        tolerance = 1e-11
      )
      diagnostics <- assemble_fixed_ltsa(X, nn_idx, ndim = 1L)
      expect_identical(diagnostics$rank_deficient_count, 0L)
      expect_identical(diagnostics$min_local_rank, 1L)
    }
  }
})

test_that("local projectors preserve represented translations", {
  n <- 6L
  nn_idx <- all_neighbor_fixture(n)
  signal <- 0.1 + 2^-20 * c(-2, -1, 0, 0, 1, 2)
  fixtures <- list(
    svd = matrix(signal, ncol = 1L),
    gram = outer(signal, c(1, -2, 0.5, 4, -0.25, 2, -4))
  )

  for (X in fixtures) {
    translated <- X + 1024
    for (normalize in c(FALSE, TRUE)) {
      for (n_assembly_threads in c(1L, 2L)) {
        original <- ltsa_fixed_B(
          X,
          nn_idx,
          ndim = 1L,
          normalize = normalize,
          n_assembly_threads = n_assembly_threads
        )
        shifted <- ltsa_fixed_B(
          translated,
          nn_idx,
          ndim = 1L,
          normalize = normalize,
          n_assembly_threads = n_assembly_threads
        )

        expect_equal(
          as.matrix(shifted),
          as.matrix(original),
          tolerance = 1e-7
        )
      }
    }
  }
})

test_that("separated tangent spaces survive extreme finite scaling", {
  n <- 6L
  nn_idx <- all_neighbor_fixture(n)
  u <- c(-2, -1, 0, 0, 1, 2)
  v <- c(1, -2, 1, 1, -2, 1)
  fixtures <- list(
    svd = cbind(u, v),
    gram = cbind(u, v, u + v, u - v, 2 * u, -3 * v, 0.5 * u + 2 * v)
  )
  tangent <-
    tcrossprod(u) / drop(crossprod(u)) + tcrossprod(v) / drop(crossprod(v))
  expected_W <- diag(n) - matrix(1 / n, n, n) - tangent
  ordinary <- n * expected_W
  expected <- list(
    `FALSE` = ordinary,
    `TRUE` = normalize_dense_ltsa(ordinary)
  )

  for (X in fixtures) {
    for (normalize in c(FALSE, TRUE)) {
      for (n_assembly_threads in c(1L, 2L)) {
        for (scale in c(1e-310, 1e-170, 1e170, 1e300)) {
          scaled <- ltsa_fixed_B(
            scale * X,
            nn_idx,
            ndim = 2L,
            normalize = normalize,
            n_assembly_threads = n_assembly_threads
          )
          expect_equal(
            as.matrix(scaled),
            expected[[as.character(normalize)]],
            tolerance = 1e-11
          )
        }
      }
    }
  }
})

test_that("reused workspaces handle changing neighborhood ranks", {
  block_size <- 6L
  u <- c(-2, -1, 0, 0, 1, 2)
  v <- c(1, -2, 1, 1, -2, 1)
  blocks <- list(
    cbind(u, v),
    matrix(0.1, block_size, 2L),
    cbind(u, 2 * u),
    cbind(v, u)
  )
  X_svd <- do.call(rbind, blocks)
  X_gram <- cbind(X_svd, matrix(0, nrow(X_svd), 5L))
  block_indices <- split(
    seq_len(nrow(X_svd)),
    rep(seq_along(blocks), each = block_size)
  )
  nn_idx <- do.call(
    rbind,
    lapply(
      block_indices,
      function(indices) {
        matrix(rep.int(indices, block_size), nrow = block_size, byrow = TRUE)
      }
    )
  )

  J <- matrix(1 / block_size, block_size, block_size)
  tangent_u <- tcrossprod(u) / drop(crossprod(u))
  tangent_v <- tcrossprod(v) / drop(crossprod(v))
  local_weights <- list(
    diag(block_size) - J - tangent_u - tangent_v,
    diag(block_size) - J,
    diag(block_size) - J - tangent_u,
    diag(block_size) - J - tangent_u - tangent_v
  )

  for (X in list(svd = X_svd, gram = X_gram)) {
    for (normalize in c(FALSE, TRUE)) {
      expected_blocks <- lapply(local_weights, function(W) block_size * W)
      expected_raw <- as.matrix(Matrix::bdiag(expected_blocks))
      expected <- if (normalize) {
        normalize_dense_ltsa(expected_raw)
      } else {
        expected_raw
      }

      for (n_assembly_threads in c(1L, 2L)) {
        expect_warning(
          B <- ltsa_fixed_B(
            X,
            nn_idx,
            ndim = 2L,
            normalize = normalize,
            n_assembly_threads = n_assembly_threads
          ),
          "numerical rank below ndim"
        )
        diagnostics <- assemble_fixed_ltsa(
          X,
          nn_idx,
          ndim = 2L,
          n_assembly_threads = n_assembly_threads
        )

        expect_equal(as.matrix(B), expected, tolerance = 1e-11)
        expect_equal(as.matrix(diagnostics$B), expected_raw, tolerance = 1e-11)
        expect_identical(diagnostics$rank_deficient_count, 2L * block_size)
        expect_identical(diagnostics$min_local_rank, 0L)
      }
    }
  }
})

test_that("Jacobi operators preserve the transformed LTSA null vector", {
  X_svd <- cbind(
    c(-3, -2, -1, 0, 2, 5, 9, 14),
    c(2, -1, 4, 0, 3, 8, 7, 11)
  )
  fixtures <- list(
    svd = X_svd,
    gram = cbind(X_svd, matrix(0, nrow(X_svd), 5L))
  )
  nn_idx <- t(vapply(
    seq_len(nrow(X_svd)),
    function(i) as.integer((i - 1L + 0:3) %% nrow(X_svd) + 1L),
    integer(4L)
  ))

  for (X in fixtures) {
    for (n_assembly_threads in c(1L, 2L)) {
      B <- ltsa_fixed_B(
        X,
        nn_idx,
        ndim = 1L,
        n_assembly_threads = n_assembly_threads
      )
      diagonal <- diag(B)
      expect_gt(max(diagonal) - min(diagonal), 1e-4)

      C <- ltsa_fixed_B(
        X,
        nn_idx,
        ndim = 1L,
        normalize = TRUE,
        n_assembly_threads = n_assembly_threads
      )
      expect_equal(
        drop(C %*% sqrt(diagonal)),
        rep(0, nrow(X)),
        tolerance = 1e-10
      )
    }
  }
})

test_that("non-finite centering arithmetic is not rank deficiency", {
  X <- matrix(
    c(-.Machine$double.xmax, .Machine$double.xmax, 0, 1, 2, 3),
    ncol = 1L
  )
  nn_idx <- all_neighbor_fixture(nrow(X))

  for (normalize in c(FALSE, TRUE)) {
    for (n_assembly_threads in c(1L, 2L)) {
      expect_error(
        ltsa_fixed_B(
          X,
          nn_idx,
          ndim = 1L,
          normalize = normalize,
          n_assembly_threads = n_assembly_threads
        ),
        "non-finite neighborhood arithmetic"
      )
    }
  }
})
