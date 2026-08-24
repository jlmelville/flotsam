spectral_contract_fixture <- function() {
  angle <- 2 * pi * seq.int(0L, 27L) / 28L
  X <- cbind(
    cos(angle),
    sin(angle),
    cos(2 * angle),
    sin(2 * angle),
    cos(3 * angle),
    sin(3 * angle)
  )
  rownames(X) <- paste0("spectral-", seq_len(nrow(X)))
  X
}

dense_spectral_reference <- function(B, spectral_dim, normalize) {
  B <- as.matrix(B)
  if (normalize) {
    diagonal <- diag(B)
    inv_sqrt_diagonal <- sqrt(1 / diagonal)
    operator <- inv_sqrt_diagonal *
      sweep(B, 2L, inv_sqrt_diagonal, "*")
    null_vector <- sqrt(diagonal)
  } else {
    inv_sqrt_diagonal <- NULL
    operator <- B
    null_vector <- rep(1, nrow(B))
  }

  decomposition <- eigen(operator, symmetric = TRUE)
  ordered <- order(decomposition$values)
  unit_null <- null_vector / sqrt(sum(null_vector * null_vector))
  null_index <- which.max(abs(crossprod(unit_null, decomposition$vectors)))
  selected <- ordered[ordered != null_index][seq_len(spectral_dim)]
  vectors <- decomposition$vectors[, selected, drop = FALSE]
  if (normalize) {
    vectors <- inv_sqrt_diagonal * vectors
  }
  vectors
}

test_that("spectral_dim defaults preserve the detailed LTSA contract", {
  X <- spectral_contract_fixture()

  for (normalize in c(FALSE, TRUE)) {
    args <- list(
      X = X,
      n_neighbors = 7L,
      ndim = 2L,
      nn_method = "exact",
      eig_method = "eig",
      eig_k = 6L,
      output = "result",
      normalize = normalize
    )
    implicit <- do.call(ltsa, args)
    explicit <- do.call(ltsa, c(args, list(spectral_dim = 2L)))

    implicit$assembly$neighbor_elapsed <- 0
    explicit$assembly$neighbor_elapsed <- 0
    expect_identical(explicit, implicit)
    expect_named(implicit, c("embedding", "eigen", "assembly"))
    expect_false("spectral_embedding" %in% names(implicit))
    expect_false("spectral" %in% names(implicit$eigen))
    expect_named(
      implicit$eigen$diagnostics,
      "near_zero_block_truncated"
    )
  }
})

test_that("spectral_dim validation preserves a displayed prefix", {
  X <- spectral_contract_fixture()
  common <- list(
    X = X,
    n_neighbors = 7L,
    ndim = 2L,
    nn_method = "exact"
  )

  for (invalid in list(0, 1.5, NA_real_, Inf, c(3, 4))) {
    expect_error(
      do.call(ltsa, c(common, list(spectral_dim = invalid))),
      "spectral_dim"
    )
  }
  expect_error(
    do.call(ltsa, c(common, list(spectral_dim = 1L))),
    "spectral_dim must be at least ndim"
  )
  expect_error(
    do.call(ltsa, c(common, list(spectral_dim = 4L))),
    'spectral_dim can exceed ndim only when output = "result"'
  )
  expect_error(
    do.call(ltsa, c(common, list(spectral_dim = 4L, output = "B"))),
    'spectral_dim can exceed ndim only when output = "result"'
  )
  expect_error(
    do.call(
      ltsa,
      c(
        common,
        list(spectral_dim = 4L, eig_k = 4L, output = "result")
      )
    ),
    "spectral_dim \\+ 1 <= eig_k < n"
  )
  expect_error(
    do.call(
      ltsa,
      c(
        common,
        list(spectral_dim = 4L, eig_k = 5L, output = "result")
      )
    ),
    "spectral_dim \\+ 2 <= eig_k < n"
  )
  expect_error(
    ltsa(
      X[seq_len(6L), ],
      n_neighbors = 4L,
      ndim = 2L,
      nn_method = "exact",
      spectral_dim = 5L,
      output = "result"
    ),
    "spectral_dim must be at most nrow\\(X\\) - 3"
  )
})

test_that("retained modes come from the displayed LTSA operator", {
  X <- spectral_contract_fixture()

  for (normalize in c(FALSE, TRUE)) {
    common <- list(
      X = X,
      n_neighbors = 7L,
      ndim = 2L,
      nn_method = "exact",
      eig_method = "eig",
      eig_k = 6L,
      output = "result",
      include_B = TRUE,
      normalize = normalize
    )
    displayed <- do.call(ltsa, common)
    retained <- do.call(ltsa, c(common, list(spectral_dim = 4L)))

    expect_named(
      retained,
      c("embedding", "spectral_embedding", "eigen", "assembly", "B")
    )
    expect_identical(dim(retained$embedding), c(28L, 2L))
    expect_identical(dim(retained$spectral_embedding), c(28L, 4L))
    expect_identical(
      retained$embedding,
      retained$spectral_embedding[, seq_len(2L), drop = FALSE]
    )
    expect_same_subspace(
      retained$embedding,
      displayed$embedding,
      tolerance = 1e-10
    )
    expect_sparse_equivalent(retained$B, displayed$B, tolerance = 0)
    expect_same_subspace(
      retained$spectral_embedding,
      dense_spectral_reference(
        retained$B,
        spectral_dim = 4L,
        normalize = normalize
      ),
      tolerance = 1e-8
    )

    eigen <- retained$eigen
    expect_length(eigen$values, 2L)
    expect_length(eigen$residuals, 2L)
    expect_identical(eigen$spectral$dimension, 4L)
    expect_length(eigen$spectral$values, 4L)
    expect_length(eigen$spectral$residuals, 4L)
    expect_equal(
      eigen$diagnostics$scaled_boundary_gap,
      (eigen$ritz_values[[3L]] - eigen$ritz_values[[2L]]) /
        max(eigen$lambda_max, 1),
      tolerance = 1e-14
    )
    expect_equal(
      eigen$spectral$diagnostics$scaled_boundary_gap,
      (eigen$ritz_values[[5L]] - eigen$ritz_values[[4L]]) /
        max(eigen$lambda_max, 1),
      tolerance = 1e-14
    )

    if (normalize) {
      diagonal <- diag(retained$B)
      expected_gram <- diag(4L)
      dimnames(expected_gram) <- list(
        colnames(retained$spectral_embedding),
        colnames(retained$spectral_embedding)
      )
      expect_equal(
        crossprod(
          retained$spectral_embedding,
          diagonal * retained$spectral_embedding
        ),
        expected_gram,
        tolerance = 1e-10
      )
      expect_equal(
        unname(as.matrix(retained$B %*% retained$spectral_embedding)),
        unname(
          diagonal *
            sweep(
              retained$spectral_embedding,
              2L,
              eigen$spectral$values,
              "*"
            )
        ),
        tolerance = 1e-10
      )
    } else {
      expect_equal(
        unname(crossprod(retained$spectral_embedding)),
        diag(4L),
        tolerance = 1e-10
      )
      expect_equal(
        unname(colSums(retained$spectral_embedding)),
        rep(0, 4L),
        tolerance = 1e-10
      )
    }
  }
})

test_that("displayed instability remains visible beside a stable retained block", {
  n <- 72L
  angle <- 2 * pi * (seq_len(n) - 1L) / n
  circle <- cbind(cos(angle), sin(angle))

  result <- expect_silent(ltsa(
    circle,
    n_neighbors = 9L,
    ndim = 1L,
    nn_method = "exact",
    eig_method = "eig",
    eig_k = 6L,
    output = "result",
    spectral_dim = 4L
  ))

  expect_identical(result$eigen$status, "warning")
  expect_true(any(grepl(
    "Weak Ritz boundary gap",
    result$eigen$messages,
    fixed = TRUE
  )))
  expect_identical(result$eigen$spectral$status, "ok")
  expect_lt(result$eigen$diagnostics$scaled_boundary_gap, 1e-12)
  expect_gt(
    result$eigen$spectral$diagnostics$scaled_boundary_gap,
    1e-3
  )
})

test_that("displayed and retained messages name their diagnostic scope", {
  n_components <- 6L
  component_size <- 6L
  n <- n_components * component_size
  component <- rep(seq_len(n_components), each = component_size)
  position <- rep(seq_len(component_size), times = n_components)
  X <- cbind(component * 100 + position, 0)
  nn_idx <- t(vapply(
    seq_len(n),
    function(index) {
      group_start <- (component[[index]] - 1L) * component_size
      local_position <- position[[index]]
      as.integer(
        group_start +
          c(
            local_position,
            local_position %% component_size + 1L,
            (local_position - 2L) %% component_size + 1L
          )
      )
    },
    integer(3L)
  ))

  result <- expect_silent(ltsa(
    X,
    ndim = 1L,
    nn_method = nn_idx,
    eig_method = "eig",
    eig_k = 6L,
    output = "result",
    spectral_dim = 3L
  ))

  expect_true(any(grepl(
    "The requested embedding cuts through",
    result$eigen$messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "The retained spectral block cuts through",
    result$eigen$spectral$messages,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "The requested embedding cuts through",
    result$eigen$spectral$messages,
    fixed = TRUE
  )))
})

test_that("all eigensolver routes return the retained public block", {
  methods <- c("auto", "rspectra", "irlba", "svdr", "eig", "eigen")
  X <- spectral_contract_fixture()

  for (method in methods) {
    set.seed(20260826)
    result <- ltsa(
      X,
      n_neighbors = 7L,
      ndim = 2L,
      nn_method = "exact",
      eig_method = method,
      eig_k = 6L,
      output = "result",
      spectral_dim = 4L
    )

    expect_identical(dim(result$embedding), c(28L, 2L))
    expect_identical(dim(result$spectral_embedding), c(28L, 4L))
    expect_identical(
      result$embedding,
      result$spectral_embedding[, seq_len(2L), drop = FALSE]
    )
    expect_false(identical(result$eigen$spectral$status, "invalid"))
  }
})
