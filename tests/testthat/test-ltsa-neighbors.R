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

exact_nn_idx <- function(X, n_neighbors, include_self) {
  nn <- rnndescent::brute_force_knn(
    data = X,
    k = min(n_neighbors + 1L, nrow(X)),
    n_threads = 0L,
    verbose = FALSE
  )
  flotsam:::canonicalize_ltsa_computed_neighbors(
    nn$idx,
    nn$dist,
    n_neighbors = n_neighbors,
    include_self = include_self
  )
}

capture_ltsa_messages <- function(...) {
  capture.output(invisible(ltsa(...)), type = "message")
}

precomputed_nn_idx <- function(n, n_neighbors, include_self) {
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

test_that("computed neighbors are canonicalized deterministically", {
  # ltsa() does not return its neighbor indices, so test the internal helper
  # directly to verify tie ordering and self handling before assembly.
  # fmt: skip
  nn_idx <- matrix(
    c(
      4L, 1L, 2L, 2L, 3L, 5L,
      4L, 1L, 2L, 6L, 3L, 3L,
      4L, 1L, 2L, 6L, 5L, 5L,
      4L, 1L, 2L, 3L, 5L, 6L,
      1L, 2L, 3L, 4L, 5L, 6L,
      1L, 2L, 3L, 4L, 5L, 6L
    ),
    nrow = 6L,
    byrow = TRUE
  )
  # fmt: skip
  nn_dist <- matrix(
    c(
      0.2, 0.0, 0.1, 0.05, 0.1, 0.1,
      0.2, 0.1, 0.0, 0.10, 0.1, 0.05,
      0.2, 0.1, 0.1, 0.10, 0.1, 0.05,
      0.0, 0.1, 0.1, 0.10, 0.1, 0.10,
      0.1, 0.1, 0.1, 0.10, 0.0, 0.10,
      0.1, 0.1, 0.1, 0.10, 0.1, 0.00
    ),
    nrow = 6L,
    byrow = TRUE
  )

  with_self <- flotsam:::canonicalize_ltsa_computed_neighbors(
    nn_idx,
    nn_dist,
    n_neighbors = 4L,
    include_self = TRUE
  )
  without_self <- flotsam:::canonicalize_ltsa_computed_neighbors(
    nn_idx,
    nn_dist,
    n_neighbors = 4L,
    include_self = FALSE
  )

  expect_identical(with_self[1L, ], c(1L, 2L, 3L, 5L))
  expect_identical(with_self[2L, ], c(2L, 3L, 1L, 6L))
  expect_identical(with_self[3L, ], c(3L, 5L, 1L, 2L))
  expect_identical(without_self[1L, ], c(1L, 2L, 3L, 5L, 4L))
  expect_identical(without_self[2L, ], c(2L, 3L, 1L, 6L, 4L))
  expect_identical(without_self[3L, ], c(3L, 5L, 1L, 2L, 6L))
})

test_that("computed neighbor canonicalization reports insufficient rows", {
  # fmt: skip
  nn_idx <- matrix(
    c(
      1L, 2L, 3L, 3L,
      2L, 1L, 1L, 1L,
      3L, 1L, 2L, 2L,
      4L, 1L, 2L, 2L
    ),
    nrow = 4L,
    byrow = TRUE
  )
  nn_dist <- matrix(0, nrow = 4L, ncol = 4L)

  expect_error(
    flotsam:::canonicalize_ltsa_computed_neighbors(
      nn_idx,
      nn_dist,
      n_neighbors = 3L,
      include_self = TRUE
    ),
    "row 2.*unique"
  )
  expect_error(
    flotsam:::canonicalize_ltsa_computed_neighbors(
      nn_idx,
      nn_dist,
      n_neighbors = 2L,
      include_self = FALSE
    ),
    "row 2.*unique nonself"
  )

  bad_idx <- nn_idx
  bad_idx[1L, ] <- c(1L, 0L, 2L, 3L)
  expect_error(
    flotsam:::canonicalize_ltsa_computed_neighbors(
      bad_idx,
      nn_dist,
      n_neighbors = 3L,
      include_self = TRUE
    ),
    "row 1.*missing or out-of-range"
  )
})

test_that("precomputed exact neighborhoods match computed exact LTSA B", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])

  for (include_self in c(TRUE, FALSE)) {
    nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = include_self)
    computed <- ltsa(
      X,
      n_neighbors = 6L,
      ndim = 2L,
      nn_method = "exact",
      include_self = include_self,
      output = "B"
    )
    precomputed <- ltsa(
      X,
      n_neighbors = NULL,
      ndim = 2L,
      nn_method = nn_idx,
      include_self = include_self,
      output = "B"
    )

    expect_sparse_equivalent(precomputed, computed, tolerance = 0)
  }
})

test_that("default nnd neighbor path returns usable public diagnostics", {
  set.seed(123)
  X <- as.matrix(iris[seq_len(30L), seq_len(4L)])

  expect_warning(
    result <- ltsa(
      X,
      n_neighbors = 8L,
      ndim = 2L,
      eig_method = "eig",
      output = "result",
      n_threads = 0L
    ),
    NA
  )

  expect_equal(dim(result$embedding), c(30L, 2L))
  expect_identical(result$assembly$neighbor_source, "nnd")
  expect_identical(result$assembly$n_neighbors, 8L)
  expect_true(is.finite(result$assembly$neighbor_elapsed))
  expect_gte(result$assembly$neighbor_elapsed, 0)
  expect_identical(result$eigen$method, "eig")
  expect_identical(result$eigen$status, "ok")
  expect_false("B" %in% names(result))
})

test_that("precomputed graph supplied as nn_method skips nearest-neighbor search", {
  set.seed(20)
  X <- matrix(rnorm(8L * 10L), nrow = 8L)
  nn_idx <- precomputed_nn_idx(8L, 4L, include_self = TRUE)

  reference <- flotsam:::assemble_ltsa_B(
    X = X,
    nn_idx = nn_idx,
    ndim = 2L,
    include_self = TRUE
  )
  candidate <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B"
  )

  expect_sparse_equivalent(candidate, reference$B, tolerance = 1e-11)
})

test_that("nn_method can carry a precomputed neighbor graph", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)
  reference <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B"
  )

  from_matrix <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B"
  )
  from_result <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = list(idx = nn_idx),
    include_self = TRUE,
    output = "B"
  )

  expect_sparse_equivalent(from_matrix, reference, tolerance = 0)
  expect_sparse_equivalent(from_result, reference, tolerance = 0)
})

test_that("detailed results report precomputed neighbor diagnostics", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  result <- ltsa(
    X,
    n_neighbors = NULL,
    ndim = 2L,
    nn_method = nn_idx,
    eig_method = "eig",
    include_self = TRUE,
    output = "result"
  )

  expect_identical(result$assembly$n_neighbors, 6L)
  expect_identical(result$assembly$neighbor_source, "precomputed")
  expect_true(is.na(result$assembly$neighbor_elapsed))
})

test_that("effective components use assembled co-membership neighborhoods", {
  X <- matrix(c(0, 1, 2, 3, 10, 11, 12, 13), ncol = 1L)
  # Queries deliberately point across the two groups. With include_self =
  # FALSE, only the final three members of each row are passed to assembly,
  # so directed query-to-neighbor connectivity must not join the groups.
  # fmt: skip
  nn_idx <- matrix(
    c(
      1L, 5L, 6L, 7L,
      2L, 5L, 7L, 8L,
      3L, 6L, 7L, 8L,
      4L, 5L, 6L, 8L,
      5L, 1L, 2L, 3L,
      6L, 1L, 3L, 4L,
      7L, 2L, 3L, 4L,
      8L, 1L, 2L, 4L
    ),
    nrow = 8L,
    byrow = TRUE
  )

  expect_warning(
    result <- ltsa(
      X,
      ndim = 1L,
      nn_method = nn_idx,
      include_self = FALSE,
      eig_method = "eig",
      output = "result"
    ),
    NA
  )

  expect_identical(result$assembly$component_count, 2L)
  expect_identical(result$assembly$component_sizes, c(4L, 4L))
  expect_identical(
    result$assembly$component_membership,
    c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  )
  warning_messages <- character()
  invisible(withCallingHandlers(
    ltsa(
      X,
      ndim = 1L,
      nn_method = nn_idx,
      include_self = FALSE,
      eig_method = "eig"
    ),
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  ))
  expect_true(any(grepl(
    "effective-neighborhood co-membership graph has 2 disconnected",
    warning_messages,
    fixed = TRUE
  )))
})

test_that("precomputed duplicate neighborhoods are rejected by row", {
  X <- matrix(seq_len(8L * 4L), nrow = 8L)

  for (include_self in c(TRUE, FALSE)) {
    bad <- precomputed_nn_idx(8L, 4L, include_self)
    bad[3L, ncol(bad)] <- bad[3L, 2L]

    expect_error(
      ltsa(
        X,
        ndim = 2L,
        nn_method = bad,
        include_self = include_self,
        output = "B"
      ),
      "row 3.*duplicate"
    )
  }
})

test_that("exact search canonicalizes identical observation neighborhoods", {
  X <- matrix(0, nrow = 8L, ncol = 4L)

  for (include_self in c(TRUE, FALSE)) {
    neighbors <- flotsam:::prepare_ltsa_neighbors(
      X = X,
      n_neighbors = 4L,
      nn_method = "exact",
      nn_idx = NULL,
      include_self = include_self,
      n_threads = 0L
    )
    nn_idx <- neighbors$nn_idx
    expected_width <- if (include_self) 4L else 5L

    expect_identical(dim(nn_idx), c(8L, expected_width))
    expect_true(all(nn_idx >= 1L & nn_idx <= nrow(X)))
    expect_true(all(vapply(
      seq_len(nrow(X)),
      function(i) anyDuplicated(nn_idx[i, ]) == 0L,
      logical(1)
    )))
    if (include_self) {
      expect_true(all(vapply(
        seq_len(nrow(X)),
        function(i) i %in% nn_idx[i, ],
        logical(1)
      )))
    } else {
      expect_identical(nn_idx[, 1L], seq_len(nrow(X)))
      expect_true(all(vapply(
        seq_len(nrow(X)),
        function(i) !(i %in% nn_idx[i, -1L]),
        logical(1)
      )))
    }

    expect_warning(
      computed <- ltsa(
        X,
        n_neighbors = 4L,
        ndim = 2L,
        nn_method = "exact",
        include_self = include_self,
        output = "B",
        n_threads = 0L
      ),
      "numerical rank"
    )
    expect_s4_class(computed, "dgCMatrix")
  }
})

test_that("precomputed neighbor graph validation rejects invalid graphs", {
  X <- as.matrix(iris[seq_len(8L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 4L, include_self = TRUE)

  expect_error(
    ltsa(X, ndim = 2L, nn_method = nn_idx[-1L, ], output = "B"),
    "one row per observation"
  )

  bad <- nn_idx
  bad[1L, 2L] <- 9L
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, output = "B"),
    "between 1 and nrow"
  )

  bad <- nn_idx
  bad[1L, 2L] <- NA_integer_
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, output = "B"),
    "finite whole-number"
  )

  bad <- nn_idx
  storage.mode(bad) <- "double"
  bad[1L, 2L] <- bad[1L, 2L] + 0.5
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, output = "B"),
    "whole-number"
  )

  expect_error(
    ltsa(X, ndim = 2L, nn_method = as.vector(nn_idx), output = "B"),
    "nn_method"
  )

  expect_error(
    ltsa(X, ndim = 2L, nn_method = nn_idx[, 1:2], output = "B"),
    "at least ndim \\+ 2"
  )

  expect_error(
    ltsa(
      X,
      n_neighbors = 5L,
      ndim = 2L,
      nn_method = nn_idx,
      output = "B"
    ),
    "ncol\\(nn_method\\)"
  )

  bad <- nn_idx
  bad[1L, ] <- c(2L, 3L, 4L, 5L)
  expect_error(
    ltsa(X, ndim = 2L, nn_method = bad, include_self = TRUE, output = "B"),
    "own row index"
  )

  nn_idx_no_self <- exact_nn_idx(X, n_neighbors = 4L, include_self = FALSE)
  bad <- nn_idx_no_self
  bad[2L, 1L] <- 1L
  expect_error(
    ltsa(
      X,
      ndim = 2L,
      nn_method = bad,
      include_self = FALSE,
      output = "B"
    ),
    "first column"
  )
})

test_that("precomputed graph with n_threads does not warn", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  expect_warning(
    invisible(ltsa(
      X,
      ndim = 2L,
      nn_method = nn_idx,
      include_self = TRUE,
      output = "B",
      n_threads = 2L
    )),
    NA
  )
})

test_that("verbose output describes computed and precomputed neighbor handling", {
  X <- as.matrix(iris[seq_len(18L), seq_len(4L)])
  nn_idx <- exact_nn_idx(X, n_neighbors = 6L, include_self = TRUE)

  expect_length(
    capture_ltsa_messages(
      X,
      n_neighbors = 6L,
      ndim = 2L,
      nn_method = "exact",
      include_self = TRUE,
      output = "B"
    ),
    0L
  )
  expect_length(
    capture_ltsa_messages(
      X,
      ndim = 2L,
      nn_method = nn_idx,
      include_self = TRUE,
      output = "B"
    ),
    0L
  )

  computed_messages <- capture_ltsa_messages(
    X,
    n_neighbors = 6L,
    ndim = 2L,
    nn_method = "exact",
    include_self = TRUE,
    output = "B",
    n_threads = 2L,
    verbose = TRUE
  )
  expect_true(any(grepl(
    "Finding nearest neighbors with method 'exact' using n_threads = 2",
    computed_messages,
    fixed = TRUE
  )))

  default_precomputed_messages <- capture_ltsa_messages(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B",
    verbose = TRUE
  )
  expect_true(any(grepl(
    "Using precomputed nearest-neighbor graph with k = 6",
    default_precomputed_messages,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "Ignoring n_threads",
    default_precomputed_messages,
    fixed = TRUE
  )))

  precomputed_messages <- capture_ltsa_messages(
    X,
    ndim = 2L,
    nn_method = nn_idx,
    include_self = TRUE,
    output = "B",
    n_threads = 2L,
    n_assembly_threads = 3L,
    verbose = TRUE
  )
  expect_true(any(grepl(
    "Using precomputed nearest-neighbor graph with k = 6",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Ignoring n_threads = 2 because precomputed nearest-neighbor graph was supplied",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "assembly workers requested/active: 3/",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Assembling LTSA matrix",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "Getting neighborhoods",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "raw staging memory estimate",
    precomputed_messages,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "row-major|fallback|modeled|memory cap",
    precomputed_messages,
    ignore.case = TRUE
  )))
})
