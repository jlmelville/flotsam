test_that("invalid methods are rejected early", {
  expect_error(
    ltsa(iris[1:10, ], nn_method = "bad"),
    "nn_method"
  )
  expect_error(
    ltsa(iris[1:10, ], nn_method = c("exact", "nnd")),
    "nn_method"
  )
  expect_error(
    ltsa(iris[1:10, ], nn_method = NA_character_),
    "nn_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = "bad"),
    "eig_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = "fullsvd"),
    "eig_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = c("eig", "eigen")),
    "eig_method"
  )
  expect_error(
    ltsa(iris[1:10, ], eig_method = NA_character_),
    "eig_method"
  )
  expect_error(
    ltsa(iris[1:10, ], output = "bad"),
    "should be one of"
  )
})

test_that("eigen is accepted as an eig alias", {
  iris10_eig <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "eig"
  )
  iris10_eigen <- ltsa(
    iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8,
    include_self = FALSE,
    eig_method = "eigen"
  )
  expect_equal(iris10_eigen, iris10_eig)
})

test_that("input data must be numeric and finite", {
  expect_error(
    ltsa(seq_len(4)),
    "matrix or data frame"
  )
  expect_error(
    ltsa(data.frame(group = letters[1:4])),
    "at least one numeric column"
  )
  mixed_df <- flotsam:::prepare_input_matrix(
    data.frame(x = 1:3, group = letters[1:3])
  )
  expect_equal(unname(mixed_df), matrix(as.double(1:3), ncol = 1))
  expect_identical(colnames(mixed_df), "x")
  expect_error(
    ltsa(matrix(letters[1:4], ncol = 2)),
    "numeric values"
  )
  expect_error(
    ltsa(matrix(c(1, 2, NA, 4), ncol = 2)),
    "finite"
  )
  expect_error(
    ltsa(matrix(numeric(), nrow = 2, ncol = 0)),
    "at least one column"
  )
  expect_error(
    ltsa(matrix(1, nrow = 1, ncol = 2)),
    "at least two observations"
  )
})

test_that("dimension and neighborhood arguments are validated", {
  expect_error(
    ltsa(iris[1:10, ], ndim = 0),
    "ndim"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 1.5),
    "ndim"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 10),
    "ndim must be less than the number of observations"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 2, eig_k = 2L),
    "eig_k"
  )
  expect_error(
    ltsa(iris[1:10, ], ndim = 2, eig_k = 10L),
    "eig_k"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 2.5),
    "n_neighbors"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 3, ndim = 2),
    "at least ndim \\+ 2"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 11, include_self = TRUE),
    "too large"
  )
  expect_error(
    ltsa(iris[1:10, ], n_neighbors = 10, include_self = FALSE),
    "too large"
  )
  expect_error(
    ltsa(iris[1:10, ], n_threads = -1),
    "n_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_threads = 1.5),
    "n_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = 0),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = -1),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = 1.5),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = NA_real_),
    "n_assembly_threads"
  )
  expect_error(
    ltsa(iris[1:10, ], n_assembly_threads = c(1, 2)),
    "n_assembly_threads"
  )
})

test_that("removed assembly memory controls are rejected through dots", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L
  )

  for (argument in c("copy_max_mib", "assembly_max_mib")) {
    args <- base_args
    args[[argument]] <- 0
    expect_error(
      do.call(ltsa, args),
      paste0(
        "^Argument `",
        argument,
        '` is not supported for eig_method = "auto"$'
      )
    )
  }
})

make_c1_precomputed_neighbors <- function(n, include_self, n_neighbors) {
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

test_that("effective neighborhood size is validated before processing", {
  set.seed(101)
  X <- matrix(stats::rnorm(6L * 4L), nrow = 6L, ncol = 4L)

  for (include_self in c(TRUE, FALSE)) {
    expect_error(
      ltsa(
        X,
        ndim = 2L,
        n_neighbors = 3L,
        nn_method = "exact",
        include_self = include_self,
        output = "B"
      ),
      "at least ndim \\+ 2.*local residual direction"
    )

    expect_error(
      ltsa(
        X,
        ndim = 2L,
        n_neighbors = 3L,
        nn_method = make_c1_precomputed_neighbors(6L, include_self, 3L),
        include_self = include_self,
        output = "B"
      ),
      "at least ndim \\+ 2.*local residual direction"
    )
  }
})

test_that("ndim plus two is the smallest accepted effective neighborhood", {
  set.seed(102)
  X <- matrix(stats::rnorm(6L * 4L), nrow = 6L, ncol = 4L)

  for (include_self in c(TRUE, FALSE)) {
    expect_warning(
      computed <- ltsa(
        X,
        ndim = 2L,
        n_neighbors = 4L,
        nn_method = "exact",
        include_self = include_self,
        output = "B"
      ),
      NA
    )
    expect_s4_class(computed, "dgCMatrix")

    expect_warning(
      precomputed <- ltsa(
        X,
        ndim = 2L,
        n_neighbors = 4L,
        nn_method = make_c1_precomputed_neighbors(6L, include_self, 4L),
        include_self = include_self,
        output = "B"
      ),
      NA
    )
    expect_s4_class(precomputed, "dgCMatrix")
  }
})

test_that("logical arguments must be scalar TRUE or FALSE", {
  expect_error(
    ltsa(iris[1:10, ], include_self = NA),
    "include_self"
  )
  expect_error(
    ltsa(iris[1:10, ], normalize = c(TRUE, FALSE)),
    "normalize"
  )
  expect_error(
    ltsa(iris[1:10, ], include_B = 1),
    "include_B"
  )
  expect_error(
    ltsa(iris[1:10, ], verbose = NA),
    "verbose"
  )
})

test_that("eigen diagnostic tolerances must be finite positive numbers", {
  expect_error(
    ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8L, resid_tol = 0),
    "resid_tol must be a finite positive number"
  )
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8L,
      resid_tol = c(1e-5, 1e-4)
    ),
    "resid_tol must be a finite positive number"
  )
  expect_error(
    ltsa(iris[1:10, ], nn_method = "exact", n_neighbors = 8L, gap_tol = Inf),
    "gap_tol must be a finite positive number"
  )
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8L,
      gap_tol = NA_real_
    ),
    "gap_tol must be a finite positive number"
  )
})

test_that("minimum neighborhood size rejects the former normalized failure", {
  expect_error(
    ltsa(
      iris[1:4, ],
      nn_method = "exact",
      n_neighbors = 3L,
      include_self = TRUE,
      ndim = 2L,
      normalize = TRUE,
      eig_method = "eig"
    ),
    "at least ndim \\+ 2.*local residual direction"
  )
})

test_that("unsupported eigen controls fail concisely", {
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8,
      include_self = FALSE,
      eig_method = "svdr",
      not_an_argument = TRUE
    ),
    '^Argument `not_an_argument` is not supported for eig_method = "svdr"$'
  )
})

test_that("public eigen control names are validated before eigenanalysis", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )

  accepted <- list(
    rspectra = list(
      tol = 1e-8,
      maxitr = 5000L,
      ncv = 8L
    ),
    irlba = list(tol = 1e-8, maxit = 5000L, reorth = TRUE),
    svdr = list(tol = 1e-8, it = 500L, extra = 2L)
  )
  for (method in names(accepted)) {
    for (argument in names(accepted[[method]])) {
      result <- expect_silent(
        do.call(
          ltsa,
          c(
            base_args,
            list(eig_method = method),
            accepted[[method]][argument]
          )
        )
      )
      expect_identical(result$eigen$method, method)
      expect_identical(result$eigen$backend$name, method)
    }
  }

  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "rspectra", tol_resid = 1e-8))
    ),
    '^Argument `tol_resid` is not supported for eig_method = "rspectra"$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "rspectra", maxit = 100L))
    ),
    '^Argument `maxit` is not supported for eig_method = "rspectra"$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "irlba", maxitr = 100L))
    ),
    '^Argument `maxitr` is not supported for eig_method = "irlba"$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "svdr", maxit = 100L))
    ),
    '^Argument `maxit` is not supported for eig_method = "svdr"$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "eig", tol = 1e-8))
    ),
    '^Argument `tol` is not supported for eig_method = "eig"$'
  )
  expect_error(
    do.call(ltsa, c(base_args, list(eig_method = "auto", tol = 1e-8))),
    '^Argument `tol` is not supported for eig_method = "auto"$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "rspectra", dense_n = 0L))
    ),
    '^Argument `dense_n` is not supported for eig_method = "rspectra"$'
  )
  expect_error(
    do.call(ltsa, c(base_args, list(eig_method = "eig", shift_eps = 1e-6))),
    '^Argument `shift_eps` is not supported for eig_method = "eig"$'
  )

  expect_error(
    do.call(
      ltsa,
      c(
        base_args,
        list(
          ndim = 2L,
          include_B = FALSE,
          normalize = FALSE,
          n_threads = 1L,
          n_assembly_threads = 1L,
          verbose = FALSE,
          eig_method = "rspectra",
          1e-8
        )
      )
    ),
    "must be named"
  )
  expect_error(
    do.call(
      ltsa,
      c(
        base_args,
        list(eig_method = "rspectra", tol = 1e-8, tol = 1e-9)
      )
    ),
    '^Argument `tol` is supplied more than once$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "rspectra", retvec = FALSE))
    ),
    '^Argument `retvec` is not supported for eig_method = "rspectra"$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "irlba", nv = 0L))
    ),
    '^Argument `nv` is not supported for eig_method = "irlba"$'
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "svdr", k = 0L))
    ),
    '^Argument `k` is not supported for eig_method = "svdr"$'
  )
})

test_that("backend control values are passed through unchanged", {
  # Direct helper coverage is intentional: asserting these values through
  # ltsa() would assert the installed backends' current validation policies.
  backend_args <- list(
    rspectra = list(tol = 0, maxitr = 2.5, ncv = NA_integer_),
    irlba = list(tol = 0, maxit = 2.5, reorth = 1),
    svdr = list(tol = 0, it = 2.5, extra = 100)
  )

  for (method in names(backend_args)) {
    controls <- flotsam:::validate_eigen_controls(
      args = backend_args[[method]],
      eig_method = method,
      output = "result"
    )
    expect_identical(controls$provider_args, backend_args[[method]])
  }
})

test_that("every eigen control belongs only to its selected mode", {
  control_values <- list(
    dense_n = 0L,
    dense_fraction = 1,
    tol = 1e-8,
    maxitr = 5000L,
    ncv = 8L,
    maxit = 5000L,
    reorth = TRUE,
    it = 500L,
    extra = 2L,
    shift_eps = 1e-6,
    resid_tol = 1e-5,
    gap_tol = 1e-4
  )
  allowed <- list(
    auto = c("dense_n", "dense_fraction", "resid_tol", "gap_tol"),
    rspectra = c(
      "tol",
      "maxitr",
      "ncv",
      "shift_eps",
      "resid_tol",
      "gap_tol"
    ),
    irlba = c("tol", "maxit", "reorth", "shift_eps", "resid_tol", "gap_tol"),
    svdr = c("tol", "it", "extra", "shift_eps", "resid_tol", "gap_tol"),
    eig = c("resid_tol", "gap_tol")
  )

  for (method in names(allowed)) {
    for (control in names(control_values)) {
      args <- control_values[control]
      if (control %in% allowed[[method]]) {
        expect_no_error(
          flotsam:::validate_eigen_controls(args, method, "result")
        )
      } else {
        expect_error(
          flotsam:::validate_eigen_controls(args, method, "result"),
          paste0(
            '^Argument `',
            control,
            '` is not supported for eig_method = "',
            method,
            '"$'
          ),
          info = paste(method, control)
        )
      }
    }
  }
})

test_that("flotsam-owned eigen controls are validated locally", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )

  expect_error(
    do.call(ltsa, c(base_args, list(dense_n = 1.5))),
    "dense_n"
  )
  expect_error(
    do.call(ltsa, c(base_args, list(dense_fraction = 1.1))),
    "dense_fraction"
  )
  expect_error(
    do.call(
      ltsa,
      c(base_args, list(eig_method = "rspectra", shift_eps = 0))
    ),
    "shift_eps"
  )
})

test_that("output B rejects all eigenanalysis controls before assembly", {
  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8L,
      output = "B",
      tol = 1e-8
    ),
    'output = "B" does not accept arguments in `...`'
  )
  expect_error(
    do.call(
      ltsa,
      list(
        X = iris[1:10, ],
        nn_method = "exact",
        n_neighbors = 8L,
        ndim = 2L,
        eig_method = "rspectra",
        eig_k = NULL,
        output = "B",
        include_B = FALSE,
        include_self = TRUE,
        normalize = FALSE,
        n_threads = 1L,
        n_assembly_threads = 1L,
        verbose = FALSE,
        1e-8
      )
    ),
    'output = "B" does not accept arguments in `...`'
  )
})

test_that("contradictory controls fail before neighbor search or assembly", {
  local_mocked_bindings(
    prepare_ltsa_neighbors = function(...) {
      stop("neighbor search reached", call. = FALSE)
    },
    assemble_ltsa_B = function(...) {
      stop("assembly reached", call. = FALSE)
    },
    .package = "flotsam"
  )

  expect_error(
    ltsa(
      iris[1:10, ],
      nn_method = "exact",
      n_neighbors = 8L,
      eig_method = "auto",
      tol = 1e-8
    ),
    '^Argument `tol` is not supported for eig_method = "auto"$'
  )
})

test_that("default auto policy uses dense eigenanalysis on a small problem", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )

  result <- expect_silent(do.call(ltsa, base_args))

  expect_identical(result$eigen$method, "auto")
  expect_identical(result$eigen$backend$name, "dense_eigen")
})

test_that("explicit iterative methods are honored on a small problem", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )

  for (method in c("rspectra", "irlba", "svdr")) {
    set.seed(42)
    result <- expect_silent(
      do.call(ltsa, c(base_args, list(eig_method = method)))
    )
    expect_identical(result$eigen$method, method)
    expect_identical(result$eigen$backend$name, method)
  }
})

test_that("backend-specific arguments pass to the explicit iterative backend", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )
  backend_arguments <- list(
    rspectra = list(tol = 1e-8, maxitr = 5000L, ncv = 8L),
    irlba = list(tol = 1e-8, maxit = 5000L, reorth = TRUE),
    svdr = list(tol = 1e-8, it = 500L, extra = 2L)
  )

  for (method in names(backend_arguments)) {
    for (argument in names(backend_arguments[[method]])) {
      result <- expect_silent(
        do.call(
          ltsa,
          c(
            base_args,
            list(eig_method = method),
            backend_arguments[[method]][argument]
          )
        )
      )
      expect_identical(result$eigen$method, method)
      expect_identical(result$eigen$backend$name, method)
    }
  }
})

test_that("shift_eps is accepted by every explicit iterative backend", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result",
    shift_eps = 1e-6
  )

  for (method in c("rspectra", "irlba", "svdr")) {
    result <- expect_silent(
      do.call(ltsa, c(base_args, list(eig_method = method)))
    )
    expect_identical(result$eigen$method, method)
    expect_identical(result$eigen$backend$name, method)
  }
})

test_that("auto thresholds select dense and RSpectra routes", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )

  dense_n_result <- expect_silent(
    do.call(ltsa, c(base_args, list(dense_n = 10L, dense_fraction = 1)))
  )
  expect_identical(dense_n_result$eigen$method, "auto")
  expect_identical(dense_n_result$eigen$backend$name, "dense_eigen")

  dense_fraction_result <- expect_silent(
    do.call(ltsa, c(base_args, list(dense_n = 0L, dense_fraction = 0.4)))
  )
  expect_identical(dense_fraction_result$eigen$method, "auto")
  expect_identical(dense_fraction_result$eigen$backend$name, "dense_eigen")

  iterative_result <- expect_silent(
    do.call(ltsa, c(base_args, list(dense_n = 0L, dense_fraction = 1)))
  )
  expect_identical(iterative_result$eigen$method, "auto")
  expect_identical(iterative_result$eigen$backend$name, "rspectra")
})

test_that("explicit dense methods report canonical policy and actual backend", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )

  for (method in c("eig", "eigen")) {
    result <- expect_silent(
      do.call(ltsa, c(base_args, list(eig_method = method)))
    )
    expect_identical(result$eigen$method, "eig")
    expect_identical(result$eigen$backend$name, "dense_eigen")
  }
})

test_that("diagnostic arguments do not change automatic routing", {
  base_args <- list(
    X = iris[1:10, ],
    nn_method = "exact",
    n_neighbors = 8L,
    include_self = FALSE,
    eig_k = 4L,
    output = "result"
  )

  for (argument_name in c("resid_tol", "gap_tol")) {
    argument <- setNames(
      list(if (argument_name == "resid_tol") 1e-5 else 1e-4),
      argument_name
    )
    dense_result <- expect_silent(do.call(ltsa, c(base_args, argument)))
    expect_identical(dense_result$eigen$method, "auto")
    expect_identical(dense_result$eigen$backend$name, "dense_eigen")

    iterative_result <- expect_silent(
      do.call(
        ltsa,
        c(base_args, list(dense_n = 0L, dense_fraction = 1), argument)
      )
    )
    expect_identical(iterative_result$eigen$method, "auto")
    expect_identical(iterative_result$eigen$backend$name, "rspectra")
  }
})
