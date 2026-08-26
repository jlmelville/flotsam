# Diagnosing suspicious LTSA results

When an LTSA map looks bent, collapsed, unstable, or unexpectedly
different from another run, first determine which stage supplied the bad
evidence. A neighbor graph can be disconnected, local neighborhoods can
be rank deficient, the eigensolver can be inaccurate, or the requested
dimensions can cut through a poorly separated low-frequency block. Those
cases call for different fixes.

This article gives a question-to-action path through the result object.
See [`?ltsa`](https://jlmelville.github.io/flotsam/reference/ltsa.md)
for the exhaustive argument contract, and see [Exploring LTSA spectral
blocks](https://jlmelville.github.io/flotsam/articles/spectral-blocks.md)
for worked examples of interpreting extra modes.

## Start with the result object

Request `output = "result"`, then inspect status and messages before
tuning anything. This deterministic Swiss-roll example uses exact
neighbors to remove neighbor-search randomness and retains two extra
modes so both the displayed and retained boundaries are available.

``` r

set.seed(20260825)
n <- 1000
phi <- stats::runif(n, min = 1.5 * pi, max = 4.5 * pi)
swiss_roll <- data.frame(
  x = phi * cos(phi),
  y = phi * sin(phi),
  z = stats::runif(n, max = 10)
)

set.seed(20260825)
diagnostic_fit <- ltsa(
  swiss_roll,
  ndim = 2,
  n_neighbors = 15,
  nn_method = "exact",
  eig_method = "rspectra",
  eig_k = 16,
  output = "result",
  spectral_dim = 4,
  n_threads = 1,
  n_assembly_threads = 1
)
```

The displayed `ndim` prefix and the retained `spectral_dim` block have
separate status, messages, and boundary gaps.

``` r

data.frame(
  scope = c("displayed", "retained"),
  status = c(
    diagnostic_fit$eigen$status,
    diagnostic_fit$eigen$spectral$status
  ),
  boundary_gap = formatC(
    c(
      diagnostic_fit$eigen$diagnostics$scaled_boundary_gap,
      diagnostic_fit$eigen$spectral$diagnostics$scaled_boundary_gap
    ),
    digits = 4,
    format = "e"
  ),
  max_scaled_residual = formatC(
    c(
      max(diagnostic_fit$eigen$residuals),
      max(diagnostic_fit$eigen$spectral$residuals)
    ),
    digits = 4,
    format = "e"
  )
)
#>       scope  status boundary_gap max_scaled_residual
#> 1 displayed warning   1.2355e-06          4.4757e-07
#> 2  retained warning   2.9753e-05          4.4757e-07

diagnostic_fit$eigen$messages
#> [1] "Weak Ritz boundary gap after the selected block: 1.235e-06 < 1e-04."
diagnostic_fit$eigen$spectral$messages
#> [1] "Weak Ritz boundary gap after the retained block: 2.975e-05 < 1e-04."
```

Now rule out structural problems in the operator.

``` r

data.frame(
  components = diagnostic_fit$assembly$component_count,
  rank_deficient_neighborhoods =
    diagnostic_fit$assembly$rank_deficient_count,
  minimum_local_rank = diagnostic_fit$assembly$min_local_rank,
  requested_method = diagnostic_fit$eigen$method,
  backend = diagnostic_fit$eigen$backend$name
)
#>   components rank_deficient_neighborhoods minimum_local_rank requested_method
#> 1          1                            0                  2         rspectra
#>    backend
#> 1 rspectra
```

Here the selected scaled residuals are below the default diagnostic
tolerance, the effective-neighborhood graph has one component, and every
local basis has the requested rank. Both statuses instead point to weak
Ritz boundary gaps. The evidence therefore points away from a
disconnected graph, deficient local geometry, or gross eigensolver error
and toward inspecting neighboring modes together. The next action is to
retain a larger block from this operator, as the example already does,
rather than changing `ndim` and rebuilding it.

## Choose the next action from the signal

| Signal | What it supports | Next action |
|----|----|----|
| Backend failure, `status = "invalid"`, or a residual message | The candidate solve is not accurate enough | Name the backend explicitly and adjust its convergence controls; use dense `"eig"` as a small-data reference |
| Weak displayed gap or `near_zero_block_truncated = TRUE` | The requested prefix is not isolated from adjacent modes | Keep `ndim` fixed, retain adjacent modes with `spectral_dim`, and interpret the span rather than individual coordinates |
| `assembly$component_count > 1` | Effective neighborhoods form disconnected groups | Increase or revise neighborhood construction, or analyze the components separately |
| Positive `rank_deficient_count` | Some neighborhoods cannot support the requested local tangent rank | Inspect repeated or degenerate observations and reconsider `n_neighbors` or `ndim` |
| Clean structural and solver diagnostics, but localized or collapsed coordinates | The estimator or selected span may not match the intended structure | Inspect support, compare ordinary and normalized LTSA on one fixed graph, and verify the scientific target |

Connectedness does not prove that a graph is well coupled, and a small
residual does not prove that a chosen boundary is meaningful. Use each
field to rule in or out its own failure class rather than treating
`status` as a single quality score.

### Keep the three dimension controls distinct

- `ndim` sets the local tangent rank, changes the LTSA operator, and
  sets the displayed prefix.
- `spectral_dim` retains more nonconstant modes from that fixed
  operator. It is the inspection control for a weak displayed boundary.
- `eig_k` gives the eigensolver and Rayleigh–Ritz step a wider candidate
  span. It can recover a missing candidate direction but does not decide
  how many modes the result retains.

A stricter tolerance cannot separate a genuinely repeated eigenspace.
When `eigen$diagnostics$near_zero_block_truncated` is `TRUE`, individual
coordinates inside the observed near-zero block are not uniquely
identifiable; reason about the subspace instead.

## When a tighter solve also needs a wider candidate span

The 10,000-sample S-curve-with-a-hole example exercises a different
failure. Its reproduction script generates the public `snedata` dataset,
builds one self-first 15-neighbor graph with one thread, and reuses that
graph for all three solves. It resets the eigensolver seed immediately
before every sampled call.

The frozen full run reports:

| Candidate request | Outcome | Minimum latent-coordinate R-squared |
|----|----|---:|
| `eig_k = 6`, `tol = 1e-6` | Warning result; bent two-dimensional span | 0.0463 |
| `eig_k = 6`, `tol = 1e-8` | RSpectra fails with 0 of 6 candidates converged | — |
| `eig_k = 18`, `tol = 1e-8` | Warning result; expected unrolled span recovered | 0.999999 |

The R-squared witness regresses each recoverable latent coordinate on
the returned two-dimensional span. It is unchanged by sign flips,
rotations, or rescaling of the embedding. The wider strict result
recovers the expected unrolled geometry; the hole is also visible in the
generated figure. Its status remains `"warning"` because the displayed
Ritz boundary is weak, so the example demonstrates successful geometric
recovery without claiming that the two-mode boundary is strongly
isolated.

![With six candidate vectors, one latent direction is largely missed and
the embedding bends.](figures/s-curve-hole-eig-k-6.png)

With six candidate vectors, one latent direction is largely missed and
the embedding bends.

![With eighteen candidate vectors and the tighter tolerance, the same
fixed graph yields the expected unrolled geometry and visible
hole.](figures/s-curve-hole-eig-k-18-tol-1e-8.png)

With eighteen candidate vectors and the tighter tolerance, the same
fixed graph yields the expected unrolled geometry and visible hole.

From the package source tree, after installing the current package and
`snedata`, one command regenerates both figures and stops if the solver,
assembly, or latent-geometry invariants fail:

``` bash
FLOTSAM_NUMERICAL_DIAGNOSTICS_FIGURE_DIR=vignettes/articles/figures \
  Rscript vignettes/articles/numerical-diagnostics-reproduction.R
```

Use the reduced route for a quick no-download execution check. It
deliberately does not overwrite the full article figures because the
scale-dependent six-vector failure does not occur on the smaller sample.

``` bash
FLOTSAM_NUMERICAL_DIAGNOSTICS_SMOKE=true \
  Rscript vignettes/articles/numerical-diagnostics-reproduction.R
```

Show the complete S-curve reproduction script (full route uses 10,000
generated samples)

``` r

# Reproduce the S-curve example in the numerical-diagnostics article.
#
# The full route uses 10,000 generated samples and writes the two article
# figures when FLOTSAM_NUMERICAL_DIAGNOSTICS_FIGURE_DIR is set. Set
# FLOTSAM_NUMERICAL_DIAGNOSTICS_SMOKE=true for a reduced no-download smoke run;
# the smoke route checks the recipe but deliberately does not write article
# figures.

library(flotsam)

if (!requireNamespace("snedata", quietly = TRUE)) {
  stop(
    "Install the article data helper with pak::pak('jlmelville/snedata').",
    call. = FALSE
  )
}

data_seed <- 42L
neighbor_seed <- 42L
solver_seed <- 42L
n_neighbors <- 15L
smoke <- tolower(Sys.getenv(
  "FLOTSAM_NUMERICAL_DIAGNOSTICS_SMOKE",
  unset = "false"
)) %in%
  c("1", "true", "yes")
figure_directory <- Sys.getenv(
  "FLOTSAM_NUMERICAL_DIAGNOSTICS_FIGURE_DIR",
  unset = ""
)

self_first_nnd <- function(X, n_neighbors, seed) {
  set.seed(seed)
  raw <- rnndescent::nnd_knn(
    data = X,
    k = n_neighbors + 1L,
    n_threads = 1L,
    verbose = FALSE
  )

  result <- matrix(NA_integer_, nrow(X), n_neighbors)
  for (row in seq_len(nrow(X))) {
    indices <- raw$idx[row, ]
    distances <- raw$dist[row, ]
    valid <- is.finite(distances) & !is.na(indices) & indices != row
    order_in_row <- order(distances[valid], indices[valid])
    neighbors <- indices[valid][order_in_row]
    neighbors <- neighbors[!duplicated(neighbors)]
    if (length(neighbors) < n_neighbors - 1L) {
      stop("Approximate-neighbor search returned too few neighbors")
    }
    result[row, ] <- c(row, neighbors[seq_len(n_neighbors - 1L)])
  }

  stopifnot(all(result[, 1L] == seq_len(nrow(result))))
  result
}

# snedata constructs x = sin(t), z = sign(t) * (cos(t) - 1), with
# t in [-3 * pi / 2, 3 * pi / 2]. Recover t from the public coordinates so
# the semantic check does not depend on an unexported generator value.
latent_coordinates <- function(data) {
  signed_t <- -sign(data$z)
  cosine <- pmax(-1, pmin(1, 1 - abs(data$z)))
  abs_t <- acos(cosine)
  outside_principal_turn <- signed_t != 0 & sign(data$x) != signed_t
  abs_t[outside_principal_turn] <-
    2 * pi - abs_t[outside_principal_turn]
  t <- signed_t * abs_t

  stopifnot(
    max(abs(sin(t) - data$x)) < 1e-10,
    max(abs(sign(t) * (cos(t) - 1) - data$z)) < 1e-10
  )
  cbind(t = t, height = data$y)
}

latent_r_squared <- function(embedding, latent) {
  vapply(
    seq_len(ncol(latent)),
    function(column) {
      summary(stats::lm(latent[, column] ~ embedding))$r.squared
    },
    numeric(1)
  )
}

run_fit <- function(X, graph, eig_k, tol) {
  set.seed(solver_seed)
  ltsa(
    X,
    n_neighbors = n_neighbors,
    nn_method = graph,
    eig_method = "rspectra",
    eig_k = eig_k,
    output = "result",
    n_threads = 1L,
    n_assembly_threads = 1L,
    tol = tol
  )
}

capture_fit <- function(X, graph, eig_k, tol) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(
      run_fit(X, graph, eig_k, tol),
      error = function(condition) condition
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

fit_summary <- function(name, captured, latent, eig_k, tol) {
  value <- captured$value
  if (inherits(value, "error")) {
    return(data.frame(
      scenario = name,
      outcome = "error",
      eig_k = eig_k,
      tol = tol,
      status = NA_character_,
      min_latent_r_squared = NA_real_,
      components = NA_integer_,
      rank_deficient_neighborhoods = NA_integer_,
      detail = conditionMessage(value)
    ))
  }

  data.frame(
    scenario = name,
    outcome = "result",
    eig_k = eig_k,
    tol = tol,
    status = value$eigen$status,
    min_latent_r_squared = min(
      latent_r_squared(value$embedding, latent)
    ),
    components = value$assembly$component_count,
    rank_deficient_neighborhoods = value$assembly$rank_deficient_count,
    detail = paste(value$eigen$messages, collapse = " | ")
  )
}

draw_embedding <- function(fit, colors, filename, title) {
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(
    file.path(figure_directory, filename),
    width = 1200L,
    height = 900L,
    res = 120L,
    type = "cairo-png",
    bg = "white"
  )
  on.exit(grDevices::dev.off())
  graphics::plot(
    fit$embedding,
    col = colors,
    pch = 16,
    cex = 0.45,
    xlab = "LTSA 1",
    ylab = "LTSA 2",
    main = title
  )
}

set.seed(data_seed)
sample_count <- if (smoke) 2000L else 10000L
s_curve <- snedata::s_curve_hole(n_samples = sample_count)
X <- s_curve[c("x", "y", "z")]
latent <- latent_coordinates(s_curve)
graph <- self_first_nnd(X, n_neighbors, neighbor_seed)

narrow <- capture_fit(X, graph, eig_k = 6L, tol = 1e-6)
strict_narrow <- capture_fit(X, graph, eig_k = 6L, tol = 1e-8)
strict_wide <- capture_fit(X, graph, eig_k = 18L, tol = 1e-8)

summary_table <- rbind(
  fit_summary("narrow", narrow, latent, eig_k = 6L, tol = 1e-6),
  fit_summary(
    "strict_narrow",
    strict_narrow,
    latent,
    eig_k = 6L,
    tol = 1e-8
  ),
  fit_summary(
    "strict_wide",
    strict_wide,
    latent,
    eig_k = 18L,
    tol = 1e-8
  )
)
print(summary_table, row.names = FALSE, digits = 7)

successful <- Filter(
  function(captured) !inherits(captured$value, "error"),
  list(narrow, strict_narrow, strict_wide)
)
stopifnot(
  all(vapply(
    successful,
    function(captured) captured$value$eigen$status != "invalid",
    logical(1)
  )),
  all(vapply(
    successful,
    function(captured) captured$value$assembly$component_count == 1L,
    logical(1)
  )),
  all(vapply(
    successful,
    function(captured) {
      captured$value$assembly$rank_deficient_count == 0L
    },
    logical(1)
  ))
)

if (smoke) {
  stopifnot(
    length(successful) == 3L,
    all(summary_table$min_latent_r_squared > 0.95)
  )
  message("Reduced smoke route passed; article figures were not written.")
} else {
  narrow_score <- min(latent_r_squared(narrow$value$embedding, latent))
  wide_score <- min(latent_r_squared(strict_wide$value$embedding, latent))
  stopifnot(
    inherits(strict_narrow$value, "error"),
    grepl(
      "failed to converge enough LTSA candidate vectors",
      conditionMessage(strict_narrow$value),
      fixed = TRUE
    ),
    narrow_score < 0.25,
    wide_score > 0.99,
    wide_score - narrow_score > 0.75
  )

  if (nzchar(figure_directory)) {
    draw_embedding(
      narrow$value,
      s_curve$color,
      "s-curve-hole-eig-k-6.png",
      "S-curve with hole, eig_k = 6"
    )
    draw_embedding(
      strict_wide$value,
      s_curve$color,
      "s-curve-hole-eig-k-18-tol-1e-8.png",
      "S-curve with hole, eig_k = 18, tol = 1e-8"
    )
  }
  message("Full S-curve semantic checks passed.")
}

versions <- vapply(
  c("flotsam", "snedata", "rnndescent", "RSpectra"),
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
print(versions)
cat("R", as.character(getRversion()), "\n")
```

## Tune the selected backend, not an unspecified route

If residuals or a hard backend failure point to solver accuracy, select
the backend explicitly before passing its controls:

``` r

strict_fit <- ltsa(
  swiss_roll,
  eig_method = "rspectra",
  eig_k = 16,
  output = "result",
  tol = 1e-8,
  maxitr = 5000,
  ncv = 40
)
```

The package default `tol = 1e-6` is looser than RSpectra’s own default.
Moving from `1e-6` to `1e-8` is therefore **tightening** the tolerance.
If that stricter solve cannot return the requested candidate span,
increasing `eig_k` may help, as in the S-curve example; simply
tightening `tol`, enlarging `ncv`, or raising the iteration limit does
not recover a direction absent from the candidate span.

Backend-specific controls passed through `...` are deliberately curated:

- `eig_method = "rspectra"`: `tol`, `maxitr`, and `ncv`;
- `eig_method = "irlba"`: `tol`, `maxit`, and `reorth`; and
- `eig_method = "svdr"`: `tol`, `it`, and `extra`.

`resid_tol` and `gap_tol` change diagnostics, not the selected route.
`dense_n` and `dense_fraction` apply only to `eig_method = "auto"`,
while `shift_eps` applies only to an explicit iterative method. Unknown
controls and controls for another backend are rejected.

For a small diagnostic case, `eig_method = "eig"` (or its `"eigen"`
alias) computes the full dense eigensystem and provides a useful
eigenvalue or eigenspace reference. It is too expensive for large
datasets, and repeated eigenspaces still make individual eigenvectors
non-unique.

## Why `eig_k` can matter

`flotsam` reflects the alignment operator as $`\mu I - B`$, turning the
lowest directions of $`B`$ into the largest algebraic directions
requested from an iterative backend. This shifted-candidate approach was
motivated by practical problems with clustered low eigenvalues and by a
[Spectra discussion](https://github.com/yixuan/spectra/issues/126).

After the candidate solve, a Rayleigh–Ritz step removes the known
constant direction and polishes the requested span. It can resolve mixed
or rotated vectors already present in that span, but it cannot recover
an entirely missing direction. That is why a wider `eig_k` can change a
result even when the requested `ndim` is unchanged.

With `eig_method = "auto"`, `flotsam` uses dense
[`base::eigen()`](https://rdrr.io/r/base/eigen.html) when the input is
small or the candidate request is a large fraction of the sample count;
otherwise it uses RSpectra. `eigen$method` records the requested policy
and `eigen$backend$name` records the backend that actually ran. Even on
the dense route, `eig_k` bounds the candidate span passed to
Rayleigh–Ritz.

Use `output = "B"` to inspect the operator without running final
eigenanalysis. With `normalize = FALSE` this returns the raw alignment
matrix; with `normalize = TRUE` it returns the normalized operator used
by the generalized estimator.

## Normalization changes the estimator

With `normalize = TRUE`, `flotsam` solves

``` math
Bv = \lambda Dv, \qquad D = \operatorname{diag}(B),
```

through the symmetric operator $`D^{-1/2} B D^{-1/2}`$ and maps the
selected vectors back. Standard LTSA uses ordinary orthogonality and
centering; normalized LTSA uses `D`-weighted constraints. This is a
different estimator, not merely a better-conditioned route to the same
coordinates, and `D` is not ordinary graph degree.

When comparing the two estimators, freeze and reuse the neighbor graph,
reset the solver seed before each fit, and keep thread and solver
settings equal. Then compare both spectral evidence and support. The
[spectral-block
article](https://jlmelville.github.io/flotsam/articles/spectral-blocks.md)
applies that protocol to MNIST and Fashion-MNIST.

## Threading and cross-implementation checks

`n_threads` controls computed-neighbor search; `n_assembly_threads`
controls construction of the alignment matrix. Serial settings are the
reproducible default. If outer LTSA assembly and the linked linear
algebra library both use multiple threads, they can oversubscribe the
machine. For MKL, one possible shell configuration is:

``` bash
export MKL_THREADING_LAYER=GNU
export MKL_NUM_THREADS=1
```

The scikit-learn LTSA implementation uses exact neighbor-search
strategies and dense eigenanalysis or ARPACK. `flotsam` can instead use
approximate neighbor descent and shifted RSpectra or `irlba` candidate
solves. A cross-implementation comparison is only diagnostic after
matching the input, neighbor graph, requested dimensions, and tolerances
as closely as the interfaces allow; otherwise a difference does not
identify which stage caused it.
