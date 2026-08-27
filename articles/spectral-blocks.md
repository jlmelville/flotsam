# Exploring LTSA spectral blocks

LTSA uses `ndim` twice: once to build the local tangent spaces and again
to choose how many global coordinates to return. Those choices can
differ. A circle is locally one-dimensional, for example, but two
coordinates show its closed path.

`spectral_dim` retains extra coordinates from the same LTSA operator.
This article shows how to inspect them on a circle, COIL-20, MNIST, and
Fashion-MNIST.

Click any figure to open a larger version; click the enlarged image or
press Escape to close it. Keyboard users can focus a figure and press
Enter or Space.

## Keep more modes from the same operator

Given a data matrix `X`, request four low-frequency modes while keeping
a two-dimensional local tangent model and a two-coordinate returned
embedding:

``` r

fit <- ltsa(
  X,
  ndim = 2,
  n_neighbors = 15,
  output = "result",
  spectral_dim = 4
)

dim(fit$embedding)                    # n by 2
dim(fit$spectral_embedding)           # n by 4
identical(
  fit$embedding,
  fit$spectral_embedding[, 1:2, drop = FALSE]
)
```

The calculation follows this pipeline:

``` text
data + graph + ndim
        |
        v
  fixed operator B
        |
        v
 eig_k candidate modes
        |
        v
spectral_dim retained modes
        |
        v
 first ndim displayed
```

| Control | Changes `B`? | Purpose |
|----|:--:|----|
| `ndim` | yes | Sets the local tangent rank and displayed prefix |
| `spectral_dim` | no | Retains extra nonconstant modes from the fixed operator |
| `eig_k` | no | Gives the eigensolver a numerical candidate span |

Most users can leave `eig_k` at its default. When retaining an expanded
block, a manual request must be at least `spectral_dim + 2`: one place
for the known constant direction and another for the mode just beyond
the retained block.

The displayed and retained boundaries have separate diagnostics:

``` r

fit$eigen$diagnostics$scaled_boundary_gap
fit$eigen$spectral$diagnostics$scaled_boundary_gap
fit$eigen$status
fit$eigen$spectral$status
```

A boundary gap compares the last included mode with the next one. A
large gap isolates the chosen subspace. A small gap suggests
interpreting the neighboring modes together and checking their geometry
and support.

## A short inspection routine

Three questions organize the diagnostics:

1.  **Is the operator structurally sound?** Check graph connectivity and
    local ranks.
2.  **Was the spectrum computed reliably?** Check status, messages, and
    residuals.
3.  **Is the displayed view well determined and useful?** Check boundary
    gaps, inspect a small retained block, and measure how broadly each
    mode is supported.

When comparing stochastic fits, reuse the same fixed in-memory neighbor
matrix, reset the eigensolver seed immediately before each fit, and keep
thread settings equal.

## A circle needs a spectral pair

A circle has a one-dimensional tangent at every point. A single
real-valued coordinate cannot follow the whole loop continuously, so the
first two low-frequency modes work as a pair. Their signs and
orientation can change between fits; treat the two-dimensional span as
the stable unit.

The example uses exact neighbors and dense eigenanalysis, making the
repeated pairs easy to see.

``` r

theta <- seq(0, 2 * pi, length.out = 721)[-721]
circle <- cbind(x = cos(theta), y = sin(theta))

circle_fit <- ltsa(
  circle,
  ndim = 1,
  n_neighbors = 9,
  nn_method = "exact",
  eig_method = "eig",
  eig_k = 8,
  output = "result",
  spectral_dim = 4
)

circle_fit$eigen$values
#> [1] 4.464673e-07
circle_fit$eigen$spectral$values
#> [1] 4.464673e-07 4.464673e-07 7.139639e-06 7.139639e-06
circle_fit$eigen$status
#> [1] "warning"
circle_fit$eigen$messages
#> [1] "Weak Ritz boundary gap after the selected block: 1.293e-17 < 1e-04."
```

![The input circle, its one-dimensional LTSA display, the first two
modes from the same operator, and the low spectrum. The mode pair
follows the complete periodic
progression.](spectral-blocks_files/figure-html/circle-figure-1.png)

The input circle, its one-dimensional LTSA display, the first two modes
from the same operator, and the low spectrum. The mode pair follows the
complete periodic progression.

Setting `ndim = 2` would build a two-dimensional local tangent model.
This fit instead keeps `ndim = 1` and uses `spectral_dim = 4`: the first
two retained modes reveal the global coordinate pair, while the next two
expose its neighboring spectral context.

## COIL-20: a real periodic example

[COIL-20](https://cave.cs.columbia.edu/repository/COIL-20) contains 72
views of each object as it rotates. One object traces a closed, locally
one-dimensional curve through pixel space, so the circle example
predicts a spectral pair.

At `k = 7`, mixing objects produces a disconnected neighbor graph, so we
fit objects separately. Objects 4 and 1 were fixed by an earlier
topology-only screen run before any embedding. The screen required
connected graphs, bidirectional local pose coverage, wraparound support,
and at most one-sixth pose-order shortcut edges. Nine objects qualified,
and ordering them by shortcut rate selected objects 4 and 1. The
companion below reproduces these fixed-object fits and their spectral
claims, not the earlier 20-object screen. Each fit uses ordinary LTSA
with `ndim = 1`, exact neighbors, and four retained modes. The
coordinates come from pixels; pose supplies the point colors and
thumbnail labels.

![COIL-20 object 4. Modes 1 and 2 form a pose-ordered orbit, with
thumbnails placed in the same radial direction as their
observations.](figures/ltsa-spectral-block-coil-object-4.png)

COIL-20 object 4. Modes 1 and 2 form a pose-ordered orbit, with
thumbnails placed in the same radial direction as their observations.

![COIL-20 object 1. The path is more distorted and still follows the
object’s rotation.](figures/ltsa-spectral-block-coil-object-1.png)

COIL-20 object 1. The path is more distorted and still follows the
object’s rotation.

The gap after mode 2 is much larger than the gap between modes 1 and 2.
For objects 4 and 1, the mode 1–2 gaps are `0.000108` and `0.000181`,
while the mode 2–3 gaps are `0.00250` and `0.00215`—about 23 and 12
times larger. This is the spectral pattern expected from a mode pair.
The second mode completes the rotational orbit already present in the
fixed operator.

## Check whether a mode represents many observations

A low-energy mode can be dominated by a few observations. For a
Euclidean-normalized coordinate vector `v`, the effective support
fraction is

``` math
\frac{1}{n\sum_i v_i^4}.
```

This is the normalized participation ratio. It equals $`1/n`$ when one
observation carries the mode, $`r/n`$ when the energy is spread equally
across exactly $`r`$ observations, and 1 when it is spread equally
across all observations. Here, a mode’s energy means its squared
coordinate values, `v^2`:

``` r

effective_support_fraction <- function(v) {
  v <- v / sqrt(sum(v^2))
  1 / (length(v) * sum(v^4))
}
```

This scalar shows whether a coordinate is carried by many observations
or only a few. The plots then reveal the variation organized by a
broadly supported coordinate.

For a retained block, first orthonormalize its columns and let each
observation’s block leverage be the sum of its squared entries across
that basis. Normalizing those leverages to sum to one gives the same
effective support formula for the whole span. This block support is
unchanged by rotating or rescaling coordinates within the span. It does
not establish where the span should end, so read it together with the
boundary gap.

## MNIST: a complementary view

The 7,877 MNIST digit-1 images show how localization changes the
interpretation of extra modes. Ordinary and normalized LTSA reuse the
same fixed neighbor matrix, `ndim = 2`, `k = 15`, and four retained
modes. Normalized LTSA changes the weighting through a generalized
eigenproblem, giving a distinct estimator over the same graph.

The ordinary first three modes are concentrated on individual
observations. Normalized modes 1 and 2 are broadly supported and arrange
handwriting by slant, curvature, stroke width, hooks, and serifed or
branched forms.

At the displayed boundary, the ordinary modes 1–2 span has block support
`0.000261` with a weak mode 2–3 gap of `2.81e-5`. The normalized span
has block support `0.480` with a mode 2–3 gap of `0.00201`. The
span-level diagnostic therefore agrees with the individual modes: the
ordinary view is almost point-localized, while the normalized view is
broadly supported.

![The ordinary MNIST modes 1–2 span is dominated by two atypical
observations; leverage then falls by several orders of
magnitude.](figures/ltsa-spectral-block-mnist-ordinary-high-leverage.png)

The ordinary MNIST modes 1–2 span is dominated by two atypical
observations; leverage then falls by several orders of magnitude.

![Normalized LTSA modes 1 and 2 for MNIST digit 1. The numbered images
cover the central plotted region and show how handwriting changes across
the map.](figures/ltsa-spectral-block-mnist-normalized-display.png)

Normalized LTSA modes 1 and 2 for MNIST digit 1. The numbered images
cover the central plotted region and show how handwriting changes across
the map.

Modes 2 and 3 give another view of the same morphology. Mode 3 has less
support, so this pair works best as a complementary diagnostic view.

![Normalized LTSA modes 2 and 3 expose a different sheet-and-branch view
of digit-1
morphology.](figures/ltsa-spectral-block-mnist-normalized-extra.png)

Normalized LTSA modes 2 and 3 expose a different sheet-and-branch view
of digit-1 morphology.

![MNIST low spectra and effective support fractions. Ordinary modes 1–3
approach the single-observation limit; normalized modes become
progressively more localized after the first
pair.](figures/ltsa-spectral-block-mnist-diagnostics.png)

MNIST low spectra and effective support fractions. Ordinary modes 1–3
approach the single-observation limit; normalized modes become
progressively more localized after the first pair.

Across all four retained modes, ordinary LTSA has block support
`0.000689` at the weak mode 4–5 gap `7.06e-5`. Normalized LTSA has block
support `0.0497` at the mode 4–5 gap `0.00196`; its modes 3 and 4 are
more localized than the displayed pair. These retained-block values
describe the exercised spans, not a claim that either boundary is
uniquely preferred.

The numbered representatives cover the central 98% of each plot axis.
Starting near the plot center, the selection repeatedly chooses the
image farthest from those already shown. This gives a compact tour of
the central region.

## Fashion-MNIST: localization depends on the estimator

Fashion-MNIST shows why a two-dimensional display should be read
together with its support diagnostics, and how estimator choice can
change both the map and the observations that dominate it. We fit
trousers, dresses, and bags separately with ordinary LTSA, using
`ndim = 2`, `k = 15`, and twelve retained modes. Each class uses one
fixed approximate-neighbor graph in memory.

Because dimensionality reduction is often judged through a
two-dimensional scatterplot, we first inspect the ordinary modes 1–2
maps. For each class, the span has block support about `0.000286`,
barely above the two-coordinate point-support benchmark
`2 / 7000 = 0.000285714`. The largest 1% of observations carry more than
99.98% of the displayed span’s leverage.

![The ordinary two-dimensional LTSA displays for three Fashion-MNIST
classes, shown at their full coordinate ranges. Orange points 1 and 2
are the two highest-leverage observations in each displayed span; the
matching numbered source images appear directly below each plot. Nearly
all other images collapse into the small gray clump near the
origin.](figures/ltsa-spectral-block-fashion-display.png)

The ordinary two-dimensional LTSA displays for three Fashion-MNIST
classes, shown at their full coordinate ranges. Orange points 1 and 2
are the two highest-leverage observations in each displayed span; the
matching numbered source images appear directly below each plot. Nearly
all other images collapse into the small gray clump near the origin.

These scatterplots faithfully show the returned coordinates, but they
are highly misleading as maps of typical within-class variation: two
observations set the visible scale while nearly 7,000 others are
compressed together. This diagnoses the ordinary-LTSA displayed prefix,
not Fashion-MNIST itself. The adjacent numbered images identify the
leverage maxima as atypical striped or noisy observations rather than
representative garments.

To isolate the effect of estimator weighting, we repeated each fit on
the same data and fixed graph, changing only to normalized LTSA. The
normalized modes 1–2 spans have block support `0.0898`, `0.0760`, and
`0.327` for trousers, dresses, and bags. Their top-1% leverage masses
fall to 27.2%, 32.3%, and 10.7%, with resolved mode 2–3 gaps of
`0.00186`, `0.000427`, and `0.00159`.

![Normalized LTSA modes 1 and 2 for the same three Fashion-MNIST classes
and graphs. Orange points 1 and 2 are the two highest-leverage
observations in each displayed span; the matching numbered source images
appear directly below each plot. Unlike the ordinary displays, thousands
of observations contribute visibly to each
map.](figures/ltsa-spectral-block-fashion-normalized-display.png)

Normalized LTSA modes 1 and 2 for the same three Fashion-MNIST classes
and graphs. Orange points 1 and 2 are the two highest-leverage
observations in each displayed span; the matching numbered source images
appear directly below each plot. Unlike the ordinary displays, thousands
of observations contribute visibly to each map.

Normalization changes the high-leverage observations as well as the map:
none of these numbered images belongs to the corresponding ordinary
pair. They are recognizable class members, but they are only the two
leverage maxima, not representatives of the normalized display. The much
lower top-1% masses show that leverage is distributed across many more
observations.

For these fixed classes and settings, normalization produces
substantially broader and visually richer displays, especially for bags.
This is not a claim that normalization is always better for
visualization: it changes the estimator’s weighting, and its usefulness
depends on the data, graph, and neighborhood. The support and boundary
diagnostics justify the interpretation here; the estimator name alone
does not.

The ordinary localization is not confined to the displayed pair. Across
the three classes, all 36 coordinates in modes 1–12 lie near the
single-observation support limit: their effective support fractions
range from `0.000142917` to `0.000146491`, compared with
`1 / 7000 = 0.000142857`. The largest 1% of observations carry at least
99.88% of every mode’s energy. The graphs are connected, local patches
have full rank, and the eigensolves converge.

![Fashion-MNIST spectra through mode 13 and per-mode effective support
through mode 12. Every retained mode is localized; the separators mark
the displayed and inspected prefix
boundaries.](figures/ltsa-spectral-block-fashion-diagnostics.png)

Fashion-MNIST spectra through mode 13 and per-mode effective support
through mode 12. Every retained mode is localized; the separators mark
the displayed and inspected prefix boundaries.

The displayed-pair thumbnails show what sets the scale of a
two-dimensional map. To ask what drives leverage across the retained
span, we instead rank observations in each full 12-mode block. These are
extreme-leverage diagnostics, not representative samples; numbering
restarts in every row.

For ordinary LTSA, twelve-mode block support is `0.001720` for trousers,
`0.001716` for dresses, and `0.001719` for bags, close to the
twelve-coordinate point-support benchmark `12 / 7000 = 0.001714`. The
largest 1% of observations carry at least 99.93% of each block’s
leverage. The mode 12–13 gaps are weak—`2.48e-5`, `1.46e-7`, and
`8.99e-7`—so the ordinary conclusion stops at the exercised block.

The normalized twelve-mode supports are broader but not uniformly as
broad as the normalized displayed pairs: `0.1269`, `0.0293`, and
`0.0222`. Their top-1% leverage masses are 13.8%, 29.6%, and 29.9%, at
resolved mode 12–13 gaps of `0.000866`, `0.00453`, and `0.000596`.

![The top 12 leverage observations for each ordinary and normalized
twelve-mode Fashion-MNIST block. Ordinary rows are dominated by striped,
patterned, or noisy images; normalized rows contain recognizable
garments. Numbering restarts for every estimator and
class.](figures/ltsa-spectral-block-fashion-high-leverage.png)

The top 12 leverage observations for each ordinary and normalized
twelve-mode Fashion-MNIST block. Ordinary rows are dominated by striped,
patterned, or noisy images; normalized rows contain recognizable
garments. Numbering restarts for every estimator and class.

The paired rankings show how the estimators differ. Only one trouser
appears in both top twelves; the dress and bag sets have no overlap. The
two normalized displayed-pair highlights remain ranks 1 and 2 in the
full trouser block and ranks 3 and 4 for bags, but fall to ranks 59 and
62 for dresses. Later normalized modes therefore shift which dresses
drive the retained block. This also explains why the normalized
full-block support can be more concentrated than the visually rich modes
1–2 map.

The eigengap describes boundary stability; block support describes how
broadly the retained span is represented across observations. Here, the
basis-dependent per-mode supports and the basis-invariant block support
lead to the same conclusion.

## What to carry into practice

The circle and COIL-20 show the clearest use: a locally one-dimensional
closed path appears as a pair in the low-frequency block. MNIST shows
how normalized LTSA can produce a coherent complementary view, while
support reveals whether each coordinate describes many observations or
only a few. Fashion-MNIST shows that estimator choice can change both a
fixed displayed pair and the observations driving a larger retained
block. Normalization broadens the displayed pair here without making
every later mode broadly supported or establishing a universal estimator
preference.

Use `spectral_dim` when topology or a weak boundary gives you a reason
to look past the displayed prefix. Read the block together with
connectivity, rank, solver, gap, and support diagnostics. For details of
the returned fields, see [Diagnosing suspicious LTSA
results](https://jlmelville.github.io/flotsam/articles/numerical-diagnostics.md).

Retaining a block supplies candidate coordinates; choosing which
coordinates to display is a separate problem. [Independent
Eigencoordinate
Selection](https://proceedings.neurips.cc/paper_files/paper/2019/file/6a10bbd480e4c5573d8f3af73ae0454b-Paper.pdf)
is one published method for selecting smooth, locally independent
coordinates from a larger spectral set. `spectral_dim` does not perform
that selection; it exposes the candidate modes needed to investigate it.

## Reproduce the real-data examples

The complete script downloads the datasets with `snedata`, reuses fixed
in-memory neighbor matrices for estimator comparisons, fits the retained
blocks, checks the claims above, and draws the thumbnail and support
plots. Install its extra dependencies with
`pak::pak("jlmelville/snedata")` and `install.packages("rnndescent")`.
Set `FLOTSAM_SPECTRAL_BLOCK_FIGURE_DIR` to a directory to write all ten
real-data article figures with the exact filenames and dimensions
declared in the script; otherwise they draw on the current device.

From the package source tree, after installing the current package, this
command regenerates the complete figure set and stops if a semantic
invariant fails:

``` bash
FLOTSAM_SPECTRAL_BLOCK_FIGURE_DIR=vignettes/articles/figures \
  Rscript vignettes/articles/spectral-blocks-reproduction.R
```

Show the complete data and plotting script (downloads three image
datasets)

``` r

# Reproduce the real-data examples in the spectral-blocks article.
#
# This script downloads COIL-20, MNIST, and Fashion-MNIST with snedata.
# The image datasets are fitted when their sections run, so expect several
# minutes of computation and substantial memory use.

library(flotsam)

mnist_neighbor_seed <- 20260824L
mnist_solver_seed <- 20260824L
fashion_neighbor_seed <- 20260825L
fashion_solver_seed <- 20260825L
figure_directory <- Sys.getenv(
  "FLOTSAM_SPECTRAL_BLOCK_FIGURE_DIR",
  unset = ""
)

article_figures <- data.frame(
  filename = c(
    "ltsa-spectral-block-coil-object-4.png",
    "ltsa-spectral-block-coil-object-1.png",
    "ltsa-spectral-block-mnist-ordinary-high-leverage.png",
    "ltsa-spectral-block-mnist-normalized-display.png",
    "ltsa-spectral-block-mnist-normalized-extra.png",
    "ltsa-spectral-block-mnist-diagnostics.png",
    "ltsa-spectral-block-fashion-display.png",
    "ltsa-spectral-block-fashion-normalized-display.png",
    "ltsa-spectral-block-fashion-diagnostics.png",
    "ltsa-spectral-block-fashion-high-leverage.png"
  ),
  width = c(
    1800L,
    1800L,
    2000L,
    1800L,
    1800L,
    1800L,
    2200L,
    2200L,
    2200L,
    2500L
  ),
  height = c(
    950L,
    950L,
    1100L,
    1000L,
    1000L,
    650L,
    1400L,
    1400L,
    1300L,
    1600L
  ),
  res = c(150L, 150L, 180L, 150L, 150L, 150L, 180L, 180L, 180L, 180L)
)

open_output_figure <- function(filename) {
  if (!nzchar(figure_directory)) {
    return(FALSE)
  }
  specification <- article_figures[article_figures$filename == filename, ]
  if (nrow(specification) != 1L) {
    stop("Unknown or duplicated article figure specification: ", filename)
  }
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(
    file.path(figure_directory, filename),
    width = specification$width,
    height = specification$height,
    res = specification$res
  )
  TRUE
}

close_output_figure <- function(opened) {
  if (opened) {
    grDevices::dev.off()
  }
  invisible(NULL)
}

drop_constant_columns <- function(X) {
  keep <- vapply(
    seq_len(ncol(X)),
    function(column) {
      values <- X[, column]
      min(values) != max(values)
    },
    logical(1)
  )
  X[, keep, drop = FALSE]
}

# Build one self-first approximate-neighbor matrix so several fits can reuse
# exactly the same graph. In flotsam, n_neighbors includes the self column.
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
  result
}

effective_support_fraction <- function(vectors) {
  vapply(
    seq_len(ncol(vectors)),
    function(column) {
      vector <- vectors[, column]
      unit <- vector / sqrt(sum(vector^2))
      1 / (length(unit) * sum(unit^4))
    },
    numeric(1)
  )
}

block_support_diagnostics <- function(vectors) {
  vectors <- as.matrix(vectors)
  column_norms <- sqrt(colSums(vectors^2))
  if (any(!is.finite(column_norms)) || any(column_norms == 0)) {
    stop("Coordinate blocks must have positive finite column norms")
  }

  scaled <- sweep(vectors, 2L, column_norms, "/")
  decomposition <- qr(scaled, tol = 1e-10, LAPACK = FALSE)
  if (decomposition$rank != ncol(scaled)) {
    stop("Coordinate block is numerically rank deficient")
  }

  basis <- qr.Q(decomposition, complete = FALSE)[,
    seq_len(ncol(scaled)),
    drop = FALSE
  ]
  leverage <- rowSums(basis^2)
  mass <- leverage / ncol(scaled)
  ordering <- order(-mass, seq_along(mass), method = "radix")
  top_count <- ceiling(0.01 * nrow(scaled))
  list(
    support = 1 / (nrow(scaled) * sum(mass^2)),
    top_one_percent_mass = sum(mass[ordering[seq_len(top_count)]]),
    mass = mass,
    ordering = ordering
  )
}

top_one_percent_mode_mass <- function(vectors) {
  vectors <- as.matrix(vectors)
  vapply(
    seq_len(ncol(vectors)),
    function(column) {
      energy <- vectors[, column]^2
      energy <- energy / sum(energy)
      top_count <- ceiling(0.01 * length(energy))
      sum(sort(energy, decreasing = TRUE)[seq_len(top_count)])
    },
    numeric(1)
  )
}

assert_fit_integrity <- function(fit, ndim, spectral_dim) {
  stopifnot(
    identical(dim(fit$embedding), c(nrow(fit$embedding), ndim)),
    identical(
      dim(fit$spectral_embedding),
      c(nrow(fit$embedding), spectral_dim)
    ),
    identical(
      fit$embedding,
      fit$spectral_embedding[, seq_len(ndim), drop = FALSE]
    ),
    fit$eigen$status != "invalid",
    fit$eigen$spectral$status != "invalid",
    fit$eigen$rank >= ndim,
    fit$eigen$spectral$rank >= spectral_dim,
    fit$assembly$component_count == 1L,
    fit$assembly$rank_deficient_count == 0L,
    fit$assembly$min_local_rank >= ndim
  )
  if (
    isTRUE(fit$eigen$backend$convergence_known) &&
      !is.null(fit$eigen$backend$nconv)
  ) {
    stopifnot(fit$eigen$backend$nconv == fit$eigen$eig_k)
  }
  invisible(fit)
}

assert_selection <- function(selected, population_size, count) {
  stopifnot(
    length(selected) == count,
    !anyDuplicated(selected),
    all(selected >= 1L),
    all(selected <= population_size)
  )
  invisible(selected)
}

central_window <- function(coordinates) {
  limits <- apply(
    coordinates,
    2L,
    stats::quantile,
    probs = c(0.01, 0.99),
    names = FALSE
  )
  spans <- apply(limits, 2L, diff)
  spans[spans == 0] <- 1
  outside <- coordinates[, 1L] < limits[1L, 1L] |
    coordinates[, 1L] > limits[2L, 1L] |
    coordinates[, 2L] < limits[1L, 2L] |
    coordinates[, 2L] > limits[2L, 2L]
  list(
    xlim = limits[, 1L] + c(-0.04, 0.04) * spans[[1L]],
    ylim = limits[, 2L] + c(-0.04, 0.04) * spans[[2L]],
    outside_count = sum(outside)
  )
}

representative_indices <- function(coordinates, count) {
  scaled <- scale(coordinates)
  scaled[!is.finite(scaled)] <- 0
  limits <- apply(
    scaled,
    2L,
    stats::quantile,
    probs = c(0.01, 0.99),
    names = FALSE
  )
  eligible <- which(
    scaled[, 1L] >= limits[1L, 1L] &
      scaled[, 1L] <= limits[2L, 1L] &
      scaled[, 2L] >= limits[1L, 2L] &
      scaled[, 2L] <= limits[2L, 2L]
  )

  center <- apply(scaled[eligible, , drop = FALSE], 2L, stats::median)
  selected <- eligible[[which.min(
    rowSums(sweep(scaled[eligible, , drop = FALSE], 2L, center)^2)
  )]]
  nearest <- rowSums(
    sweep(scaled[eligible, , drop = FALSE], 2L, scaled[selected, ])^2
  )

  while (length(selected) < count) {
    nearest[eligible %in% selected] <- -Inf
    next_position <- which.max(nearest)
    next_index <- eligible[[next_position]]
    selected <- c(selected, next_index)
    next_distance <- rowSums(
      sweep(
        scaled[eligible, , drop = FALSE],
        2L,
        scaled[next_index, ]
      )^2
    )
    nearest <- pmin(nearest, next_distance)
  }
  selected
}

pixel_raster <- function(pixel_row, height, transpose = TRUE) {
  image <- matrix(as.numeric(pixel_row), nrow = height)
  if (transpose) {
    image <- t(image)
  }
  maximum <- max(image)
  if (maximum > 0) {
    image <- image / maximum
  }
  as.raster(matrix(
    grDevices::gray(1 - image),
    nrow = nrow(image),
    ncol = ncol(image)
  ))
}

draw_image_panel <- function(
  X,
  selected,
  height,
  columns,
  transpose = TRUE,
  number_cex = 1.35
) {
  rows <- ceiling(length(selected) / columns)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, columns), ylim = c(0, rows), asp = 1)
  for (index in seq_along(selected)) {
    column <- (index - 1L) %% columns
    row <- rows - 1L - (index - 1L) %/% columns
    graphics::rasterImage(
      pixel_raster(X[selected[[index]], ], height, transpose),
      column + 0.04,
      row + 0.04,
      column + 0.96,
      row + 0.96,
      interpolate = FALSE
    )
    graphics::text(
      column + 0.14,
      row + 0.86,
      labels = index,
      col = "#D55E00",
      font = 2,
      cex = number_cex
    )
  }
}

draw_ranked_thumbnail_row <- function(images, selected, title) {
  count <- length(selected)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, count), ylim = c(-0.26, 1), asp = 1)
  for (position in seq_along(selected)) {
    row <- selected[[position]]
    left <- position - 1L
    graphics::rasterImage(
      pixel_raster(images[row, ], 28L),
      left + 0.04,
      0.04,
      left + 0.96,
      0.96,
      interpolate = FALSE
    )
    graphics::text(
      left + 0.08,
      0.86,
      labels = position,
      col = "#D55E00",
      font = 2,
      cex = 0.8,
      adj = c(0, 0.5)
    )
  }
  graphics::title(main = title, line = 0.2)
}

draw_ranked_thumbnail_grid <- function(
  images,
  selected,
  columns = 4L
) {
  rows <- ceiling(length(selected) / columns)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, columns), ylim = c(0, rows), asp = 1)
  for (position in seq_along(selected)) {
    column <- (position - 1L) %% columns
    row <- rows - 1L - (position - 1L) %/% columns
    selected_row <- selected[[position]]
    graphics::rasterImage(
      pixel_raster(images[selected_row, ], 28L),
      column + 0.04,
      row + 0.04,
      column + 0.96,
      row + 0.96,
      interpolate = FALSE
    )
    graphics::text(
      column + 0.1,
      row + 0.84,
      labels = position,
      col = "#D55E00",
      font = 2,
      cex = 0.9
    )
  }
}

draw_embedding_panel <- function(
  coordinates,
  selected,
  main,
  xlab,
  ylab,
  text_cex = 2
) {
  window <- central_window(coordinates)
  graphics::plot(
    coordinates,
    pch = 16,
    cex = 0.38,
    col = grDevices::adjustcolor("grey25", alpha.f = 0.20),
    asp = 1,
    xlim = window$xlim,
    ylim = window$ylim,
    xlab = xlab,
    ylab = ylab,
    main = main
  )
  graphics::text(
    coordinates[selected, , drop = FALSE],
    labels = seq_along(selected),
    col = "#D55E00",
    font = 2,
    cex = text_cex
  )
}

point_segment_distance <- function(point, start, end) {
  direction <- end - start
  squared_length <- sum(direction^2)
  if (squared_length == 0) {
    return(sqrt(sum((point - start)^2)))
  }
  position <- sum((point - start) * direction) / squared_length
  position <- max(0, min(1, position))
  closest <- start + position * direction
  sqrt(sum((point - closest)^2))
}

annotation_penalty <- function(distance, threshold, weight) {
  weight * max(0, threshold - distance)^2
}

choose_annotation_geometry <- function(highlighted, plot_limits) {
  stopifnot(nrow(highlighted) == 2L)
  x_span <- diff(plot_limits[1:2])
  y_span <- diff(plot_limits[3:4])
  scaled_points <- cbind(
    (highlighted[, 1L] - plot_limits[[1L]]) / x_span,
    (highlighted[, 2L] - plot_limits[[3L]]) / y_span
  )
  directions <- rbind(
    c(1, 0),
    c(-1, 0),
    c(0, 1),
    c(0, -1),
    c(1, 1),
    c(1, -1),
    c(-1, 1),
    c(-1, -1)
  )
  directions <- directions / sqrt(rowSums(directions^2))
  choices <- expand.grid(
    first = seq_len(nrow(directions)),
    second = seq_len(nrow(directions))
  )
  best_score <- Inf
  best_labels <- NULL
  label_distance <- 0.09

  for (choice in seq_len(nrow(choices))) {
    selected_directions <- directions[
      c(choices$first[[choice]], choices$second[[choice]]),
      ,
      drop = FALSE
    ]
    labels <- scaled_points + label_distance * selected_directions
    if (any(labels < 0.045 | labels > 0.955)) {
      next
    }

    score <- annotation_penalty(
      sqrt(sum((labels[1L, ] - labels[2L, ])^2)),
      0.1,
      3000
    )
    for (index in seq_len(2L)) {
      other <- 3L - index
      score <- score +
        annotation_penalty(
          sqrt(sum((labels[index, ] - scaled_points[other, ])^2)),
          0.12,
          4000
        )
      score <- score +
        annotation_penalty(
          point_segment_distance(
            scaled_points[other, ],
            scaled_points[index, ],
            labels[index, ]
          ),
          0.05,
          5000
        )
      score <- score +
        annotation_penalty(
          point_segment_distance(
            labels[other, ],
            scaled_points[index, ],
            labels[index, ]
          ),
          0.05,
          3000
        )
    }
    score <- score + 0.01 * sum((labels - 0.5)^2)

    if (score < best_score) {
      best_score <- score
      best_labels <- labels
    }
  }

  if (is.null(best_labels)) {
    stop("Could not place highlighted-observation annotations")
  }
  directions <- best_labels - scaled_points
  directions <- directions / sqrt(rowSums(directions^2))
  segment_starts <- scaled_points + 0.022 * directions
  segment_ends <- best_labels - 0.026 * directions
  to_user_coordinates <- function(points) {
    cbind(
      plot_limits[[1L]] + points[, 1L] * x_span,
      plot_limits[[3L]] + points[, 2L] * y_span
    )
  }
  list(
    labels = to_user_coordinates(best_labels),
    segment_starts = to_user_coordinates(segment_starts),
    segment_ends = to_user_coordinates(segment_ends)
  )
}

draw_full_embedding_panel <- function(
  coordinates,
  selected,
  main,
  estimator
) {
  graphics::plot(
    coordinates,
    pch = 16,
    cex = 0.52,
    col = grDevices::adjustcolor("grey20", alpha.f = 0.34),
    asp = 1,
    xlab = paste(estimator, "LTSA mode 1"),
    ylab = paste(estimator, "LTSA mode 2"),
    main = main
  )
  highlighted <- coordinates[selected, , drop = FALSE]
  annotation <- choose_annotation_geometry(
    highlighted,
    graphics::par("usr")
  )
  graphics::segments(
    annotation$segment_starts[, 1L],
    annotation$segment_starts[, 2L],
    annotation$segment_ends[, 1L],
    annotation$segment_ends[, 2L],
    col = "#D55E00",
    lwd = 1.4
  )
  graphics::points(
    highlighted,
    pch = 21,
    cex = 3.8,
    col = "white",
    bg = NA,
    lwd = 4.5
  )
  graphics::points(
    highlighted,
    pch = 21,
    cex = 3,
    col = "#D55E00",
    bg = NA,
    lwd = 2.4
  )
  graphics::text(
    annotation$labels,
    labels = seq_along(selected),
    col = "white",
    font = 2,
    cex = 1.8
  )
  graphics::text(
    annotation$labels,
    labels = seq_along(selected),
    col = "#D55E00",
    font = 2,
    cex = 1.2
  )
}

draw_high_leverage_pair <- function(images, selected) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 2), ylim = c(0, 1), asp = 1)
  for (position in seq_along(selected)) {
    left <- position - 1L
    graphics::rasterImage(
      pixel_raster(images[selected[[position]], ], 28L),
      left + 0.06,
      0.06,
      left + 0.94,
      0.94,
      interpolate = FALSE
    )
    graphics::points(
      left + 0.17,
      0.83,
      pch = 21,
      cex = 3.4,
      col = "white",
      bg = "white"
    )
    graphics::text(
      left + 0.17,
      0.83,
      labels = position,
      col = "#D55E00",
      font = 2,
      cex = 1.35
    )
  }
  graphics::title(
    main = "Source images 1 and 2",
    line = 0.2,
    cex.main = 1.05
  )
}

plot_pair_with_images <- function(
  coordinates,
  X,
  height,
  count,
  columns,
  pair,
  main,
  image_title,
  transpose = TRUE,
  image_number_cex = 1.35
) {
  selected <- representative_indices(coordinates, count)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::layout(matrix(c(1L, 1L, 1L, 2L, 2L), nrow = 1L))
  graphics::par(mar = c(4.4, 4.5, 4.2, 1.2))
  window <- central_window(coordinates)
  draw_embedding_panel(
    coordinates,
    selected,
    paste0(
      main,
      "\nRepresentative coverage; central 98% per axis; ",
      window$outside_count,
      " observations outside frame"
    ),
    paste0("same-operator LTSA mode ", pair[[1L]]),
    paste0("same-operator LTSA mode ", pair[[2L]])
  )
  graphics::par(mar = c(1, 1, 4.2, 1))
  draw_image_panel(
    X,
    selected,
    height,
    columns,
    transpose,
    image_number_cex
  )
  graphics::title(main = image_title)
  invisible(selected)
}

draw_low_spectrum <- function(
  fit,
  main,
  boundaries = c(2.5, 4.5),
  boundary_labels = c("displayed boundary", "retained boundary")
) {
  values <- fit$eigen$ritz_values[1:5]
  graphics::plot(
    seq_along(values),
    values,
    type = "b",
    pch = 21,
    bg = c("#0072B2", "#0072B2", "grey65", "grey65", "white"),
    xaxt = "n",
    xlab = "nonconstant Ritz mode",
    ylab = "Ritz value",
    main = main
  )
  graphics::axis(1, at = seq_along(values))
  boundary_lty <- seq_along(boundaries) + 1L
  graphics::abline(
    v = boundaries,
    lty = boundary_lty,
    col = "grey45"
  )
  graphics::legend(
    "topleft",
    legend = boundary_labels,
    lty = boundary_lty,
    col = "grey45",
    bty = "n",
    cex = 0.84
  )
}

draw_extended_spectrum <- function(fit, main, retained_modes) {
  values <- fit$eigen$ritz_values[seq_len(retained_modes + 1L)]
  graphics::plot(
    seq_along(values),
    values,
    type = "b",
    pch = 21,
    bg = c(rep("#0072B2", retained_modes), "white"),
    col = "grey25",
    xaxt = "n",
    xlab = "nonconstant Ritz mode",
    ylab = "Ritz value",
    main = main
  )
  graphics::axis(1, at = seq_along(values))
  graphics::abline(
    v = c(2.5, 4.5, 8.5, retained_modes + 0.5),
    lty = c(2, 3, 3, 2),
    col = "grey45"
  )
}

plot_support <- function(
  fits,
  labels,
  colors,
  main,
  legend_position = "left",
  boundary_markers = numeric(),
  xlim = NULL
) {
  mode_count <- ncol(fits[[1L]]$spectral_embedding)
  support <- vapply(
    fits,
    function(fit) effective_support_fraction(fit$spectral_embedding),
    numeric(mode_count)
  )
  graphics::matplot(
    seq_len(mode_count),
    support,
    type = "b",
    log = "y",
    pch = seq.int(21L, length.out = length(fits)),
    col = colors,
    bg = colors,
    xlim = xlim,
    ylim = range(c(support, 0.01)),
    xaxt = "n",
    xlab = "nonconstant mode",
    ylab = "effective support fraction (log scale)",
    main = main
  )
  graphics::axis(1, at = seq_len(mode_count))
  if (length(boundary_markers) > 0L) {
    graphics::abline(v = boundary_markers, lty = 3, col = "grey70")
  }
  graphics::abline(h = 0.01, lty = 3, col = "grey45")
  graphics::legend(
    legend_position,
    legend = c(labels, "1% support"),
    col = c(colors, "grey45"),
    pch = c(seq.int(21L, length.out = length(fits)), NA),
    pt.bg = c(colors, NA),
    lty = c(rep(1, length(fits)), 3),
    bty = "n"
  )
}

# COIL-20 ---------------------------------------------------------------

coil <- snedata::download_coil20(as = "list")
coil_objects <- c(4L, 1L)
coil_fits <- lapply(coil_objects, function(object) {
  rows <- coil$meta$object == object
  X <- coil$data[rows, , drop = FALSE]
  poses <- as.integer(coil$meta$pose[rows])
  fit <- ltsa(
    X,
    ndim = 1,
    n_neighbors = 7,
    nn_method = "exact",
    eig_method = "eig",
    eig_k = 8,
    output = "result",
    spectral_dim = 4
  )
  assert_fit_integrity(fit, ndim = 1L, spectral_dim = 4L)
  pair_gaps <- diff(fit$eigen$ritz_values[seq_len(3L)])
  stopifnot(pair_gaps[[2L]] > 10 * pair_gaps[[1L]])
  list(object = object, X = X, poses = poses, fit = fit)
})

plot_coil_orbit <- function(fit, X, poses, object) {
  coordinates <- fit$spectral_embedding[, 1:2, drop = FALSE]
  x_range <- range(coordinates[, 1L])
  y_range <- range(coordinates[, 2L])
  x_span <- diff(range(coordinates[, 1L]))
  y_span <- diff(range(coordinates[, 2L]))
  xlim <- x_range + c(-0.46, 0.46) * x_span
  ylim <- y_range + c(-0.46, 0.46) * y_span
  center <- colMeans(coordinates)
  selected <- match(seq.int(0L, 63L, by = 9L), poses)
  direction <- cbind(
    (coordinates[selected, 1L] - center[[1L]]) / x_span,
    (coordinates[selected, 2L] - center[[2L]]) / y_span
  )
  angle <- atan2(direction[, 2L], direction[, 1L])
  plot_radius <- c(0.40 * diff(xlim), 0.40 * diff(ylim))
  centers <- cbind(
    center[[1L]] + plot_radius[[1L]] * cos(angle),
    center[[2L]] + plot_radius[[2L]] * sin(angle)
  )

  colors <- grDevices::hcl.colors(72L, "Spectral", rev = TRUE)
  graphics::plot(
    coordinates,
    type = "n",
    asp = 1,
    xlim = xlim,
    ylim = ylim,
    xlab = "LTSA mode 1",
    ylab = "LTSA mode 2",
    main = paste0("COIL-20 object ", object, ": rotational orbit")
  )
  graphics::points(
    coordinates,
    pch = 21,
    bg = colors[poses + 1L],
    col = "grey20",
    cex = 1.25,
    lwd = 0.45
  )
  for (index in seq_along(selected)) {
    row <- selected[[index]]
    graphics::segments(
      coordinates[row, 1L],
      coordinates[row, 2L],
      centers[index, 1L],
      centers[index, 2L],
      col = grDevices::adjustcolor("grey30", alpha.f = 0.7),
      lwd = 0.8
    )
  }
  for (index in seq_along(selected)) {
    row <- selected[[index]]
    width <- 0.19 * x_span
    height <- 0.19 * y_span
    graphics::rect(
      centers[index, 1L] - width / 2,
      centers[index, 2L] - height / 2,
      centers[index, 1L] + width / 2,
      centers[index, 2L] + height / 2,
      col = "white",
      border = "grey35",
      lwd = 0.8
    )
    graphics::rasterImage(
      matrix(as.numeric(X[row, ]), nrow = 128L, ncol = 128L),
      centers[index, 1L] - width / 2,
      centers[index, 2L] - height / 2,
      centers[index, 1L] + width / 2,
      centers[index, 2L] + height / 2,
      interpolate = FALSE
    )
    graphics::text(
      centers[index, 1L],
      centers[index, 2L] - 0.62 * height,
      labels = paste("pose", poses[[row]]),
      cex = 0.9,
      xpd = NA
    )
  }
}

draw_coil_low_spectrum <- function(fit, object) {
  values <- fit$eigen$ritz_values[1:5]
  graphics::plot(
    seq_along(values),
    values,
    type = "b",
    pch = 21,
    bg = c("#D55E00", "#D55E00", "grey65", "grey65", "white"),
    col = "grey25",
    xaxt = "n",
    xlab = "nonconstant Ritz mode",
    ylab = "Ritz value",
    main = paste0("Object ", object, ": low spectrum")
  )
  graphics::axis(1, at = seq_along(values))
  graphics::abline(
    v = c(1.5, 2.5, 4.5),
    lty = c(3, 2, 3),
    col = "grey45"
  )
  graphics::legend(
    "topleft",
    legend = c(
      "displayed 1D boundary",
      "periodic-pair boundary",
      "retained block boundary"
    ),
    lty = c(3, 2, 3),
    col = "grey45",
    bty = "n",
    cex = 0.86
  )
}

for (coil_result in coil_fits) {
  coil_figure_open <- open_output_figure(paste0(
    "ltsa-spectral-block-coil-object-",
    coil_result$object,
    ".png"
  ))
  old_par <- graphics::par(no.readonly = TRUE)
  graphics::layout(matrix(c(1L, 1L, 2L), nrow = 1L))
  graphics::par(mar = c(4.4, 4.6, 3.2, 1.2))
  plot_coil_orbit(
    coil_result$fit,
    coil_result$X,
    coil_result$poses,
    coil_result$object
  )
  graphics::par(mar = c(4.4, 4.5, 3.2, 1.2))
  draw_coil_low_spectrum(
    coil_result$fit,
    coil_result$object
  )
  graphics::par(old_par)
  close_output_figure(coil_figure_open)
}

# MNIST digit 1 ---------------------------------------------------------

mnist <- snedata::download_mnist(as = "list")
mnist_labels <- as.integer(as.character(mnist$meta$label))
mnist_images <- mnist$data[mnist_labels == 1L, , drop = FALSE]
mnist_X <- drop_constant_columns(mnist_images)
mnist_graph <- self_first_nnd(mnist_X, 15L, mnist_neighbor_seed)

fit_mnist <- function(normalize) {
  set.seed(mnist_solver_seed)
  ltsa(
    mnist_X,
    ndim = 2,
    n_neighbors = 15,
    nn_method = mnist_graph,
    eig_method = "rspectra",
    eig_k = 8,
    normalize = normalize,
    output = "result",
    spectral_dim = 4
  )
}

mnist_ordinary <- fit_mnist(FALSE)
mnist_normalized <- fit_mnist(TRUE)
assert_fit_integrity(mnist_ordinary, ndim = 2L, spectral_dim = 4L)
assert_fit_integrity(mnist_normalized, ndim = 2L, spectral_dim = 4L)
mnist_ordinary_display <- block_support_diagnostics(
  mnist_ordinary$spectral_embedding[, 1:2, drop = FALSE]
)
mnist_ordinary_retained <- block_support_diagnostics(
  mnist_ordinary$spectral_embedding
)
mnist_normalized_display <- block_support_diagnostics(
  mnist_normalized$spectral_embedding[, 1:2, drop = FALSE]
)
mnist_normalized_retained <- block_support_diagnostics(
  mnist_normalized$spectral_embedding
)

stopifnot(
  mnist_ordinary_display$support < 0.001,
  mnist_normalized_display$support > 0.1,
  mnist_normalized_display$support > 100 * mnist_ordinary_display$support,
  mnist_ordinary_display$top_one_percent_mass > 0.99,
  mnist_normalized_display$top_one_percent_mass < 0.1,
  mnist_ordinary_retained$support < 0.002,
  mnist_normalized_retained$support > 0.01,
  mnist_normalized_retained$support > mnist_ordinary_retained$support,
  mnist_ordinary$eigen$diagnostics$scaled_boundary_gap < 1e-4,
  mnist_normalized$eigen$diagnostics$scaled_boundary_gap > 1e-4,
  mnist_ordinary$eigen$spectral$diagnostics$scaled_boundary_gap < 1e-4,
  mnist_normalized$eigen$spectral$diagnostics$scaled_boundary_gap > 1e-4
)

print(data.frame(
  estimator = rep(c("ordinary", "normalized"), each = 2L),
  prefix = rep(c(2L, 4L), 2L),
  support = c(
    mnist_ordinary_display$support,
    mnist_ordinary_retained$support,
    mnist_normalized_display$support,
    mnist_normalized_retained$support
  ),
  top_one_percent_mass = c(
    mnist_ordinary_display$top_one_percent_mass,
    mnist_ordinary_retained$top_one_percent_mass,
    mnist_normalized_display$top_one_percent_mass,
    mnist_normalized_retained$top_one_percent_mass
  ),
  boundary_gap = c(
    mnist_ordinary$eigen$diagnostics$scaled_boundary_gap,
    mnist_ordinary$eigen$spectral$diagnostics$scaled_boundary_gap,
    mnist_normalized$eigen$diagnostics$scaled_boundary_gap,
    mnist_normalized$eigen$spectral$diagnostics$scaled_boundary_gap
  )
))

mnist_selected <- mnist_ordinary_display$ordering[seq_len(16L)]
assert_selection(mnist_selected, nrow(mnist_images), 16L)
mnist_figure_open <- open_output_figure(
  "ltsa-spectral-block-mnist-ordinary-high-leverage.png"
)
old_par <- graphics::par(no.readonly = TRUE)
graphics::layout(matrix(c(1L, 2L, 2L), nrow = 1L))
graphics::par(mar = c(4.2, 4.4, 3.3, 1.0))
graphics::plot(
  seq_along(mnist_ordinary_display$mass),
  mnist_ordinary_display$mass[mnist_ordinary_display$ordering],
  type = "l",
  log = "y",
  xlab = "observation leverage rank",
  ylab = "normalized block leverage mass",
  main = "Ordinary MNIST modes 1--2"
)
graphics::abline(v = 16.5, lty = 2, col = "#D55E00")
graphics::par(mar = c(0.5, 0.5, 3.3, 0.5))
draw_ranked_thumbnail_grid(
  mnist_images,
  mnist_selected,
  columns = 4L
)
graphics::title(main = "Top 16 observations by displayed-block leverage")
graphics::par(old_par)
close_output_figure(mnist_figure_open)

mnist_normalized_display_open <- open_output_figure(
  "ltsa-spectral-block-mnist-normalized-display.png"
)
mnist_normalized_display_selected <- plot_pair_with_images(
  mnist_normalized$spectral_embedding[, 1:2, drop = FALSE],
  mnist_images,
  height = 28L,
  count = 16L,
  columns = 4L,
  pair = 1:2,
  main = "Normalized LTSA digit 1: modes 1 and 2",
  image_title = "Representative digits from the map"
)
close_output_figure(mnist_normalized_display_open)
assert_selection(
  mnist_normalized_display_selected,
  nrow(mnist_images),
  16L
)

mnist_normalized_extra_open <- open_output_figure(
  "ltsa-spectral-block-mnist-normalized-extra.png"
)
mnist_normalized_extra_selected <- plot_pair_with_images(
  mnist_normalized$spectral_embedding[, 2:3, drop = FALSE],
  mnist_images,
  height = 28L,
  count = 16L,
  columns = 4L,
  pair = 2:3,
  main = "Normalized LTSA digit 1: modes 2 and 3",
  image_title = "Representative digits from the map"
)
close_output_figure(mnist_normalized_extra_open)
assert_selection(mnist_normalized_extra_selected, nrow(mnist_images), 16L)

mnist_diagnostics_open <- open_output_figure(
  "ltsa-spectral-block-mnist-diagnostics.png"
)
old_par <- graphics::par(no.readonly = TRUE)
graphics::par(mfrow = c(1L, 3L), mar = c(4.2, 4.5, 3.3, 1))
draw_low_spectrum(
  mnist_ordinary,
  "Ordinary LTSA low spectrum",
  boundary_labels = c("displayed 2D boundary", "retained 4D boundary")
)
draw_low_spectrum(
  mnist_normalized,
  "Normalized LTSA low spectrum",
  boundary_labels = c("displayed 2D boundary", "retained 4D boundary")
)
mnist_ordinary_support <- effective_support_fraction(
  mnist_ordinary$spectral_embedding
)
mnist_normalized_support <- effective_support_fraction(
  mnist_normalized$spectral_embedding
)
graphics::plot(
  seq_len(4L),
  mnist_ordinary_support,
  type = "b",
  log = "y",
  pch = 21,
  bg = "#D55E00",
  col = "#D55E00",
  ylim = range(c(mnist_ordinary_support, mnist_normalized_support)),
  xlim = c(1, 5.7),
  xaxt = "n",
  xlab = "nonconstant mode",
  ylab = "effective support fraction (log scale)",
  main = "Low modes need not be global"
)
graphics::axis(1, at = seq_len(4L))
graphics::lines(
  seq_len(4L),
  mnist_normalized_support,
  type = "b",
  pch = 21,
  bg = "#0072B2",
  col = "#0072B2"
)
graphics::abline(h = 0.01, lty = 3, col = "grey45")
graphics::legend(
  "topright",
  legend = c("ordinary", "normalized", "1% support"),
  col = c("#D55E00", "#0072B2", "grey45"),
  pch = c(21, 21, NA),
  pt.bg = c("#D55E00", "#0072B2", NA),
  lty = c(1, 1, 3),
  bty = "n"
)
graphics::par(old_par)
close_output_figure(mnist_diagnostics_open)

mnist_selection_scopes <- data.frame(
  panel = c(
    "ordinary-high-leverage",
    "normalized-display",
    "normalized-extra"
  ),
  estimator = c("ordinary", "normalized", "normalized"),
  modes = c("1:2", "1:2", "2:3"),
  rule = c(
    "top displayed-block leverage",
    "central-region representative coverage",
    "central-region representative coverage"
  ),
  count = c(
    length(mnist_selected),
    length(mnist_normalized_display_selected),
    length(mnist_normalized_extra_selected)
  )
)
print(mnist_selection_scopes, row.names = FALSE)

# Fashion-MNIST ---------------------------------------------------------

fashion <- snedata::download_fashion_mnist(as = "list")
fashion_labels <- as.integer(as.character(fashion$meta$label))
fashion_classes <- c(Trouser = 1L, Dress = 3L, Bag = 8L)

fashion_fits <- lapply(seq_along(fashion_classes), function(index) {
  label <- fashion_classes[[index]]
  rows <- fashion_labels == label
  images <- fashion$data[rows, , drop = FALSE]
  X <- drop_constant_columns(images)
  graph <- self_first_nnd(X, 15L, fashion_neighbor_seed)
  fit_fashion <- function(normalize) {
    set.seed(fashion_solver_seed)
    ltsa(
      X,
      ndim = 2,
      n_neighbors = 15,
      nn_method = graph,
      eig_method = "rspectra",
      eig_k = 16,
      normalize = normalize,
      output = "result",
      spectral_dim = 12
    )
  }
  fit <- fit_fashion(FALSE)
  normalized_fit <- fit_fashion(TRUE)
  assert_fit_integrity(fit, ndim = 2L, spectral_dim = 12L)
  assert_fit_integrity(normalized_fit, ndim = 2L, spectral_dim = 12L)
  list(
    fit = fit,
    normalized_fit = normalized_fit,
    images = images,
    display_block = block_support_diagnostics(
      fit$spectral_embedding[, 1:2, drop = FALSE]
    ),
    normalized_display_block = block_support_diagnostics(
      normalized_fit$spectral_embedding[, 1:2, drop = FALSE]
    ),
    block = block_support_diagnostics(fit$spectral_embedding),
    normalized_block = block_support_diagnostics(
      normalized_fit$spectral_embedding
    )
  )
})
names(fashion_fits) <- names(fashion_classes)

fashion_ordinary_mode_support <- lapply(
  fashion_fits,
  function(result) effective_support_fraction(result$fit$spectral_embedding)
)
fashion_ordinary_mode_top_one_percent <- lapply(
  fashion_fits,
  function(result) {
    top_one_percent_mode_mass(result$fit$spectral_embedding)
  }
)
ordinary_mode_support <- unlist(
  fashion_ordinary_mode_support,
  use.names = FALSE
)
ordinary_mode_top_one_percent <- unlist(
  fashion_ordinary_mode_top_one_percent,
  use.names = FALSE
)
stopifnot(
  min(ordinary_mode_support) >= 1 / nrow(fashion_fits[[1L]]$images) - 1e-12,
  max(ordinary_mode_support) < 0.00015,
  min(ordinary_mode_top_one_percent) > 0.998
)

fashion_display_summary <- do.call(
  rbind,
  lapply(names(fashion_fits), function(class_name) {
    result <- fashion_fits[[class_name]]
    data.frame(
      class = class_name,
      estimator = c("ordinary", "normalized"),
      prefix = 2L,
      support = c(
        result$display_block$support,
        result$normalized_display_block$support
      ),
      top_one_percent_mass = c(
        result$display_block$top_one_percent_mass,
        result$normalized_display_block$top_one_percent_mass
      ),
      boundary_gap = c(
        result$fit$eigen$diagnostics$scaled_boundary_gap,
        result$normalized_fit$eigen$diagnostics$scaled_boundary_gap
      )
    )
  })
)
rownames(fashion_display_summary) <- NULL
print(fashion_display_summary)
ordinary_display_rows <- fashion_display_summary$estimator == "ordinary"
normalized_display_rows <- !ordinary_display_rows
stopifnot(
  all(fashion_display_summary$support[ordinary_display_rows] < 0.001),
  all(fashion_display_summary$support[normalized_display_rows] > 0.05),
  all(
    fashion_display_summary$support[normalized_display_rows] >
      fashion_display_summary$support[ordinary_display_rows]
  ),
  all(
    fashion_display_summary$top_one_percent_mass[ordinary_display_rows] > 0.999
  ),
  all(
    fashion_display_summary$top_one_percent_mass[normalized_display_rows] < 0.4
  ),
  all(
    fashion_display_summary$top_one_percent_mass[normalized_display_rows] <
      fashion_display_summary$top_one_percent_mass[ordinary_display_rows]
  ),
  all(fashion_display_summary$boundary_gap[ordinary_display_rows] < 1e-4),
  all(fashion_display_summary$boundary_gap[normalized_display_rows] > 1e-4)
)

fashion_display_open <- open_output_figure(
  "ltsa-spectral-block-fashion-display.png"
)
old_par <- graphics::par(no.readonly = TRUE)
graphics::layout(
  matrix(seq_len(6L), nrow = 2L, byrow = TRUE),
  heights = c(3, 1.2)
)
graphics::par(mar = c(4.2, 4.4, 3.2, 1.0))
fashion_display_selected <- list()
for (class_name in names(fashion_fits)) {
  result <- fashion_fits[[class_name]]
  coordinates <- result$fit$embedding
  selected <- result$display_block$ordering[seq_len(2L)]
  assert_selection(selected, nrow(result$images), 2L)
  fashion_display_selected[[class_name]] <- selected
  draw_full_embedding_panel(
    coordinates,
    selected,
    class_name,
    "ordinary"
  )
}
graphics::par(mar = c(0.8, 0.4, 2.4, 0.4))
for (class_name in names(fashion_fits)) {
  result <- fashion_fits[[class_name]]
  selected <- result$display_block$ordering[seq_len(2L)]
  draw_high_leverage_pair(result$images, selected)
}
graphics::layout(matrix(1L))
graphics::par(old_par)
close_output_figure(fashion_display_open)

fashion_normalized_display_open <- open_output_figure(
  "ltsa-spectral-block-fashion-normalized-display.png"
)
old_par <- graphics::par(no.readonly = TRUE)
graphics::layout(
  matrix(seq_len(6L), nrow = 2L, byrow = TRUE),
  heights = c(3, 1.2)
)
graphics::par(mar = c(4.2, 4.4, 3.2, 1.0))
fashion_normalized_display_selected <- list()
for (class_name in names(fashion_fits)) {
  result <- fashion_fits[[class_name]]
  coordinates <- result$normalized_fit$embedding
  selected <- result$normalized_display_block$ordering[seq_len(2L)]
  assert_selection(selected, nrow(result$images), 2L)
  fashion_normalized_display_selected[[class_name]] <- selected
  draw_full_embedding_panel(
    coordinates,
    selected,
    class_name,
    "normalized"
  )
}
graphics::par(mar = c(0.8, 0.4, 2.4, 0.4))
for (class_name in names(fashion_fits)) {
  result <- fashion_fits[[class_name]]
  selected <- result$normalized_display_block$ordering[seq_len(2L)]
  draw_high_leverage_pair(result$images, selected)
}
graphics::layout(matrix(1L))
graphics::par(old_par)
close_output_figure(fashion_normalized_display_open)

fashion_spectrum_open <- open_output_figure(
  "ltsa-spectral-block-fashion-diagnostics.png"
)
old_par <- graphics::par(no.readonly = TRUE)
graphics::par(mfrow = c(2L, 2L), mar = c(4.2, 4.5, 3.2, 1.0))
for (class_name in names(fashion_fits)) {
  draw_extended_spectrum(
    fashion_fits[[class_name]]$fit,
    paste(class_name, "modes 1--13"),
    retained_modes = 12L
  )
}
plot_support(
  lapply(fashion_fits, function(x) x$fit),
  names(fashion_fits),
  c("#0072B2", "#D55E00", "#009E73"),
  "Per-mode support through mode 12",
  boundary_markers = c(2.5, 4.5, 8.5, 12.5),
  xlim = c(1, 13)
)
graphics::par(old_par)
close_output_figure(fashion_spectrum_open)

fashion_block_summary <- do.call(
  rbind,
  lapply(names(fashion_fits), function(class_name) {
    result <- fashion_fits[[class_name]]
    data.frame(
      class = class_name,
      estimator = c("ordinary", "normalized"),
      prefix = 12L,
      support = c(result$block$support, result$normalized_block$support),
      top_one_percent_mass = c(
        result$block$top_one_percent_mass,
        result$normalized_block$top_one_percent_mass
      ),
      boundary_gap = c(
        result$fit$eigen$spectral$diagnostics$scaled_boundary_gap,
        result$normalized_fit$eigen$spectral$diagnostics$scaled_boundary_gap
      )
    )
  })
)
rownames(fashion_block_summary) <- NULL
print(fashion_block_summary)
ordinary_block_rows <- fashion_block_summary$estimator == "ordinary"
normalized_block_rows <- !ordinary_block_rows
stopifnot(
  all(fashion_block_summary$support[ordinary_block_rows] > 0.0015),
  all(fashion_block_summary$support[ordinary_block_rows] < 0.002),
  all(fashion_block_summary$support[normalized_block_rows] > 0.02),
  all(
    fashion_block_summary$support[normalized_block_rows] >
      10 * fashion_block_summary$support[ordinary_block_rows]
  ),
  all(
    fashion_block_summary$top_one_percent_mass[ordinary_block_rows] > 0.999
  ),
  all(
    fashion_block_summary$top_one_percent_mass[normalized_block_rows] < 0.31
  ),
  all(fashion_block_summary$boundary_gap[ordinary_block_rows] < 1e-4),
  all(fashion_block_summary$boundary_gap[normalized_block_rows] > 1e-4)
)

old_par <- graphics::par(no.readonly = TRUE)
fashion_thumbnail_open <- open_output_figure(
  "ltsa-spectral-block-fashion-high-leverage.png"
)
graphics::par(mfrow = c(6L, 1L), mar = c(2.4, 0.5, 2.1, 0.5))
fashion_block_selected <- list()
for (class_name in names(fashion_fits)) {
  result <- fashion_fits[[class_name]]
  blocks <- list(
    ordinary = result$block,
    normalized = result$normalized_block
  )
  for (estimator in names(blocks)) {
    selected <- blocks[[estimator]]$ordering[seq_len(12L)]
    assert_selection(selected, nrow(result$images), 12L)
    fashion_block_selected[[paste(class_name, estimator, sep = "-")]] <-
      selected
    draw_ranked_thumbnail_row(
      result$images,
      selected,
      paste(class_name, estimator, "top 12 from modes 1--12")
    )
  }
}
graphics::par(old_par)
close_output_figure(fashion_thumbnail_open)

display_overlap <- vapply(
  names(fashion_fits),
  function(class_name) {
    length(intersect(
      fashion_display_selected[[class_name]],
      fashion_normalized_display_selected[[class_name]]
    ))
  },
  integer(1)
)
block_overlap <- vapply(
  names(fashion_fits),
  function(class_name) {
    length(intersect(
      fashion_block_selected[[paste(class_name, "ordinary", sep = "-")]],
      fashion_block_selected[[paste(class_name, "normalized", sep = "-")]]
    ))
  },
  integer(1)
)
normalized_display_ranks <- lapply(
  names(fashion_fits),
  function(class_name) {
    result <- fashion_fits[[class_name]]
    sort(match(
      fashion_normalized_display_selected[[class_name]],
      result$normalized_block$ordering
    ))
  }
)
names(normalized_display_ranks) <- names(fashion_fits)

stopifnot(
  identical(unname(display_overlap), c(0L, 0L, 0L)),
  identical(unname(block_overlap), c(1L, 0L, 0L)),
  identical(normalized_display_ranks$Trouser, c(1L, 2L)),
  identical(normalized_display_ranks$Dress, c(59L, 62L)),
  identical(normalized_display_ranks$Bag, c(3L, 4L))
)

fashion_selection_scopes <- data.frame(
  class = rep(names(fashion_fits), each = 4L),
  estimator = rep(
    c("ordinary", "normalized", "ordinary", "normalized"),
    times = 3L
  ),
  modes = rep(c("1:2", "1:2", "1:12", "1:12"), times = 3L),
  rule = rep(
    c(
      "top displayed-block leverage",
      "top displayed-block leverage",
      "top retained-block leverage",
      "top retained-block leverage"
    ),
    times = 3L
  )
)
print(fashion_selection_scopes, row.names = FALSE)

if (nzchar(figure_directory)) {
  figure_paths <- file.path(figure_directory, article_figures$filename)
  stopifnot(
    all(file.exists(figure_paths)),
    all(file.info(figure_paths)$size > 0)
  )
  message(
    "Wrote all ",
    length(figure_paths),
    " spectral-block article figures."
  )
}
```

## Data sources

The image examples use
[COIL-20](https://cave.cs.columbia.edu/repository/COIL-20) (Sameer A.
Nene, Shree K. Nayar, and Hiroshi Murase),
[MNIST](http://yann.lecun.com/exdb/mnist/) (Yann LeCun, Corinna Cortes,
and Christopher J. C. Burges), and
[Fashion-MNIST](https://arxiv.org/abs/1708.07747) (Han Xiao, Kashif
Rasul, and Roland Vollgraf).
