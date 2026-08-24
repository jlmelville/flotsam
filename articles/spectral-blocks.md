# When the first LTSA coordinates are not the whole story

LTSA usually returns as many coordinates as the local tangent dimension
requested with `ndim`. That is a useful default, but it can hide a small
ambiguity: the first `ndim` eigenvectors are an ordered slice of a
low-frequency spectral block, and the boundary of that slice may be
weak. A nearby mode can contain a complementary view of the same fixed
LTSA operator.

This article shows when retaining a few extra modes is useful, and just
as importantly, when it is not. The examples are pedagogical subsets
rather than a claim that LTSA embeds an entire image collection well in
one shot.

## One operator, two dimensions to choose

`ndim` has two jobs in an ordinary call: it sets the local tangent rank
used to construct the alignment operator `B`, and it sets the number of
coordinates in the returned embedding. Changing `ndim` therefore
rebuilds `B`.

`spectral_dim` does something narrower. With `output = "result"`, it
retains more nonconstant eigenvectors after constructing the same
operator:

``` text
data + neighborhoods + ndim  --->  fixed LTSA operator B
                                         |
                                         v
                              low nonconstant modes
                               |               |
                               v               v
                    embedding[, 1:ndim]   spectral_embedding
```

`eig_k` is different again: it is the wider numerical candidate request
passed to the eigensolver and Rayleigh–Ritz postprocessing. It does not
set the local tangent rank or the number of modes returned. When
`spectral_dim > ndim`, the request must leave room for the known
constant null direction and one mode beyond the retained block, so
`eig_k` must be at least `spectral_dim + 2`. That extra boundary mode is
what makes the retained-block gap diagnosable.

For example:

``` r

fit <- ltsa(
  X,
  ndim = 2,
  n_neighbors = 15,
  output = "result",
  spectral_dim = 4
)

dim(fit$embedding)                   # n by 2
dim(fit$spectral_embedding)          # n by 4
identical(
  fit$embedding,
  fit$spectral_embedding[, 1:2, drop = FALSE]
)
```

The displayed and retained boundaries have separate diagnostics:

``` r

fit$eigen$diagnostics$scaled_boundary_gap
fit$eigen$spectral$diagnostics$scaled_boundary_gap
fit$eigen$status
fit$eigen$spectral$status
```

The first gap describes the displayed `ndim`-dimensional prefix. The
second describes the boundary after the retained `spectral_dim` modes. A
weak displayed gap says that the displayed subspace is poorly separated
from the next direction. It does **not** say that the next direction is
useful.

## A circle needs a spectral pair

A circle is locally one-dimensional, so `ndim = 1` is the appropriate
local tangent rank. Globally, however, no single continuous real-valued
coordinate can go once around a closed loop without identifying or
cutting points. Two low-frequency coordinates can represent the periodic
progression together.

The example below uses exact neighbors and dense eigenanalysis so the
small example is deterministic and the repeated spectral pairs are
recovered completely.

``` r

theta <- seq(0, 2 * pi, length.out = 721)[-721]
circle <- cbind(x = cos(theta), y = sin(theta))

circle_fit <- suppressWarnings(ltsa(
  circle,
  ndim = 1,
  n_neighbors = 9,
  nn_method = "exact",
  eig_method = "eig",
  eig_k = 8,
  output = "result",
  spectral_dim = 4
))

circle_fit$eigen$values
#> [1] 4.464673e-07
circle_fit$eigen$spectral$values
#> [1] 4.464673e-07 4.464673e-07 7.139639e-06 7.139639e-06
circle_fit$eigen$status
#> [1] "warning"
circle_fit$eigen$messages
#> [1] "Weak Ritz boundary gap after the selected block: 1.293e-17 < 1e-04."
```

![A locally one-dimensional circle. One returned coordinate cannot show
the closed progression, while the first same-operator spectral pair
does. The low spectrum contains repeated
pairs.](spectral-blocks_files/figure-html/circle-figure-1.png)

A locally one-dimensional circle. One returned coordinate cannot show
the closed progression, while the first same-operator spectral pair
does. The low spectrum contains repeated pairs.

The orientation and signs of the pair are arbitrary; its two-dimensional
span is the identifiable object. Increasing `ndim` to two would not be
an equivalent way to obtain this picture: it would estimate two local
tangent directions and construct a different operator.

## COIL-20: the same issue in real images

[COIL-20](https://cave.cs.columbia.edu/repository/COIL-20) contains 72
views of each of 20 objects as it rotates. Each object was fitted
separately, avoiding the disconnected graph obtained by mixing objects
at a small neighborhood. Before any embedding was inspected, all objects
were screened at `k = 7` using only their pixel-neighbor graphs. The
first two objects under a fixed rule were objects 4 and 1: both direct
and effective graphs had to be connected, most neighbors had to be
locally adjacent in pose, each pose needed neighbors on both sides, and
the graph needed wraparound edges.

The fits used ordinary LTSA, `ndim = 1`, exact self-first neighbors,
four retained modes, and `eig_k = 8`. No supervised rotation was
applied. Pose is used only for color and for choosing eight evenly
spaced explanatory images.

![COIL-20 object 4. Raw modes 1 and 2 form a pose-ordered oval; the
first spectral pair is much more clearly separated from mode 3 than
either member is from the
other.](figures/ltsa-spectral-block-coil-object-4.png)

COIL-20 object 4. Raw modes 1 and 2 form a pose-ordered oval; the first
spectral pair is much more clearly separated from mode 3 than either
member is from the other.

![COIL-20 object 1. The path is more distorted, but the independently
selected object still gives a coherent periodic progression rather than
a hand-picked geometric
circle.](figures/ltsa-spectral-block-coil-object-1.png)

COIL-20 object 1. The path is more distorted, but the independently
selected object still gives a coherent periodic progression rather than
a hand-picked geometric circle.

The scaled mode 1/2 gaps were `0.000108` and `0.000181`, while the pair
2/3 gaps were `0.00250` and `0.00215`. In both examples, the second
coordinate turns a one-dimensional slice into an intelligible rotational
orbit already present in the fixed operator.

The essential reproduction is short; the graph-only screening and
thumbnail layout are article-specific presentation code:

``` r

coil <- snedata::download_coil20(as = "list")
rows <- coil$meta$object == 4
object_4 <- coil$data[rows, , drop = FALSE]

coil_fit <- ltsa(
  object_4,
  ndim = 1,
  n_neighbors = 7,
  nn_method = "exact",
  eig_method = "eig",
  eig_k = 8,
  output = "result",
  spectral_dim = 4
)

plot(coil_fit$spectral_embedding[, 1:2])
```

This is evidence for inspecting a retained block, not evidence that
COIL-20 as a whole has a useful one-shot LTSA embedding.

## MNIST: useful variation and a warning

The complete 7,877-observation digit-1 subset of
[MNIST](https://github.com/fgnt/mnist) gives a less tidy, more realistic
example. The subset and the use of normalized LTSA were already
suggested by earlier exploratory work, so this is disclosed pedagogical
evidence rather than an independent benchmark.

Ordinary and normalized LTSA were fit to the same materialized
approximate neighbor graph with `ndim = 2`, `k = 15`, four retained
modes, and `eig_k = 8`. Normalized LTSA solves a different generalized
eigenproblem; it is not a numerical repair for ordinary LTSA.

The neighbor graph was generated once with NND, seed `20260824`, and one
neighbor thread. Both fits reused that self-first graph, used one
assembly thread and RSpectra, and reset the solver seed to `20260824`
immediately before each fit.

The ordinary first three modes were almost point-localized despite a
connected effective graph, full local ranks, and acceptable solver
residuals. Normalized modes 1 and 2 were broadly supported and gave a
coherent view of slant, curvature, stroke width, hooks, and unusual
serifed or branched forms.

![Normalized LTSA modes 1 and 2 for MNIST digit 1. The numbered images
were selected to cover the central plotted view; they are explanatory
representatives, not a quality
sample.](figures/ltsa-spectral-block-mnist-normalized-display.png)

Normalized LTSA modes 1 and 2 for MNIST digit 1. The numbered images
were selected to cover the central plotted view; they are explanatory
representatives, not a quality sample.

Modes 2 and 3 give a complementary view of that morphology, but mode 3
is already noticeably localized. The picture is interesting; the
diagnostic is the reason not to promote it uncritically as a better
embedding.

![Normalized LTSA modes 2 and 3. The additional mode exposes a different
sheet-and-branch view, but its limited support makes this a qualified
diagnostic view rather than an automatic
improvement.](figures/ltsa-spectral-block-mnist-normalized-extra.png)

Normalized LTSA modes 2 and 3. The additional mode exposes a different
sheet-and-branch view, but its limited support makes this a qualified
diagnostic view rather than an automatic improvement.

![The MNIST spectra and effective support fractions. Ordinary modes 1–3
are near the single-observation support limit; normalized modes become
progressively more localized after the first
pair.](figures/ltsa-spectral-block-mnist-diagnostics.png)

The MNIST spectra and effective support fractions. Ordinary modes 1–3
are near the single-observation support limit; normalized modes become
progressively more localized after the first pair.

A compact reproduction of the normalized retained block is:

``` r

mnist <- snedata::download_mnist(as = "list")
labels <- as.integer(as.character(mnist$meta$label))
digit_1 <- mnist$data[labels == 1, , drop = FALSE]
nonconstant <- vapply(
  seq_len(ncol(digit_1)),
  function(index) {
    pixel <- digit_1[, index]
    min(pixel) != max(pixel)
  },
  logical(1)
)
digit_1 <- digit_1[, nonconstant, drop = FALSE]

set.seed(20260824)
mnist_fit <- ltsa(
  digit_1,
  ndim = 2,
  n_neighbors = 15,
  nn_method = "nnd",
  n_threads = 1,
  eig_method = "rspectra",
  eig_k = 8,
  normalize = TRUE,
  output = "result",
  spectral_dim = 4
)
```

For a direct comparison of ordinary and normalized LTSA, materialize one
neighbor graph and pass that same self-first index matrix as `nn_method`
to both calls. At minimum, reset the seed immediately before each
sampled call and keep the thread settings identical.

The representative rule standardizes the plotted pair, keeps the central
98% per axis, starts nearest the component-wise median, and greedily
adds the point farthest from those already selected. It covers the view
deliberately. It does not measure embedding quality.

The mode-support diagnostic used here is also deliberately local to the
article. For a Euclidean-normalized coordinate vector `v`, its effective
support fraction is

``` math
\frac{1}{n\sum_i v_i^4}.
```

It equals one for a perfectly uniform-magnitude vector and approaches
`1 / n` when one observation dominates:

``` r

effective_support_fraction <- function(v) {
  v <- v / sqrt(sum(v^2))
  1 / (length(v) * sum(v^4))
}

apply(mnist_fit$spectral_embedding, 2, effective_support_fraction)
```

This is not a generic embedding-quality score. It answers the narrower
question raised by these plots: is a coordinate broadly supported, or is
an attractive central cloud produced after a few observations absorb
almost all its energy?

## Fashion-MNIST: extra modes need not help

Three fixed, previously inspected classes from
[Fashion-MNIST](https://github.com/zalandoresearch/fashion-mnist)—trouser,
dress, and bag—supply the negative example. Every class was fitted with
ordinary LTSA using the same declared `ndim = 2`, `k = 15`, and four
retained modes. No class, neighborhood, normalization setting, or mode
pair was substituted after seeing the results.

Each class used its own self-first NND graph, with seed `20260825` reset
before graph construction and one neighbor thread. Every RSpectra fit
used `eig_k = 8`, one assembly thread, and the same seed reset
immediately before eigenanalysis.

![Displayed pairs for the three fixed Fashion-MNIST classes. Central
percentile framing makes the majority cloud visible, but the
representatives do not establish a semantic
ordering.](figures/ltsa-spectral-block-fashion-comparison.png)

Displayed pairs for the three fixed Fashion-MNIST classes. Central
percentile framing makes the majority cloud visible, but the
representatives do not establish a semantic ordering.

All twelve retained coordinates were essentially at the
single-observation support limit, and more than 99.96% of each
coordinate’s energy lay in the largest 1% of observations. The graphs
were connected, local patches had full rank, and the eigensolves
converged. The low modes themselves were localized.

![Fashion-MNIST low spectra and effective support. Weak boundaries
motivate caution, but retaining extra modes does not rescue globally
unsupported
coordinates.](figures/ltsa-spectral-block-fashion-diagnostics.png)

Fashion-MNIST low spectra and effective support. Weak boundaries
motivate caution, but retaining extra modes does not rescue globally
unsupported coordinates.

This is the crucial counterexample: a weak boundary is a reason not to
treat the first displayed plane as uniquely privileged. It is not a
promise that a better raw pair exists in the retained block.

## A practical inspection sequence

When an LTSA result is scientifically or visually important:

1.  Inspect `eigen$status`, `eigen$messages`, and the displayed boundary
    gap.
2.  Check `assembly$component_count` and the local-rank diagnostics. A
    spectral block cannot repair a disconnected effective graph or
    deficient patches.
3.  If the displayed boundary is weak or the data have periodic
    structure, retain a small same-operator block with `spectral_dim`
    rather than increasing `ndim`.
4.  Inspect raw pair views without assuming that one is optimal. Treat
    signs and rotations inside a nearly repeated eigenspace as
    arbitrary.
5.  Check whether each plotted mode is broadly supported. An accurately
    solved, low-energy vector may still describe only a few
    observations.
6.  If comparing ordinary and normalized LTSA, name them as different
    estimators and keep the graph and stochastic settings fixed.

The result may be a useful complementary view, as for the circle,
COIL-20, and the qualified MNIST example. It may instead show that the
low-frequency block does not contain a defensible global chart, as in
the Fashion-MNIST subsets. Both outcomes are informative. `spectral_dim`
exposes the evidence; it does not choose the interpretation.

For details of the returned diagnostics and eigensolver controls, see
[Numerical diagnostics and solver
notes](https://jlmelville.github.io/flotsam/articles/numerical-diagnostics.md).

## Data sources and attribution

The image examples come from
[COIL-20](https://cave.cs.columbia.edu/repository/COIL-20) (Sameer A.
Nene, Shree K. Nayar, and Hiroshi Murase),
[MNIST](https://yann.lecun.com/exdb/mnist/) (Yann LeCun, Corinna Cortes,
and Christopher J. C. Burges), and
[Fashion-MNIST](https://arxiv.org/abs/1708.07747) (Han Xiao, Kashif
Rasul, and Roland Vollgraf). The package contains only the small,
low-resolution samples inside these static explanatory figures;
article-development scripts download the datasets separately.
