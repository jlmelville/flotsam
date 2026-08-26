# Local Tangent Space Alignment

LTSA turns overlapping local tangent coordinates into a global
embedding. This article explains that path and, in particular, which
controls change the neighborhood graph, the alignment operator, the
candidate eigensolve, or only the returned view. Use
[`?ltsa`](https://jlmelville.github.io/flotsam/reference/ltsa.md) for
exhaustive defaults and validation rules.

## From neighborhoods to an embedding

Let `X` be an $`N \times p`$ data matrix: $`N`$ observations in $`p`$
input features. Let $`k`$ be the neighborhood size and let $`d`$ be the
requested local tangent dimension, supplied as `ndim`.

1.  Find a $`k`$-observation neighborhood for each row of `X`.

2.  Center each $`k \times p`$ neighborhood matrix and compute its
    singular value decomposition.

3.  Retain the first $`d`$**left** singular vectors in $`U_i`$. These
    give local tangent coordinates for the observations in neighborhood
    $`i`$.

4.  Form the local residual projector

    ``` math
    W_i = I_k - \mathbf{1}\mathbf{1}^{\mathsf T}/k - U_iU_i^{\mathsf T}.
    ```

5.  Add the entries of $`W_i`$ into the corresponding rows and columns
    of an initially zero $`N \times N`$ matrix $`B`$, then repeat for
    every neighborhood.

6.  Remove the known constant null direction and use the lowest
    nonconstant directions of $`B`$ for the global coordinates.

The assembled $`B`$ is a symmetric positive-semidefinite alignment
operator. Its constant vector is a known null direction, but its
nullspace need not be one-dimensional. Disconnected effective
neighborhoods or other degeneracies can therefore make “the smallest
eigenvector” non-unique.

This construction explains why `ndim` is more than a plotting choice. It
sets the local tangent rank used inside every $`W_i`$, so changing it
rebuilds $`B`$. It also sets the number of coordinates in the ordinary
returned `embedding`.

## Follow the controls through the pipeline

``` text
data + neighborhood controls + ndim
                  |
                  v
         alignment matrix B
                  |
          normalize, if requested
                  |
                  v
      eig_method + eig_k candidates
                  |
                  v
      spectral_dim retained modes
                  |
                  v
       first ndim modes displayed
```

| Control | Stage it changes | Consequence |
|----|----|----|
| `n_neighbors`, `nn_method`, `include_self` | Neighborhoods and assembly | Change which local projectors contribute to $`B`$ |
| `ndim` | Local tangent bases, $`B`$, and displayed prefix | Changes the estimator as well as the embedding dimension |
| `normalize` | Eigenproblem built from $`B`$ | Selects ordinary or generalized LTSA |
| `eig_method` and backend controls | Candidate computation | Change how the fixed eigenproblem is solved |
| `eig_k` | Candidate span | Gives Rayleigh–Ritz more directions without changing $`B`$ |
| `spectral_dim` | Retained fixed-operator block | Keeps extra modes without changing $`B`$ or `ndim` |
| `output`, `include_B` | Returned evidence | Select an embedding, detailed result, operator, or included raw $`B`$ |

The most consequential distinction is among `ndim`, `spectral_dim`, and
`eig_k`:

- Increase `ndim` only when the local manifold model itself should have
  a higher tangent dimension. This rebuilds the alignment operator.
- Increase `spectral_dim` with `output = "result"` when you want to
  inspect more modes from the existing operator. The ordinary
  `embedding` remains `ndim`-dimensional; the larger block is returned
  as `spectral_embedding`.
- Increase `eig_k` when the numerical candidate span may be too narrow.
  It does not retain more output by itself.

For a manual expanded request, `eig_k` must leave room for the known
constant direction and a mode beyond the retained boundary, so it must
be at least `spectral_dim + 2`. The default chooses a suitable candidate
width; see
[`?ltsa`](https://jlmelville.github.io/flotsam/reference/ltsa.md) for
its exact formula and all dimensional limits.

## Choose neighborhoods before solver settings

`nn_method = "nnd"` uses approximate nearest-neighbor descent and is the
large-data default. `nn_method = "exact"` exhaustively compares
observations and is useful for smaller deterministic cases. A
precomputed 1-based neighbor index matrix, or an object with an `idx`
matrix, lets several fits reuse one fixed graph. The reference page
documents the required self-first shapes.

`n_neighbors` is the effective neighborhood size used in assembly and
must be at least `ndim + 2`: after removing the constant and local
tangent directions, each local projector needs a remaining residual
direction. Too-small or poorly overlapping neighborhoods can produce
disconnected effective-neighborhood components or weakly separated
global directions. Larger neighborhoods change the local approximation
rather than merely making the same calculation more accurate.

With `include_self = FALSE`, a precomputed self-first matrix therefore
has `n_neighbors + 1` columns: assembly removes the self column and uses
the remaining `n_neighbors` observations. This changes the effective
neighborhood definition; [work by Zhang and
co-workers](https://doi.org/10.1109/JSTARS.2017.2682189) relates that
choice to Hessian Locally Linear Embedding. Treat it as an estimator
choice, not a reproducibility switch.

## Choose the eigenproblem, then its computation

Ordinary LTSA diagonalizes $`B`$. With `normalize = TRUE`, `flotsam`
instead solves

``` math
Bv = \lambda Dv, \qquad D = \operatorname{diag}(B),
```

through $`D^{-1/2}BD^{-1/2}`$ and maps the selected coordinates back.
The normalized formulation uses `D`-weighted orthogonality and
centering. It is a different generalized estimator, not a
conditioning-only route to the same coordinates, and `D` is not ordinary
graph degree.

After that choice, `eig_method = "auto"` uses dense
[`base::eigen()`](https://rdrr.io/r/base/eigen.html) for small or
proportionally wide requests and otherwise uses RSpectra. Explicit
`"rspectra"`, `"irlba"`, `"svdr"`, and dense `"eig"` requests always use
the named backend. The detailed result keeps the requested policy in
`eigen$method` and the backend that actually ran in
`eigen$backend$name`.

Backend convergence controls, dense-route thresholds, and diagnostic
tolerances belong on the troubleshooting path rather than in this
conceptual overview. See [Diagnosing suspicious LTSA
results](https://jlmelville.github.io/flotsam/articles/numerical-diagnostics.md)
for the action associated with each status, message, gap, component, and
rank field.

## Decide what evidence to return

The default `output = "embedding"` returns the $`N \times`$`ndim`
coordinate matrix. Use `output = "result"` when you need the evidence
required to assess that map:

``` r

fit <- ltsa(X, output = "result", spectral_dim = 4)

fit$embedding
fit$spectral_embedding
fit$eigen$status
fit$eigen$messages
fit$assembly$component_count
fit$assembly$rank_deficient_count
```

When `spectral_dim > ndim`, `eigen` describes the displayed prefix and
`eigen$spectral` describes the retained block. The displayed embedding
is exactly the first `ndim` columns of `spectral_embedding`.

Use `output = "B"` to stop after the operator is ready. It returns raw
$`B`$ for ordinary LTSA and the normalized operator used by
eigenanalysis when `normalize = TRUE`. Use `include_B = TRUE` with
`output = "result"` when a detailed result should also carry the
assembled, unnormalized $`B`$.

Thread counts change resource use rather than the intended mathematical
stage, with one important caveat: multithreaded approximate neighbor
search need not reproduce an identical graph. `n_threads` controls
computed-neighbor search; `n_assembly_threads` controls construction of
$`B`$. Serial settings are the reproducible default. The [threading
guidance](https://jlmelville.github.io/flotsam/articles/numerical-diagnostics.html#threading-and-cross-implementation-checks)
also covers threaded-BLAS oversubscription.

For a figure-led example of why local dimension and useful global
coordinate count can differ, continue to [Exploring LTSA spectral
blocks](https://jlmelville.github.io/flotsam/articles/spectral-blocks.md).
Its circle example has a one-dimensional local tangent but needs a
two-mode span to show the closed path.
