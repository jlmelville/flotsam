# Numerical diagnostics and solver notes

This note collects the LTSA implementation details that are useful when
a result looks numerically suspicious, but too long for the README or
[`?ltsa`](https://jlmelville.github.io/flotsam/reference/ltsa.md).

## What `ltsa()` is solving

LTSA builds a sparse alignment matrix `B` from local neighborhoods. For
each observation, the package finds a neighborhood, computes a local
tangent-space basis, forms a small local weight matrix, and scatters
those weights into the global sparse matrix. The embedding comes from
the lowest nonconstant eigenvectors of that matrix.

This is deliberately close to other spectral manifold methods: if the
neighborhood graph is disconnected, weakly connected, or has a
low-energy cluster that cuts across the requested `ndim`, the final
eigenspace may be ambiguous.

## Fixed-width eigenanalysis

The iterative backends do not ask directly for only `ndim` output
vectors. They request a wider candidate block controlled by `eig_k`,
remove the known constant null direction, and run a small Rayleigh-Ritz
extraction inside that candidate span. The final embedding uses the
first `ndim` nonconstant vectors from that projected solve.

The default `eig_k` is intentionally wider than `ndim + 1`. Increasing
`eig_k` can help when diagnostics suggest that the candidate span did
not cover enough of a clustered low-energy eigenspace, but it also costs
more work.

## Reading result diagnostics

Use `output = "result"` when you need to inspect a solve:

``` r

res <- ltsa(X, output = "result")
res$eigen[c("status", "eig_k", "rank")]
res$eigen$messages
```

The most useful public fields are:

- `status`: solver-neutral classification of the requested solve.
- `messages`: short explanations for warning or invalid classifications.
- `eig_k`: number of backend candidate vectors requested.
- `values`: selected embedding eigenvalues.
- `ritz_values`: values available after null-vector projection inside
  the candidate span.
- `residuals`: scaled residuals for the selected vectors.
- `rank`: post-null rank available in the candidate span.
- `backend`: backend metadata such as convergence counts where
  available.

These diagnostics classify the requested solve; they are not
mathematical certificates for the entire low-energy eigenspace.

## Backend notes

`eig_method = "rspectra"` is the default. It reports RSpectra
convergence metadata, so hard backend convergence failures can be
distinguished from post-hoc residual or boundary diagnostics.

`eig_method = "irlba"` and `eig_method = "svdr"` use irlba routines to
produce candidate vectors. They do not expose the same convergence
metadata as RSpectra, so `flotsam` relies on residual and eigenspace
diagnostics after the candidate block is returned.

`eig_method = "eig"` and `"eigen"` use dense
[`base::eigen()`](https://rdrr.io/r/base/eigen.html) and are intended as
small diagnostic references.

The RSpectra path uses shifted largest-algebraic solves rather than
asking ARPACK/Spectra to work directly at zero. That approach was
motivated by practical issues around small clustered eigenvalues and by
discussion between Aaron Lun, Kevin Doherty, and Yixuan Qiu in
<https://github.com/yixuan/spectra/issues/126>.

## Normalization

`normalize = TRUE` solves a normalized LTSA formulation. It can be
useful for comparison, but it is a different spectral objective from the
default unnormalized alignment matrix. If a downstream task depends on
geometry or clustering in the embedding, compare diagnostics and
behavior before relying on the normalized result.

## Performance notes

Nearest-neighbor search is controlled by `n_threads`. Sparse
alignment-matrix assembly is controlled separately by
`n_assembly_threads`.

When `n_assembly_threads > 1`, avoid oversubscribing the machine with an
already-multithreaded BLAS/LAPACK setup. In practice that often means
using outer LTSA assembly threads while keeping BLAS threads low, or
keeping LTSA assembly serial while BLAS is free to use threads.

`copy_max_mib` caps an optional row-major copy of `X` used by the
high-dimensional local Gram route. Set it to `0` to disable that copy
when peak memory matters more than local-weight speed.

## Why spectral methods can still be finicky

Some spectral manifold-learning failures are not implementation bugs.
Very long, thin, uneven, disconnected, or weakly connected manifolds can
make the normalization and eigenspace selection problem ill-conditioned.
This is one reason `flotsam` keeps compact diagnostics visible instead
of returning only an embedding and hoping for the best.

For more background on spectral caveats, see Goldberg, Zakai, Kushnir,
and Ritov (2008), *Manifold learning: The price of normalization*, and
the other references listed in the README.
