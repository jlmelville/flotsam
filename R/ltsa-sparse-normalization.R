ltsa_normalize_sparse_operator <- function(L) {
  diagonal <- diag(L)
  if (any(!is.finite(diagonal) | diagonal <= 0)) {
    stop(
      "Cannot normalize the LTSA matrix because its diagonal contains ",
      "non-positive or non-finite entries. Increase n_neighbors or set ",
      "normalize = FALSE.",
      call. = FALSE
    )
  }

  Dinvs <- sqrt(1 / diagonal)
  L_scaled <- L
  L_scaled@x <- spm_times_scalar(L_scaled@p, L_scaled@x, Dinvs)
  list(
    Lsym = Dinvs * L_scaled,
    mass = diagonal,
    Dinvs = Dinvs,
    nullvec = 1 / Dinvs
  )
}

ltsa_normalized_details <- function(
  B,
  mass,
  symmetric_embedding,
  embedding,
  values,
  lambda_max,
  reverse_occurrence,
  component_membership,
  component_sizes
) {
  symmetric_embedding <- as.matrix(symmetric_embedding)
  embedding <- as.matrix(embedding)

  # For u = sqrt(D) v and Lsym = D^(-1/2) B D^(-1/2),
  #
  #   B v - lambda D v = sqrt(D) (Lsym u - lambda u).
  #
  # The generalized absolute residual is ||B v - lambda D v||_2. Its
  # reported scaled form divides by
  #
  #   sqrt(max(D)) * max(lambda_max, 1),
  #
  # the operator norm of the map-back residual scaling times the existing
  # symmetric-problem residual scale. The existing symmetric residuals remain
  # in eigen$residuals.
  BV <- as.matrix(B %*% embedding)
  generalized_residual <- BV - mass * sweep(embedding, 2L, values, "*")
  generalized_absolute_residuals <- sqrt(
    colSums(generalized_residual * generalized_residual)
  )
  generalized_residual_scale <- sqrt(max(mass)) *
    ltsa_residual_scale(lambda_max)

  weighted_gram_error <- crossprod(embedding, mass * embedding) -
    diag(ncol(embedding))
  weighted_centering_error <- crossprod(embedding, mass)
  map_back_difference <- symmetric_embedding - sqrt(mass) * embedding
  # Each scalar algebra error below is the maximum absolute entry of the
  # displayed matrix/vector identity.

  # reverse_occurrence[j] counts the local projector diagonal terms that
  # contribute to mass[j]. Each contribution is its residual leverage
  # W_i[ell, ell], so the count describes participation but is not D itself.
  reverse_mass_correlation <- ltsa_defined_correlation(
    reverse_occurrence,
    mass
  )

  list(
    mass = mass,
    mass_summary = ltsa_mass_summary(mass),
    symmetric_embedding = symmetric_embedding,
    generalized_absolute_residuals = generalized_absolute_residuals,
    generalized_residual_scale = generalized_residual_scale,
    generalized_residuals = generalized_absolute_residuals /
      generalized_residual_scale,
    weighted_orthogonality_error = max(abs(weighted_gram_error)),
    weighted_centering_error = max(abs(weighted_centering_error)),
    map_back_error = max(abs(map_back_difference)),
    reverse_occurrence = list(
      counts = as.integer(reverse_occurrence),
      quantiles = ltsa_bounded_quantiles(reverse_occurrence),
      correlation_with_mass = reverse_mass_correlation
    ),
    # In symmetric coordinates this compares u with sqrt(D)-weighted
    # component indicators after removing sqrt(D) 1. It is equivalent to
    # comparing v with component indicators in the D inner product.
    component_embedding_overlap = ltsa_component_embedding_overlap(
      symmetric_embedding,
      component_membership,
      component_sizes,
      point_weights = mass
    )
  )
}

ltsa_bounded_quantiles <- function(x) {
  stats::setNames(
    as.numeric(stats::quantile(
      x,
      probs = c(0, 0.25, 0.5, 0.75, 1),
      names = FALSE
    )),
    c("minimum", "lower_quartile", "median", "upper_quartile", "maximum")
  )
}

ltsa_mass_summary <- function(mass, index_limit = 8L) {
  quantiles <- ltsa_bounded_quantiles(mass)
  smallest_order <- order(mass, seq_along(mass), method = "radix")
  largest_order <- order(-mass, seq_along(mass), method = "radix")
  displayed <- min(length(mass), index_limit)

  list(
    quantiles = quantiles,
    # Ratios are returned directly without clipping. The log10 range remains
    # finite when max / min overflows, so callers can see the true scale span.
    max_to_min_ratio = max(mass) / min(mass),
    min_to_median_ratio = min(mass) / quantiles[["median"]],
    log10_max_to_min_ratio = log10(max(mass)) - log10(min(mass)),
    log10_min_to_median_ratio = log10(min(mass)) -
      log10(quantiles[["median"]]),
    index_limit = as.integer(index_limit),
    smallest_mass_indices = head(smallest_order, displayed),
    largest_mass_indices = head(largest_order, displayed)
  )
}

ltsa_defined_correlation <- function(left, right) {
  if (
    length(left) != length(right) ||
      length(left) < 2L ||
      any(!is.finite(left)) ||
      any(!is.finite(right)) ||
      all(left == left[[1L]]) ||
      all(right == right[[1L]])
  ) {
    return(NA_real_)
  }

  correlation <- suppressWarnings(stats::cor(left, right))
  if (length(correlation) != 1L || !is.finite(correlation)) {
    return(NA_real_)
  }
  as.numeric(correlation)
}
