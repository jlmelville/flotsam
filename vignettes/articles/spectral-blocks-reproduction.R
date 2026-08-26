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
