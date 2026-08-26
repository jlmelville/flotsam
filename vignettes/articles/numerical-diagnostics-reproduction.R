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
