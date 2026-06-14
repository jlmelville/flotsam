#' flotsam: Local Tangent Space Alignment
#'
#' `flotsam` implements Local Tangent Space Alignment (LTSA) for nonlinear
#' dimensionality reduction.
#'
#' The main entry point is [ltsa()]. It accepts numeric matrices and data frames
#' with observations in rows.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import Matrix
#' @useDynLib flotsam, .registration = TRUE
## usethis namespace: end
.onUnload <- function(libpath) {
  library.dynam.unload("flotsam", libpath)
}
