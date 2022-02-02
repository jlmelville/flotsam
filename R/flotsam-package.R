## usethis namespace: start
#' @import Matrix
#' @useDynLib flotsam, .registration = TRUE
## usethis namespace: end
.onUnload <- function(libpath) {
  library.dynam.unload("flotsam", libpath)
}
