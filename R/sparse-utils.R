lsym_norm <- function(M, D) {
  M@x <- spm_times_scalar(M@p, M@x, D)
  D * M
}

norm_lsym_L <- function(L) {
  Dinvs <- sqrt(1 / diag(L))
  list(
    Lsym = lsym_norm(L, Dinvs),
    Dinvs = Dinvs,
    nullvec = 1 / Dinvs
  )
}
