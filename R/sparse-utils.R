# get indices of the diagonal from the sparse matrix m
diag_spm <- function(m) {
  is <- m@i
  ps <- m@p
  n <- nrow(m)
  sapply(
    1:n,
    function(i, is, ps) {
      begin <- ps[i] + 1
      end <- ps[i + 1]
      begin + which(is[begin:end] == i - 1)
    },
    is,
    ps
  ) -
    1
}

lsym_norm <- function(M, D) {
  M@x <- spm_times_scalar(M@p, M@x, D)
  D * M
}

# vI - m
shift_lap <- function(m, v = 2.0) {
  x <- m@x

  x <- -x
  ds <- diag_spm(m)
  x[ds] <- x[ds] + v

  m@x <- x
  m
}

norm_and_shift_L <- function(L) {
  nres <- norm_lsym_L(L)
  list(Lshift = shift_lap(nres$Lsym, 2.0), Dinvs = nres$Dinvs)
}

norm_lsym_L <- function(L) {
  Dinvs <- sqrt(1 / diag(L))
  list(
    Lsym = lsym_norm(L, Dinvs),
    Dinvs = Dinvs,
    nullvec = 1 / Dinvs
  )
}
