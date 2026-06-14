# Legacy sparse-slot-search assembly retained for internal regression tests.
assemble_ltsa_B_slot_search_reference <- function(X,
                                                  nn_idx,
                                                  ndim,
                                                  include_self,
                                                  verbose = FALSE) {
  n <- nrow(X)
  B <- create_sparse(nn_idx, verbose = verbose)
  Bx <- B@x
  weight_idx <- ltsa_weight_neighborhoods(nn_idx, include_self)

  prev_time <- Sys.time()
  rank_deficient_count <- 0L
  min_local_rank <- ndim
  for (i in seq_len(n)) {
    if (verbose) {
      curr_time <- Sys.time()
      if (difftime(curr_time, prev_time, units = "secs") > 60) {
        tsmessage("processed ", i, " / ", n)
        prev_time <- curr_time
      }
    }

    nni <- weight_idx[i, ]
    local <- ltsa_local_weights(X, nni, ndim)
    if (local$rank < ndim) {
      rank_deficient_count <- rank_deficient_count + 1L
      min_local_rank <- min(min_local_rank, local$rank)
    }

    spi <- find_in_spsq(B, nni)
    Bx[spi] <- Bx[spi] + local$Wi
  }
  B@x <- Bx

  list(
    B = B,
    rank_deficient_count = rank_deficient_count,
    min_local_rank = min_local_rank
  )
}

create_sparse <- function(nn_idx, verbose = FALSE) {
  n <- nrow(nn_idx)
  ij <- nbrhood_triplets(t(nn_idx), ncol(nn_idx))
  m <- methods::new(
    "dgTMatrix",
    i = ij$i,
    j = ij$j,
    x = rep(0, length(ij$i)),
    Dim = as.integer(c(n, n))
  )
  m <- m + Matrix::t(m)
  methods::as(m, "CsparseMatrix")
}

find_in_spsq <- function(m, is) {
  sparse_idxs(m@i, m@p, is)
}
