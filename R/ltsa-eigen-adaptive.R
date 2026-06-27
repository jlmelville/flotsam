# Shared adaptive Ritz driver and public-facing iterative wrappers for LTSA
# eigenanalysis.

# Adaptive driver shared by rspectra, irlba, and svdr. It over-requests
# candidates, re-extracts LTSA Ritz vectors, expands only while diagnostics
# justify it, and returns the first residual-good/rank-good weak-gap result
# after the configured diagnostic expansion budget is exhausted.
ltsa_adaptive_ritz_eig <- function(
  B,
  ndim,
  provider,
  provider_args = list(),
  nullvec = ltsa_default_null_vector(nrow(B)),
  initial_extra = 4L,
  max_extra = 40L,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  gap_expansion_steps = 1L,
  strict_rescue = TRUE,
  strict_rescue_arg_mapper = function(provider_args, ...) provider_args,
  strict_rescue_controls = list(),
  strict_rescue_extra = 5L,
  width_first_rescue = FALSE,
  width_first_rescue_max_expansions = 2L,
  attempt_reference_vectors = NULL,
  retain_attempt_candidate_spaces = FALSE,
  bank_backend = NULL,
  verbose = FALSE
) {
  n <- nrow(B)
  if (!is.function(provider)) {
    stop("LTSA candidate provider must be a function", call. = FALSE)
  }
  if (is.null(provider_args)) {
    provider_args <- list()
  }
  if (!is.list(provider_args)) {
    stop("LTSA candidate provider arguments must be a list", call. = FALSE)
  }
  if (!is.function(strict_rescue_arg_mapper)) {
    stop("LTSA strict-rescue argument mapper must be a function", call. = FALSE)
  }
  if (is.null(strict_rescue_controls)) {
    strict_rescue_controls <- list()
  }
  if (!is.list(strict_rescue_controls)) {
    stop("LTSA strict-rescue controls must be a list", call. = FALSE)
  }
  eig_k <- ltsa_initial_ritz_candidate_k(ndim, n, initial_extra)
  max_k <- ltsa_max_ritz_candidate_k(ndim, n, max_extra)
  gap_expansion_steps <- as.integer(gap_expansion_steps)
  if (
    length(gap_expansion_steps) != 1L ||
      is.na(gap_expansion_steps) ||
      gap_expansion_steps < 0L
  ) {
    stop("gap_expansion_steps must be a non-negative integer", call. = FALSE)
  }
  strict_rescue <- as.logical(strict_rescue)
  if (length(strict_rescue) != 1L || is.na(strict_rescue)) {
    stop("strict_rescue must be TRUE or FALSE", call. = FALSE)
  }
  strict_rescue_extra <- as.integer(strict_rescue_extra)
  if (
    length(strict_rescue_extra) != 1L ||
      is.na(strict_rescue_extra) ||
      strict_rescue_extra < 1L
  ) {
    stop("strict_rescue_extra must be a positive integer", call. = FALSE)
  }
  width_first_rescue <- as.logical(width_first_rescue)
  if (length(width_first_rescue) != 1L || is.na(width_first_rescue)) {
    stop("width_first_rescue must be TRUE or FALSE", call. = FALSE)
  }
  width_first_rescue_max_expansions <-
    as.integer(width_first_rescue_max_expansions)
  if (
    length(width_first_rescue_max_expansions) != 1L ||
      is.na(width_first_rescue_max_expansions) ||
      width_first_rescue_max_expansions < 0L
  ) {
    stop(
      "width_first_rescue_max_expansions must be a non-negative integer",
      call. = FALSE
    )
  }
  width_first_rescue <- isTRUE(width_first_rescue) && isTRUE(strict_rescue)
  lambda_max <- NULL
  attempts <- list()
  best <- NULL
  weak_gap_candidate <- NULL
  first_weak_gap_eig_k <- NA_integer_
  last_error <- NULL
  gap_expansions <- 0L
  width_rescue_trigger_eig_k <- NA_integer_
  width_rescue_incumbent <- NULL

  finish_width_first_rescue <- function(selected, unresolved = FALSE) {
    selected$acceptance$width_first_rescue_trigger_eig_k <-
      width_rescue_trigger_eig_k
    selected$acceptance$width_first_rescue_ordinary_expansions <-
      ltsa_width_rescue_extra_count(
        attempts = selected$attempts,
        trigger_eig_k = width_rescue_trigger_eig_k
      )
    ordinary_eig_ks <- ltsa_ordinary_attempt_eig_ks(selected$attempts)
    if (length(ordinary_eig_ks) > 0L) {
      selected$acceptance$diagnostic_final_eig_k <- max(ordinary_eig_ks)
    }
    if (isTRUE(unresolved)) {
      selected$acceptance$return_reason <-
        "width_first_rescue_unresolved_partial"
      selected <- ltsa_maybe_warn_width_first_rescue_unresolved(
        selected,
        ndim = ndim
      )
    } else {
      selected <- ltsa_maybe_warn_partial_near_zero_block(
        selected,
        ndim = ndim,
        strict_rescue = strict_rescue
      )
    }
    selected <- ltsa_maybe_message_width_first_rescue_decision(
      selected,
      ndim = ndim,
      verbose = verbose
    )
    selected <- ltsa_maybe_warn_spectral_ambiguity(selected, ndim)
    ltsa_finalize_attempt_observability(
      selected,
      reference_vectors = attempt_reference_vectors,
      nullvec = nullvec,
      retain_candidate_spaces = retain_attempt_candidate_spaces
    )
  }

  finish_width_rescue_if_exhausted <- function() {
    if (is.null(width_rescue_incumbent)) {
      return(NULL)
    }

    extra_count <- ltsa_width_rescue_extra_count(
      attempts = attempts,
      trigger_eig_k = width_rescue_trigger_eig_k
    )
    if (length(attempts) > 0L) {
      updated_attempts <- attempts
      last_attempt <- length(updated_attempts)
      updated_attempts[[last_attempt]]$width_first_rescue_expansion <-
        extra_count
      attempts <<- updated_attempts
    }
    if (
      extra_count < width_first_rescue_max_expansions &&
        eig_k < max_k
    ) {
      return(NULL)
    }

    selected <- width_rescue_incumbent
    selected$attempts <- attempts
    selected <- ltsa_maybe_staged_strict_rescue(
      B = B,
      selected = selected,
      ndim = ndim,
      nullvec = nullvec,
      lambda_max = selected$lambda_max %||% lambda_max,
      provider = provider,
      provider_args = provider_args,
      strict_rescue_arg_mapper = strict_rescue_arg_mapper,
      strict_rescue_controls = strict_rescue_controls,
      resid_tol = resid_tol,
      gap_tol = gap_tol,
      bank_backend = bank_backend,
      verbose = verbose
    )
    unresolved <- ltsa_partial_near_zero_block(selected, ndim) &&
      !isTRUE(selected$acceptance$width_first_rescue_used)
    finish_width_first_rescue(selected, unresolved = unresolved)
  }

  repeat {
    candidate_elapsed <- NA_real_
    candidate_timing <- system.time({
      eig_res <- tryCatch(
        ltsa_call_candidate_provider(
          provider = provider,
          B = B,
          eig_k = eig_k,
          provider_args = provider_args,
          lambda_max = lambda_max,
          verbose = verbose
        ),
        error = function(e) e
      )
    })
    candidate_elapsed <- unname(candidate_timing[["elapsed"]])
    if (inherits(eig_res, "error")) {
      last_error <- conditionMessage(eig_res)
      attempts[[length(attempts) + 1L]] <- ltsa_attempt_error_summary(
        eig_k = eig_k,
        error = last_error,
        candidate_elapsed = candidate_elapsed
      )
      if (width_first_rescue && !is.null(width_rescue_incumbent)) {
        rescued <- finish_width_rescue_if_exhausted()
        if (!is.null(rescued)) {
          return(rescued)
        }
      }
    } else {
      lambda_max <- eig_res$lambda_max
      ritz_elapsed <- NA_real_
      ritz_timing <- system.time({
        rr <- tryCatch(
          ltsa_ritz_select(
            B = eig_res$matrix,
            vectors = eig_res$vectors,
            ndim = ndim,
            nullvec = nullvec,
            lambda_max = lambda_max
          ),
          error = function(e) e
        )
      })
      ritz_elapsed <- unname(ritz_timing[["elapsed"]])

      if (inherits(rr, "error")) {
        last_error <- conditionMessage(rr)
        attempts[[length(attempts) + 1L]] <- ltsa_attempt_error_summary(
          eig_k = eig_k,
          eig_res = eig_res,
          error = last_error,
          candidate_elapsed = candidate_elapsed,
          ritz_elapsed = ritz_elapsed,
          capture_candidate_space = TRUE
        )
        if (width_first_rescue && !is.null(width_rescue_incumbent)) {
          rescued <- finish_width_rescue_if_exhausted()
          if (!is.null(rescued)) {
            return(rescued)
          }
        }
      } else {
        acceptance <- accept_ltsa_ritz(
          rr = rr,
          ndim = ndim,
          lambda_max = lambda_max,
          resid_tol = resid_tol,
          gap_tol = gap_tol,
          require_gap = FALSE
        )
        attempts[[length(attempts) + 1L]] <- ltsa_attempt_summary(
          eig_k = eig_k,
          eig_res = eig_res,
          rr = rr,
          acceptance = acceptance,
          candidate_elapsed = candidate_elapsed,
          ritz_elapsed = ritz_elapsed,
          capture_candidate_space = TRUE
        )
        tsmessage(
          "Rayleigh-Ritz LTSA postprocess: post-null rank = ",
          rr$rank_after_null,
          ", max selected scaled residual = ",
          signif(acceptance$max_scaled_residual, 4),
          ", global boundary gap = ",
          signif(rr$global_gap, 4),
          ", local boundary gap = ",
          signif(rr$local_gap, 4),
          ", zero_tol = ",
          signif(rr$zero_tol, 4),
          ", near-zero nonconstant modes = ",
          rr$near_zero_nonconstant_count,
          ", threshold counts = [",
          ltsa_format_named_counts(rr$near_zero_nonconstant_counts),
          "], first Ritz values = [",
          ltsa_format_ritz_values(rr$reported_ritz_values),
          "]"
        )
        best <- ltsa_with_ritz(eig_res, rr, acceptance, attempts)

        if (width_first_rescue) {
          if (
            is.null(width_rescue_incumbent) &&
              ltsa_partial_block_rescue_allowed(best, ndim)
          ) {
            width_rescue_trigger_eig_k <- eig_k
            width_rescue_incumbent <- best
            attempts[[length(attempts)]]$width_first_rescue <- TRUE
            attempts[[length(attempts)]]$width_first_rescue_trigger <- TRUE
            attempts[[length(attempts)]]$width_first_rescue_decision <-
              "partial_incumbent"
            attempts[[length(attempts)]]$width_first_rescue_accept <- FALSE
          }

          if (!is.null(width_rescue_incumbent)) {
            attempts[[length(attempts)]]$width_first_rescue <- TRUE
            if (
              isTRUE(best$acceptance$rank_ok) &&
                isTRUE(best$acceptance$resid_ok)
            ) {
              if (ltsa_partial_near_zero_block(best, ndim)) {
                attempts[[length(attempts)]]$width_first_rescue_decision <-
                  "partial_incumbent"
                attempts[[length(attempts)]]$width_first_rescue_accept <-
                  FALSE
                replacement <- ltsa_width_rescue_partial_incumbent(
                  candidate = best,
                  incumbent = width_rescue_incumbent,
                  ndim = ndim
                )
                width_rescue_incumbent <- replacement
                width_rescue_incumbent$attempts <- attempts
              } else {
                decision <- ltsa_energy_rescue_accept_replacement(
                  incumbent = width_rescue_incumbent,
                  candidate = best,
                  ndim = ndim
                )
                attempts[[length(attempts)]] <-
                  ltsa_mark_width_rescue_attempt(
                    attempts[[length(attempts)]],
                    decision = decision
                  )
                best$attempts <- attempts
                if (isTRUE(decision$accept)) {
                  best$acceptance$return_reason <-
                    "width_first_ordinary_rescue"
                  best$acceptance$width_first_rescue_used <- TRUE
                  best$acceptance$width_first_rescue_ordinary <- TRUE
                  return(finish_width_first_rescue(best))
                }
                width_rescue_incumbent$attempts <- attempts
              }
            } else {
              width_rescue_incumbent$attempts <- attempts
            }

            rescued <- finish_width_rescue_if_exhausted()
            if (!is.null(rescued)) {
              return(rescued)
            }
          }
        }

        if (
          is.null(width_rescue_incumbent) &&
            acceptance$rank_ok &&
            acceptance$resid_ok &&
            acceptance$gap_ok
        ) {
          best$acceptance$return_reason <- "residual_rank_gap_ok"
          if (strict_rescue) {
            best <- ltsa_maybe_strict_rescue(
              B = B,
              selected = best,
              ndim = ndim,
              nullvec = nullvec,
              lambda_max = lambda_max,
              provider = provider,
              provider_args = provider_args,
              strict_rescue_arg_mapper = strict_rescue_arg_mapper,
              strict_rescue_controls = strict_rescue_controls,
              resid_tol = resid_tol,
              gap_tol = gap_tol,
              strict_rescue_extra = strict_rescue_extra,
              bank_backend = bank_backend,
              verbose = verbose
            )
          }
          best <- ltsa_maybe_warn_partial_near_zero_block(
            best,
            ndim = ndim,
            strict_rescue = strict_rescue
          )
          best <- ltsa_maybe_warn_spectral_ambiguity(best, ndim)
          return(ltsa_finalize_attempt_observability(
            best,
            reference_vectors = attempt_reference_vectors,
            nullvec = nullvec,
            retain_candidate_spaces = retain_attempt_candidate_spaces
          ))
        }
        if (is.null(width_rescue_incumbent) && acceptance$gap_issue_only) {
          if (is.null(weak_gap_candidate)) {
            first_weak_gap_eig_k <- eig_k
          }
          weak_gap_candidate <- ltsa_rescue_candidate(
            candidate = best,
            incumbent = weak_gap_candidate,
            ndim = ndim
          )
          weak_gap_candidate$acceptance$return_reason <-
            "weak_lowest_energy_residual_rank_good"
          weak_gap_candidate$acceptance$first_weak_gap_eig_k <-
            first_weak_gap_eig_k
          if (gap_expansions >= gap_expansion_steps || eig_k >= max_k) {
            selected <- weak_gap_candidate %||% best
            selected$attempts <- attempts
            selected$acceptance$return_reason <-
              "weak_lowest_energy_residual_rank_good"
            selected$acceptance$first_weak_gap_eig_k <- first_weak_gap_eig_k
            selected$acceptance$gap_expansions <- gap_expansions
            selected$acceptance$diagnostic_final_eig_k <- eig_k
            if (strict_rescue) {
              selected <- ltsa_maybe_strict_rescue(
                B = B,
                selected = selected,
                ndim = ndim,
                nullvec = nullvec,
                lambda_max = lambda_max,
                provider = provider,
                provider_args = provider_args,
                strict_rescue_arg_mapper = strict_rescue_arg_mapper,
                strict_rescue_controls = strict_rescue_controls,
                resid_tol = resid_tol,
                gap_tol = gap_tol,
                strict_rescue_extra = strict_rescue_extra,
                bank_backend = bank_backend,
                verbose = verbose
              )
            }
            if (isTRUE(selected$acceptance$strict_rescue_used)) {
              tsmessage(
                "LTSA strict rescue found lower-energy residual-good Ritz ",
                "vectors after an incomplete near-zero selected block; ",
                "returning strict-rescue result with selected values = [",
                ltsa_format_ritz_values(selected$values),
                "]"
              )
            } else {
              selected <- ltsa_maybe_warn_partial_near_zero_block(
                selected,
                ndim = ndim,
                strict_rescue = strict_rescue
              )
              selected <- ltsa_maybe_warn_spectral_ambiguity(selected, ndim)
              tsmessage(
                "LTSA Ritz boundary gap ",
                if (acceptance$gap_available) "remains below tolerance" else
                  "is unavailable",
                " after diagnostic expansion to ",
                eig_k,
                " candidates; returning residual-good, rank-good Ritz vectors."
              )
            }
            if (isTRUE(selected$acceptance$strict_rescue_used)) {
              selected <- ltsa_maybe_warn_partial_near_zero_block(
                selected,
                ndim = ndim,
                strict_rescue = strict_rescue
              )
              selected <- ltsa_maybe_warn_spectral_ambiguity(selected, ndim)
            }
            return(ltsa_finalize_attempt_observability(
              selected,
              reference_vectors = attempt_reference_vectors,
              nullvec = nullvec,
              retain_candidate_spaces = retain_attempt_candidate_spaces
            ))
          }

          gap_expansions <- gap_expansions + 1L
          attempts[[length(attempts)]]$weak_gap_expansion <- gap_expansions
        }
      }
    }

    if (eig_k >= max_k) {
      break
    }
    eig_k <- ltsa_next_ritz_candidate_k(eig_k, max_k)
  }

  if (is.null(best)) {
    stop(
      "LTSA Rayleigh-Ritz selection failed after requesting up to ",
      max_k,
      " candidate vectors: ",
      last_error %||% "no usable candidate subspace",
      call. = FALSE
    )
  }

  if (!best$acceptance$rank_ok) {
    stop(
      "LTSA Rayleigh-Ritz selection did not find enough post-null rank after ",
      "requesting up to ",
      best$eig_k,
      " candidate vectors",
      call. = FALSE
    )
  }
  if (strict_rescue) {
    best <- ltsa_maybe_strict_rescue(
      B = B,
      selected = best,
      ndim = ndim,
      nullvec = nullvec,
      lambda_max = best$lambda_max %||% lambda_max,
      provider = provider,
      provider_args = provider_args,
      strict_rescue_arg_mapper = strict_rescue_arg_mapper,
      strict_rescue_controls = strict_rescue_controls,
      resid_tol = resid_tol,
      gap_tol = gap_tol,
      strict_rescue_extra = strict_rescue_extra,
      bank_backend = bank_backend,
      verbose = verbose
    )
  }
  best <- ltsa_maybe_warn_partial_near_zero_block(
    best,
    ndim = ndim,
    strict_rescue = strict_rescue
  )
  best <- ltsa_maybe_warn_spectral_ambiguity(best, ndim)
  if (!best$acceptance$resid_ok) {
    warning(
      "LTSA Rayleigh-Ritz residuals remain above tolerance after requesting ",
      best$eig_k,
      " candidate vectors: max scaled residual = ",
      signif(best$acceptance$max_scaled_residual, 4),
      ", tolerance = ",
      signif(resid_tol, 4),
      call. = FALSE
    )
  } else if (!best$acceptance$gap_ok) {
    tsmessage(
      "LTSA Ritz boundary gap ",
      if (best$acceptance$gap_available) "remains below tolerance" else
        "is unavailable",
      " after requesting ",
      best$eig_k,
      " candidates: global gap = ",
      signif(best$acceptance$rel_gap, 4),
      ", local gap = ",
      signif(best$acceptance$local_gap, 4),
      ", tolerance = ",
      signif(gap_tol, 4),
      "; returning residual-good, rank-good Ritz vectors."
    )
  }

  ltsa_finalize_attempt_observability(
    best,
    reference_vectors = attempt_reference_vectors,
    nullvec = nullvec,
    retain_candidate_spaces = retain_attempt_candidate_spaces
  )
}

# Public-facing iterative wrappers. These are retained as fixed-width shims for
# internal callers; adaptive/rescue policy is no longer part of the runtime
# surface.
ltsa_rspectra_ritz_eig <- function(
  B,
  ndim,
  ...,
  nullvec = ltsa_default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE
) {
  provider_args <- list(...)
  ltsa_reject_rescue_policy_args(provider_args)

  ltsa_fixed_ritz_eig(
    B = B,
    ndim = ndim,
    provider = ltsa_rspectra_candidate_provider,
    provider_args = provider_args,
    nullvec = nullvec,
    eig_k = eig_k,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    verbose = verbose
  )
}

ltsa_irlba_ritz_eig <- function(
  B,
  ndim,
  ...,
  nullvec = ltsa_default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE
) {
  provider_args <- list(...)
  ltsa_reject_rescue_policy_args(provider_args)

  ltsa_fixed_ritz_eig(
    B = B,
    ndim = ndim,
    provider = ltsa_irlba_candidate_provider,
    provider_args = provider_args,
    nullvec = nullvec,
    eig_k = eig_k,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    verbose = verbose
  )
}

ltsa_svdr_ritz_eig <- function(
  B,
  ndim,
  ...,
  nullvec = ltsa_default_null_vector(nrow(B)),
  eig_k = NULL,
  resid_tol = 1e-5,
  gap_tol = 1e-4,
  verbose = FALSE
) {
  provider_args <- list(...)
  ltsa_reject_rescue_policy_args(provider_args)

  ltsa_fixed_ritz_eig(
    B = B,
    ndim = ndim,
    provider = ltsa_svdr_candidate_provider,
    provider_args = provider_args,
    nullvec = nullvec,
    eig_k = eig_k,
    resid_tol = resid_tol,
    gap_tol = gap_tol,
    verbose = verbose
  )
}
