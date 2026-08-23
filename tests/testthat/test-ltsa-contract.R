test_that("LTSA orchestration contract is stable across public route choices", {
  angle <- 2 * pi * seq.int(0L, 27L) / 28L
  X <- cbind(
    cos(angle),
    sin(angle),
    cos(2 * angle),
    sin(2 * angle),
    cos(3 * angle),
    sin(3 * angle)
  )
  rownames(X) <- paste0("sample-", seq_len(nrow(X)))

  for (include_self in c(TRUE, FALSE)) {
    nn_idx <- flotsam:::prepare_neighbors(
      X = X,
      n_neighbors = 7L,
      nn_method = "exact",
      nn_idx = NULL,
      include_self = include_self,
      n_threads = 1L
    )$nn_idx

    for (normalize in c(FALSE, TRUE)) {
      canonical <- NULL
      for (source in c("computed", "precomputed")) {
        for (n_assembly_threads in c(1L, 2L)) {
          nn_method <- if (source == "computed") "exact" else nn_idx
          common_args <- list(
            X = X,
            n_neighbors = 7L,
            ndim = 2L,
            nn_method = nn_method,
            include_self = include_self,
            normalize = normalize,
            n_threads = 1L,
            n_assembly_threads = n_assembly_threads
          )

          B <- expect_silent(
            do.call(ltsa, c(common_args, list(output = "B")))
          )
          embedding <- expect_silent(do.call(
            ltsa,
            c(
              common_args,
              list(eig_method = "eig", eig_k = 6L, output = "embedding")
            )
          ))
          result <- expect_silent(do.call(
            ltsa,
            c(
              common_args,
              list(
                eig_method = "eig",
                eig_k = 6L,
                output = "result",
                include_B = TRUE
              )
            )
          ))

          expect_s4_class(B, "dgCMatrix")
          expect_identical(dim(B), c(28L, 28L))
          expect_identical(
            dimnames(embedding),
            list(rownames(X), c("LTSA1", "LTSA2"))
          )
          expect_named(result, c("embedding", "eigen", "assembly", "B"))
          expect_same_subspace(result$embedding, embedding, tolerance = 1e-12)
          expect_identical(result$eigen$normalized, normalize)
          expect_identical(result$assembly$n_neighbors, 7L)
          expect_identical(result$assembly$include_self, include_self)
          expect_identical(
            result$assembly$neighbor_source,
            if (source == "computed") "exact" else "precomputed"
          )
          expect_identical(
            result$assembly$assembly_route,
            if (n_assembly_threads == 1L) {
              "serial_triangular"
            } else {
              "parallel_triangular_two_pass"
            }
          )
          expect_identical(
            result$assembly$requested_assembly_threads,
            n_assembly_threads
          )

          if (normalize) {
            normalization <- flotsam:::normalize_sparse_operator(result$B)
            expect_sparse_equivalent(
              B,
              normalization$normalized_operator,
              tolerance = 0
            )
          } else {
            expect_sparse_equivalent(B, result$B, tolerance = 0)
          }

          if (source == "computed" && n_assembly_threads == 1L) {
            canonical <- list(B = B, embedding = embedding, result = result)
          } else {
            expect_sparse_equivalent(B, canonical$B, tolerance = 1e-11)
            expect_sparse_equivalent(
              result$B,
              canonical$result$B,
              tolerance = 1e-11
            )
            expect_same_subspace(
              embedding,
              canonical$embedding,
              tolerance = 1e-10
            )
            expect_equal(
              result$eigen$values,
              canonical$result$eigen$values,
              tolerance = 1e-10
            )
            stable_assembly_fields <- c(
              "n_neighbors",
              "include_self",
              "rank_deficient_count",
              "min_local_rank",
              "component_count",
              "component_sizes",
              "component_membership"
            )
            expect_identical(
              result$assembly[stable_assembly_fields],
              canonical$result$assembly[stable_assembly_fields]
            )
          }
        }
      }
    }
  }
})
