test_that("iris", {
  iris_nn50_embedding <-
    ltsa(
      iris,
      nn_method = "exact",
      n_neighbors = 50,
      include_self = FALSE,
      eig_method = "eig"
    )
  expect_same_subspace(
    iris_nn50_embedding,
    iris_nn50_embedding_expected,
    tolerance = 1e-2
  )
})
