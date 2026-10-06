test_that("tree helpers collect leaves and non-leaf groups", {
  tree_df <- small_tree()
  expect_equal(tree_leaf_ids(tree_df), 1:3)
  expect_equal(gather_leaves(5, tree_df), 1:3)
  expect_equal(gather_direct_latent_nodes(5, tree_df), 4)

  info <- gather_leaf_nodes_per_non_leaf(tree_df)
  expect_true(all(c("node", "depth", "leaves", "latent_children", "weight") %in% names(info$df_result)))
  expect_equal(info$df_result$node[info$df_result$depth == 0], 5)
})

test_that("find_coarest handles root-zero trees", {
  tree_df <- small_tree(root_weight = 0)
  info <- gather_leaf_nodes_per_non_leaf(tree_df)
  coarest <- find_coarest(tree_df, info$df_result)

  expect_true(all(seq_len(3) %in% sort(unique(unlist(coarest$leaves, use.names = FALSE)))))
  expect_true(any(coarest$weight == 1))
})

test_that("linear solver and grid return expected shapes", {
  set.seed(1)
  tree_df <- small_tree()
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  beta <- c(1, 1, -0.5)
  Y <- as.numeric(X %*% beta + rnorm(10, sd = 0.1))

  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)$list_result
  fit <- acc_prox_simple_linear(Y, X, tree_result, lambda = 0.01, max_iter = 1000)
  expect_type(fit, "double")
  expect_length(fit, 3)
  expect_true(all(is.finite(fit)))

  path <- grid.simple_linear(Y, X, tree_df, true_beta = beta, seqc = c(0.02, 0.01), max_iter = 1000)
  expect_equal(dim(path$beta), c(3L, 2L))
  expect_length(path$loss, 2)
  expect_equal(path$lambda, c(0.02, 0.01))
})
