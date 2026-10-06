test_that("find_max_param.linear agrees with an independent reference calculation", {
  set.seed(303)
  tree_df <- deep_tree()
  info <- gather_leaf_nodes_per_non_leaf(tree_df)
  coarest <- find_coarest(tree_df, info$df_result)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  Y <- rnorm(10)

  expect_equal(
    find_max_param.linear(Y, X, coarest),
    reference_max_param(Y, X, coarest),
    tolerance = 1e-12
  )
})

test_that("grid.simple_linear matches explicit sequential solver calls", {
  set.seed(404)
  tree_df <- deep_tree()
  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)$list_result
  X <- matrix(rnorm(48), nrow = 12, ncol = 4)
  true_beta <- c(1, 1, -0.5, 0.2)
  Y <- as.numeric(X %*% true_beta + rnorm(12, sd = 0.15))
  seqc <- c(0.04, 0.02, 0.01)

  path <- grid.simple_linear(
    Y,
    X,
    tree_df,
    true_beta = true_beta,
    seqc = seqc,
    thresh = 1e-6,
    max_iter = 5000
  )

  manual <- matrix(0, nrow = 4, ncol = length(seqc))
  manual[, 1] <- acc_prox_simple_linear(
    Y,
    X,
    tree_result,
    lambda = seqc[1],
    thresh = 1e-6,
    max_iter = 5000
  )
  for (i in 2:length(seqc)) {
    manual[, i] <- acc_prox_simple_linear(
      Y,
      X,
      tree_result,
      lambda = seqc[i],
      warm_start = TRUE,
      init_beta = manual[, i - 1],
      thresh = 1e-6,
      max_iter = 5000
    )
  }
  manual_loss <- colSums((matrix(rep(true_beta, length(seqc)), nrow = 4) - manual)^2) / 4

  expect_equal(path$beta, manual, tolerance = 1e-12)
  expect_equal(path$loss, manual_loss, tolerance = 1e-12)
})

test_that("cv.simple_linear returns a valid selected parameter and fitted beta", {
  set.seed(505)
  tree_df <- deep_tree()
  X <- matrix(rnorm(72), nrow = 18, ncol = 4)
  Y <- as.numeric(X %*% c(0.5, 0.5, -0.3, 0.1) + rnorm(18, sd = 0.2))
  seqc <- c(0.05, 0.02)

  fit <- cv.simple_linear(
    Y,
    X,
    tree_df,
    folds = 3,
    seqc = seqc,
    thresh = 1e-5,
    max_iter = 3000
  )

  expect_true(fit$selected.param %in% seqc)
  expect_length(fit$beta, 4)
  expect_length(fit$valid.error, length(seqc))
  expect_true(all(is.finite(fit$valid.error)))
})

test_that("input validation catches common invalid calls", {
  tree_df <- small_tree()
  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)$list_result
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  Y <- rnorm(10)

  expect_error(gather_leaves(99, tree_df), "node is not present")
  expect_error(acc_prox_simple_linear(Y, X, tree_result, lambda = -0.1), "lambda")
  expect_error(acc_prox_simple_linear(Y, X, tree_result, lambda = 0.1, ridge_param = -1), "ridge_param")
  expect_error(grid.simple_linear(Y, X, tree_df, true_beta = c(1, 2)), "true_beta")

  bad_tree <- tree_df
  bad_tree$parent[1] <- 99
  expect_error(gather_leaf_nodes_per_non_leaf(bad_tree), "Every non-NA parent")
})
