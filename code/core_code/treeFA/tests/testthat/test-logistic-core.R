test_that("negative logistic likelihood handles vector and matrix beta", {
  set.seed(606)
  X <- matrix(rnorm(24), nrow = 8, ncol = 3)
  beta <- c(0.6, -0.4, 0.2)
  eta <- as.numeric(X %*% beta)
  Y <- rbinom(8, size = 1, prob = 1 / (1 + exp(-eta)))

  vec_loss <- negtv_lglkh(Y, X, beta)
  mat_loss <- negtv_lglkh(Y, X, cbind(beta, 0.5 * beta))

  expect_length(vec_loss, 1)
  expect_length(mat_loss, 2)
  expect_true(all(is.finite(vec_loss)))
  expect_true(all(is.finite(mat_loss)))
  expect_equal(vec_loss, mat_loss[1])
})

test_that("logistic solver and grid return expected shapes", {
  set.seed(707)
  tree_df <- small_tree()
  X <- matrix(rnorm(45), nrow = 15, ncol = 3)
  true_beta <- c(1, 1, -0.8)
  prob <- 1 / (1 + exp(-as.numeric(X %*% true_beta)))
  Y <- rbinom(15, size = 1, prob = prob)

  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)$list_result
  fit <- acc_prox_simple_logistic(Y, X, tree_result, lambda = 0.01, max_iter = 1000)
  expect_type(fit, "double")
  expect_length(fit, 3)
  expect_true(all(is.finite(fit)))

  path <- grid.logistic(Y, X, tree_df, true_beta = true_beta, seqc = c(0.02, 0.01), max_iter = 1000)
  expect_equal(dim(path$beta), c(3L, 2L))
  expect_length(path$loss, 2)
  expect_equal(path$lambda, c(0.02, 0.01))
  expect_true(all(is.finite(negtv_lglkh(Y, X, path$beta))))
})

test_that("cv.logistic returns a valid selected parameter and fitted beta", {
  set.seed(808)
  tree_df <- deep_tree()
  X <- matrix(rnorm(96), nrow = 24, ncol = 4)
  true_beta <- c(0.7, 0.7, -0.5, 0.2)
  prob <- 1 / (1 + exp(-as.numeric(X %*% true_beta)))
  Y <- rbinom(24, size = 1, prob = prob)
  seqc <- c(0.05, 0.02)

  fit <- cv.logistic(
    Y,
    X,
    tree_df,
    folds = 3,
    seqc = seqc,
    thresh = 1e-5,
    max_iter = 1000
  )

  expect_true(fit$selected.param %in% seqc)
  expect_length(fit$beta, 4)
  expect_length(fit$valid.error, length(seqc))
  expect_true(all(is.finite(fit$valid.error)))
})

test_that("logistic input validation catches invalid calls", {
  tree_df <- small_tree()
  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)$list_result
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  Y <- rbinom(10, size = 1, prob = 0.5)

  expect_error(acc_prox_simple_logistic(Y + 0.1, X, tree_result, lambda = 0.1), "binary")
  expect_error(acc_prox_simple_logistic(Y, X, tree_result, lambda = -0.1), "lambda")
  expect_error(grid.logistic(Y, X, tree_df, true_beta = c(1, 2)), "true_beta")
  expect_error(negtv_lglkh(Y, X, matrix(0, nrow = 2, ncol = 1)), "beta")
})
