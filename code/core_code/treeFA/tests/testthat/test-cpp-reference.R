test_that("prox_tree_cpp matches the pure R reference implementation", {
  tree_result <- gather_leaf_nodes_per_non_leaf(deep_tree())$list_result
  eta <- c(2.5, -0.2, 0.7, -1.4)
  lambda <- 0.37

  cpp <- treeFA:::prox_tree_cpp(eta, lambda = lambda, tree_list = tree_result)
  ref <- reference_prox_tree(eta, lambda = lambda, tree_result = tree_result)

  expect_equal(as.numeric(cpp), ref, tolerance = 1e-12)
})

test_that("C++ solver matches the pure R reference implementation", {
  set.seed(101)
  tree_df <- deep_tree()
  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)$list_result
  X <- matrix(rnorm(48), nrow = 12, ncol = 4)
  beta <- c(1.2, 1.2, -0.8, 0.4)
  Y <- as.numeric(X %*% beta + rnorm(12, sd = 0.2))

  ref <- reference_acc_core(
    Y = Y,
    X = X,
    tree_result = tree_result,
    lambda = 0.03,
    init_beta = rep(1, 4),
    ridge_param = 0.05,
    thresh = 1e-7,
    max_iter = 5000
  )
  cpp <- treeFA:::acc_prox_simple_linear_cpp(
    Y = Y,
    X = X,
    tree_result = tree_result,
    lambda = 0.03,
    init_beta = rep(1, 4),
    ridge_param = 0.05,
    thresh = 1e-7,
    max_iter = 5000L,
    stop_rule = 0L
  )

  expect_true(ref$converged)
  expect_true(cpp$converged)
  expect_equal(as.numeric(cpp$beta), ref$beta, tolerance = 1e-9)
})

test_that("R wrapper preserves warm starts and intercept reconstruction", {
  set.seed(202)
  tree_df <- deep_tree()
  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)$list_result
  X <- matrix(rnorm(60, mean = 0.4), nrow = 15, ncol = 4)
  beta <- c(0.6, 0.6, -0.2, 1.1)
  Y <- as.numeric(2.5 + X %*% beta + rnorm(15, sd = 0.1))
  init <- c(0.3, 0.3, -0.1, 0.8)

  fit <- acc_prox_simple_linear(
    Y,
    X,
    tree_result,
    lambda = 0.02,
    warm_start = TRUE,
    init_beta = init,
    intercept = TRUE,
    ridge_param = 0.01,
    thresh = 1e-7,
    max_iter = 5000
  )

  Xc <- scale(X, scale = FALSE)
  Yc <- Y - mean(Y)
  ref <- reference_acc_core(
    Y = Yc,
    X = Xc,
    tree_result = tree_result,
    lambda = 0.02,
    init_beta = init,
    ridge_param = 0.01,
    thresh = 1e-7,
    max_iter = 5000
  )

  expect_true(ref$converged)
  expect_equal(fit$beta1, ref$beta, tolerance = 1e-9)
  expect_equal(fit$beta0, as.numeric(mean(Y) - sum(colMeans(X) * ref$beta)), tolerance = 1e-9)
})
