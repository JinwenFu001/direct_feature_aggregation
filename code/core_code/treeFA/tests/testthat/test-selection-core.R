selection_two_leaf_tree <- function(weight = 1) {
  data.frame(node = 1:3, parent = c(3, 3, NA), weight = c(0, 0, weight))
}

selection_soft_reference <- function(x, level) sign(x) * pmax(abs(x) - level, 0)

test_that("selection reduces to the analytic orthogonal LASSO with ridge", {
  n <- 3
  X <- sqrt(n) * diag(n)
  signal <- c(2, -0.5, 0.1)
  Y <- as.numeric(X %*% signal)
  ridge <- 0.6
  tau <- 0.4
  fit <- dy_prox_simple_linear(
    Y, X, small_tree(), lambda = 0, tau = tau, ridge_param = ridge,
    abs_tol = 1e-10, rel_tol = 1e-10
  )
  expected <- selection_soft_reference(signal, tau) / (1 + ridge / n)
  expect_true(fit$converged)
  expect_equal(fit$beta, expected, tolerance = 1e-8)
  expect_identical(fit$beta[3], 0)
  expected_objective <- sum((Y - X %*% expected)^2) / (2 * n) +
    ridge * sum(expected^2) / (2 * n) + tau * sum(abs(expected))
  expect_equal(fit$objective, expected_objective, tolerance = 1e-10)
})

test_that("joint two-leaf fits agree with the analytic fused LASSO solution", {
  X <- sqrt(2) * diag(2)
  signal <- c(2, 0.4)
  for (lambda in c(0.6, 2)) {
    tau <- 0.5
    fit <- dy_prox_simple_linear(
      as.numeric(X %*% signal), X, selection_two_leaf_tree(),
      lambda = lambda, tau = tau, abs_tol = 1e-10, rel_tol = 1e-10
    )
    # For two leaves, Omega(beta) = abs(beta[1] - beta[2]) / sqrt(2).
    difference <- max(signal[1] - signal[2] - sqrt(2) * lambda, 0)
    fused <- mean(signal) + c(1, -1) * difference / 2
    expected <- selection_soft_reference(fused, tau)
    expect_true(fit$converged)
    expect_equal(fit$beta, expected, tolerance = 1e-8)
    objective <- sum((signal - fit$beta)^2) / 2 +
      lambda * abs(diff(fit$beta)) / sqrt(2) + tau * sum(abs(fit$beta))
    expect_equal(fit$objective, objective, tolerance = 1e-10)
  }
})

test_that("selection returns u and its matching state even at the iteration limit", {
  fit <- suppressWarnings(dy_prox_simple_linear(
    c(2, -1), diag(2), selection_two_leaf_tree(), lambda = 0.2, tau = 0.3,
    step_size = 0.7, init_z = c(0.9, -0.1), max_iter = 1,
    abs_tol = 1e-12, rel_tol = 1e-12, keep_history = TRUE
  ))
  expect_equal(fit$beta, selection_soft_reference(fit$state$z, 0.7 * 0.3),
               tolerance = 1e-14)
  expect_equal(fit$state$step_size, fit$step_size)
  expect_equal(fit$residual_scaled, fit$residual / fit$step_size,
               tolerance = 1e-12)
  expect_false(fit$converged)
  expect_true(is.finite(fit$objective))
})

test_that("a tiny step cannot report convergence from a tiny unscaled residual", {
  fit <- suppressWarnings(dy_prox_simple_linear(
    sqrt(2) * c(1, 2), sqrt(2) * diag(2), selection_two_leaf_tree(),
    lambda = 0, tau = 0, step_size = 1e-12, max_iter = 3,
    abs_tol = 1e-7, rel_tol = 1e-6
  ))
  expect_false(fit$converged)
  expect_lt(fit$residual, 1e-9)
  expect_gt(fit$residual_scaled, 1)
  expect_gt(fit$objective, 2)
})

test_that("selection accepts legacy tree input and resumes a saved z state", {
  set.seed(914)
  X <- matrix(rnorm(60), 20, 3)
  Y <- as.numeric(1.3 + X %*% c(1, 1, -0.4) + rnorm(20, sd = 0.1))
  tree <- small_tree()
  tree_list <- gather_leaf_nodes_per_non_leaf(tree)$list_result
  cold <- dy_prox_simple_linear(
    Y, X, tree, lambda = 0.12, tau = 0.06, intercept = TRUE,
    abs_tol = 1e-10, rel_tol = 1e-10
  )
  early <- suppressWarnings(dy_prox_simple_linear(
    Y, X, tree_list, lambda = 0.12, tau = 0.06, intercept = TRUE,
    max_iter = 3, abs_tol = 1e-10, rel_tol = 1e-10
  ))
  resumed <- dy_prox_simple_linear(
    Y, X, tree_list, lambda = 0.12, tau = 0.06, intercept = TRUE,
    init_z = early$state$z, step_size = early$step_size,
    abs_tol = 1e-10, rel_tol = 1e-10
  )
  expect_true(cold$converged)
  expect_true(resumed$converged)
  expect_equal(resumed$beta, cold$beta, tolerance = 1e-7)
  expect_equal(resumed$beta0, cold$beta0, tolerance = 1e-7)
  expect_equal(cold$beta0, mean(Y) - sum(colMeans(X) * cold$beta),
               tolerance = 1e-12)
})

test_that("zero designs and inactive penalties have well-defined solutions", {
  X <- matrix(0, 6, 2)
  Y <- c(-2, -1, 0, 1, 2, 3)
  fit <- dy_prox_simple_linear(Y, X, selection_two_leaf_tree(0),
                               lambda = 0, tau = 0, intercept = TRUE)
  expect_true(fit$converged)
  expect_equal(fit$beta, c(0, 0))
  expect_equal(fit$beta0, mean(Y))
  expect_true(is.finite(fit$step_size) && fit$step_size > 0)
})

test_that("overflow when reconstructing an intercept is not reported as convergence", {
  tree <- data.frame(node = 1L, parent = NA, weight = 0)
  fit <- suppressWarnings(dy_prox_simple_linear(
    c(0, 0), matrix(1e200, 2, 1), tree, lambda = 0, tau = 0,
    intercept = TRUE, init_z = 1e109
  ))
  expect_false(fit$converged)
  expect_identical(fit$status, "nonfinite")
  expect_false(is.finite(fit$beta0))
})

test_that("logistic intercept is fitted without centering binary responses", {
  Y <- c(rep(0, 7), rep(1, 3))
  X <- matrix(0, length(Y), 2)
  fit <- dy_prox_simple_logistic(
    Y, X, selection_two_leaf_tree(), lambda = 0.3, tau = 0.2,
    intercept = TRUE, abs_tol = 1e-10, rel_tol = 1e-10
  )
  expect_true(fit$converged)
  expect_equal(fit$beta, c(0, 0), tolerance = 1e-10)
  expect_equal(fit$beta0, qlogis(mean(Y)), tolerance = 1e-7)
  expect_equal(fit$objective,
    -mean(Y) * log(mean(Y)) - (1 - mean(Y)) * log1p(-mean(Y)),
    tolerance = 1e-10)

  for (response in list(rep(0, 10), rep(1, 10))) {
    expect_error(dy_prox_simple_logistic(
      response, X, selection_two_leaf_tree(), lambda = 0.3, tau = 0.2,
      intercept = TRUE
    ))
  }
})

test_that("unpenalized logistic fits agree with independent IRLS when finite", {
  x <- rep(-2:2, each = 8)
  X <- cbind(x, x^2 - mean(x^2))
  Y <- unlist(lapply(c(2, 3, 3, 5, 6), function(k) c(rep(1, k), rep(0, 8 - k))))
  reference <- stats::glm.fit(cbind(1, X), Y, family = stats::binomial())
  fit <- dy_prox_simple_logistic(
    Y, X, selection_two_leaf_tree(), lambda = 0, tau = 0, intercept = TRUE,
    abs_tol = 1e-10, rel_tol = 1e-10
  )
  expect_true(reference$converged)
  expect_true(fit$converged)
  expect_equal(as.numeric(c(fit$beta0, fit$beta)), unname(reference$coefficients),
               tolerance = 2e-6)
})

test_that("selection validates finite tuning parameters and valid step sizes", {
  fit <- function(...) dy_prox_simple_linear(
    c(1, -1), sqrt(2) * diag(2), selection_two_leaf_tree(), ...
  )
  for (bad in c(-1, NA_real_, Inf)) {
    expect_error(fit(lambda = bad, tau = 0))
    expect_error(fit(lambda = 0, tau = bad))
  }
  for (bad in c(-1, 0, 2, NA_real_, Inf)) {
    expect_error(fit(lambda = 0, tau = 0, step_size = bad))
  }
  expect_error(fit(lambda = 0, tau = 0, init_z = 1))
  expect_error(fit(lambda = 0, tau = 0, init_z = c(0, Inf)))
})

test_that("selection rejects malformed trees before entering compiled code", {
  base <- selection_two_leaf_tree()
  wrong_id <- base
  wrong_id$node[2] <- wrong_id$node[1]
  cycle <- base
  cycle$parent[3] <- 1
  negative <- base
  negative$weight[3] <- -1
  nonfinite <- base
  nonfinite$weight[3] <- Inf
  for (tree in list(wrong_id, cycle, negative, nonfinite)) {
    expect_error(dy_prox_simple_linear(c(1, 2), diag(2), tree, lambda = 1, tau = 1))
  }
  malformed <- list(groups = list(list(c(1L, 3L))), weights = list(1))
  expect_error(dy_prox_simple_linear(c(1, 2), diag(2), malformed, lambda = 1, tau = 1))
})
