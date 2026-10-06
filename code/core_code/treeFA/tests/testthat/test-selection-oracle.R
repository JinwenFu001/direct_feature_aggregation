# An independent consensus ADMM oracle. It explicitly splits every centered
# group operator and solves a linear system; it does not call either tree prox
# or the Davis--Yin implementation under test.
selection_admm_oracle <- function(Y, X, groups, weights, lambda, tau, ridge = 0) {
  n <- nrow(X)
  p <- ncol(X)
  operators <- lapply(groups, function(ind) {
    m <- length(ind)
    (diag(m) - matrix(1 / m, m, m)) %*% diag(p)[ind, , drop = FALSE]
  })
  B <- do.call(rbind, c(list(diag(p)), operators))
  ranges <- split(seq.int(p + 1L, nrow(B)), rep(seq_along(groups), lengths(groups)))
  factor <- chol(crossprod(X) / n + diag(ridge / n, p) + crossprod(B))
  rhs_loss <- as.numeric(crossprod(X, Y)) / n
  split_value <- dual <- numeric(nrow(B))
  beta <- numeric(p)
  residual <- Inf
  for (iteration in seq_len(20000L)) {
    rhs <- rhs_loss + as.numeric(crossprod(B, split_value - dual))
    beta <- as.numeric(backsolve(factor, forwardsolve(t(factor), rhs)))
    Bbeta <- as.numeric(B %*% beta)
    value <- Bbeta + dual
    previous <- split_value
    split_value[seq_len(p)] <- sign(value[seq_len(p)]) * pmax(abs(value[seq_len(p)]) - tau, 0)
    for (j in seq_along(groups)) {
      block <- ranges[[j]]
      norm <- sqrt(sum(value[block]^2))
      split_value[block] <- if (norm == 0) 0 else value[block] * max(1 - lambda * weights[j] / norm, 0)
    }
    dual <- dual + Bbeta - split_value
    residual <- max(sqrt(sum((Bbeta - split_value)^2)),
                    sqrt(sum(crossprod(B, split_value - previous)^2)))
    if (residual < 1e-11) break
  }
  list(beta = beta, residual = residual)
}

test_that("nested joint Gaussian estimates agree with an independent convex oracle", {
  set.seed(808)
  X <- matrix(rnorm(64), 16, 4)
  X <- sweep(X, 2, c(0.7, 1.1, 1.4, 0.9), "*")
  Y <- as.numeric(X %*% c(1, 1, 0, -0.7) + rnorm(16, sd = 0.15))
  groups <- list(1:2, 3:4, 1:4)
  weights <- c(1, 0.8, 0.5)
  for (penalties in list(c(0.15, 0.09), c(0.2, 0), c(0, 0.15))) {
    lambda <- penalties[1]
    tau <- penalties[2]
    oracle <- selection_admm_oracle(Y, X, groups, weights, lambda, tau, ridge = 0.2)
    fit <- dy_prox_simple_linear(
      Y, X, deep_tree(), lambda, tau, ridge_param = 0.2,
      abs_tol = 1e-10, rel_tol = 1e-10
    )
    expect_lt(oracle$residual, 1e-9)
    expect_true(fit$converged)
    expect_equal(fit$beta, oracle$beta, tolerance = 2e-7)
    tree_value <- sum(vapply(seq_along(groups), function(j) {
      block <- fit$beta[groups[[j]]]
      weights[j] * sqrt(sum((block - mean(block))^2))
    }, numeric(1)))
    objective <- mean((Y - X %*% fit$beta)^2) / 2 +
      0.2 * sum(fit$beta^2) / (2 * nrow(X)) + lambda * tree_value + tau * sum(abs(fit$beta))
    expect_equal(fit$objective, objective, tolerance = 1e-10)
  }
})

test_that("lambda zero matches glmnet with its standardization disabled", {
  skip_if_not_installed("glmnet")
  set.seed(300)
  X <- matrix(rnorm(120), 40, 3)
  Y <- as.numeric(0.8 + X %*% c(1.4, 0, -0.7) + rnorm(40, sd = 0.2))
  tau <- 0.13
  reference <- glmnet::glmnet(X, Y, family = "gaussian", alpha = 1, lambda = tau,
                              intercept = TRUE, standardize = FALSE, thresh = 1e-14)
  fit <- dy_prox_simple_linear(Y, X, small_tree(), lambda = 0, tau = tau,
                               intercept = TRUE, abs_tol = 1e-10, rel_tol = 1e-10)
  expect_true(fit$converged)
  expect_equal(c(fit$beta0, fit$beta), as.numeric(stats::coef(reference)), tolerance = 1e-6)
})

test_that("logistic LASSO plus ridge matches glmnet with the same objective", {
  skip_if_not_installed("glmnet")
  set.seed(316)
  X <- matrix(rnorm(180), 60, 3)
  Y <- rbinom(60, 1, plogis(0.7 + X %*% c(1.1, 0, -0.8)))
  tau <- 0.08
  ridge <- 0.6
  # glmnet uses lambda * ((1 - alpha) * ||beta||^2 / 2 + alpha * ||beta||_1).
  glmnet_lambda <- tau + ridge / nrow(X)
  reference <- glmnet::glmnet(
    X, Y, family = "binomial", lambda = glmnet_lambda,
    alpha = tau / glmnet_lambda, intercept = TRUE, standardize = FALSE,
    thresh = 1e-14
  )
  fit <- dy_prox_simple_logistic(
    Y, X, small_tree(), lambda = 0, tau = tau, ridge_param = ridge,
    intercept = TRUE, abs_tol = 1e-10, rel_tol = 1e-10
  )
  expect_true(fit$converged)
  expect_equal(c(fit$beta0, fit$beta), as.numeric(stats::coef(reference)), tolerance = 2e-6)
})
