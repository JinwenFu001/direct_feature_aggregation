test_that("selection grids preserve candidate order, duplicates, and coefficient mapping", {
  set.seed(207)
  X <- matrix(rnorm(72), 24, 3)
  Y <- as.numeric(0.4 + X %*% c(1, 1, -0.5) + rnorm(24, sd = 0.1))
  candidates <- data.frame(lambda = c(0.2, 0, 0.2, 0.07), tau = c(0.05, 0.15, 0.05, 0))
  path <- grid.select_linear(
    Y, X, small_tree(), param_grid = candidates, intercept = TRUE,
    abs_tol = 1e-9, rel_tol = 1e-9, warm_start = TRUE
  )
  expect_equal(path$param_grid, candidates, ignore_attr = TRUE)
  expect_equal(dim(path$beta), c(3L, 4L))
  expect_length(path$beta0, 4)
  expect_true(all(path$diagnostics$converged))
  expect_equal(path$n_nonzero, as.integer(colSums(path$beta != 0)))
  for (j in seq_len(nrow(candidates))) {
    single <- dy_prox_simple_linear(
      Y, X, small_tree(), candidates$lambda[j], candidates$tau[j], intercept = TRUE,
      abs_tol = 1e-9, rel_tol = 1e-9
    )
    expect_equal(path$beta[, j], single$beta, tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(path$beta0[j], single$beta0, tolerance = 1e-6, ignore_attr = TRUE)
  }
  expect_equal(path$beta[, 1], path$beta[, 3], tolerance = 1e-7)
})

test_that("explicit sequences form the requested Cartesian selection grid", {
  X <- diag(3)
  path <- grid.select_linear(
    c(2, 1, -1), X, small_tree(), lambda_seq = c(0.2, 0), tau_seq = c(0.1, 0.3)
  )
  expect_equal(path$param_grid, expand.grid(lambda = c(0.2, 0), tau = c(0.1, 0.3)),
               ignore_attr = TRUE)
  expect_error(grid.select_linear(c(2, 1, -1), X, small_tree(), lambda_seq = c(0.2, 0)))
  for (grid in list(
    data.frame(lambda = -1, tau = 0),
    data.frame(lambda = 0, tau = Inf),
    data.frame(lambda = "bad", tau = 0),
    data.frame(lambda = numeric(0), tau = numeric(0))
  )) {
    expect_error(grid.select_linear(c(2, 1, -1), X, small_tree(), param_grid = grid))
  }
})

test_that("selection CV uses the specified folds and refits the selected pair", {
  set.seed(619)
  X <- matrix(rnorm(36), 12, 3)
  Y <- as.numeric(0.6 + X %*% c(0.7, 0.7, -0.4) + rnorm(12, sd = 0.2))
  foldid <- rep(1:3, c(3, 4, 5))
  candidates <- data.frame(lambda = c(0, 0.15), tau = c(0.1, 0.04))
  result <- cv.select_linear(
    Y, X, small_tree(), param_grid = candidates, foldid = foldid, intercept = TRUE,
    abs_tol = 1e-9, rel_tol = 1e-9
  )
  expect_equal(result$foldid, foldid, ignore_attr = TRUE)
  expect_equal(result$param_grid, candidates, ignore_attr = TRUE)
  expect_true(all(result$eligible))
  manual <- matrix(NA_real_, 3, nrow(candidates))
  for (fold in 1:3) {
    train <- foldid != fold
    for (j in seq_len(nrow(candidates))) {
      fit <- dy_prox_simple_linear(
        Y[train], X[train, , drop = FALSE], small_tree(),
        candidates$lambda[j], candidates$tau[j], intercept = TRUE,
        abs_tol = 1e-9, rel_tol = 1e-9
      )
      prediction <- fit$beta0 + X[!train, , drop = FALSE] %*% fit$beta
      manual[fold, j] <- mean((Y[!train] - prediction)^2)
    }
  }
  expect_equal(result$fold.error, manual, tolerance = 2e-6, ignore_attr = TRUE)
  expected <- colSums(manual * as.numeric(table(foldid))) / length(Y)
  expect_equal(result$cv.error, expected, tolerance = 2e-6, ignore_attr = TRUE)
  expect_equal(result$selected.index, which.min(expected))
  selected <- candidates[result$selected.index, ]
  expect_equal(result$selected.lambda, selected$lambda)
  expect_equal(result$selected.tau, selected$tau)
  reference <- dy_prox_simple_linear(
    Y, X, small_tree(), selected$lambda, selected$tau, intercept = TRUE,
    abs_tol = 1e-9, rel_tol = 1e-9
  )
  expect_true(result$fit$converged)
  expect_equal(result$beta, reference$beta, tolerance = 1e-6)
  expect_equal(result$beta0, reference$beta0, tolerance = 1e-6)
})

test_that("selection CV rejects unavailable candidates and malformed fold IDs", {
  X <- rbind(diag(3), 2 * diag(3), -diag(3))
  Y <- as.numeric(X %*% c(1, -1, 0.5))
  grid <- data.frame(lambda = 0, tau = 0)
  expect_error(suppressWarnings(cv.select_linear(
    Y, X, small_tree(), param_grid = grid, foldid = rep(1:3, each = 3),
    step_size = 1e-12, max_iter = 1, retry_max_iter = 0
  )))
  expect_error(cv.select_linear(Y, X, small_tree(), param_grid = grid, foldid = 1:2))
  expect_error(cv.select_linear(Y, X, small_tree(), param_grid = grid, foldid = rep(1, 9)))
  expect_error(cv.select_linear(Y, X, small_tree(), param_grid = grid,
                                foldid = c(NA, rep(1:2, 4))))
})

test_that("logistic selection paths and CV retain the unpenalized intercept", {
  X <- matrix(0, 12, 2)
  Y <- c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1)
  foldid <- rep(1:3, each = 4)
  tree <- data.frame(node = 1:3, parent = c(3, 3, NA), weight = c(0, 0, 1))
  candidates <- data.frame(lambda = c(0.1, 0.3), tau = c(0.1, 0.4))
  path <- grid.select_logistic(
    Y, X, tree, param_grid = candidates, intercept = TRUE,
    abs_tol = 1e-10, rel_tol = 1e-10
  )
  expect_true(all(path$diagnostics$converged))
  expect_equal(as.numeric(path$beta), rep(0, 4))
  expect_equal(as.numeric(path$beta0), rep(qlogis(mean(Y)), 2), tolerance = 1e-7)
  result <- cv.select_logistic(
    Y, X, tree, param_grid = candidates, foldid = foldid,
    intercept = TRUE, abs_tol = 1e-10, rel_tol = 1e-10
  )
  expected <- vapply(1:3, function(fold) {
    probability <- mean(Y[foldid != fold])
    -mean(Y[foldid == fold] * log(probability) +
            (1 - Y[foldid == fold]) * log1p(-probability))
  }, numeric(1))
  expect_equal(result$fold.error, cbind(expected, expected),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(result$beta0, qlogis(mean(Y)), tolerance = 1e-7)
  expect_true(result$fit$converged)
})
