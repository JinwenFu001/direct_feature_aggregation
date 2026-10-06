.qr_residual <- function(X, Y) {
  qr.resid(qr(X), Y)
}

find_max_param.linear <- function(Y, X, coarest_set) {
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X

  n <- nrow(X)
  p <- ncol(X)
  all_leaves <- seq_len(p)

  .check_required_columns(coarest_set, c("leaves", "weight"), "coarest_set")
  valid_set <- coarest_set[!is.na(coarest_set$weight) & coarest_set$weight != 0, , drop = FALSE]

  if (!nrow(valid_set)) {
    r <- .qr_residual(X, Y)
    g <- as.numeric(crossprod(X, r))
    return(max(abs(g)) / n)
  }

  if (any(valid_set$weight < 0)) {
    stop("coarest_set weights must be non-negative.", call. = FALSE)
  }

  grp_idx <- valid_set$leaves
  covered <- sort(unique(unlist(grp_idx, use.names = FALSE)))
  covered <- covered[!is.na(covered)]

  if (!all(covered %in% all_leaves)) {
    stop("coarest_set contains leaf indices outside 1:p.", call. = FALSE)
  }

  total_group_size <- sum(vapply(grp_idx, function(v) length(unique(v)), integer(1)))
  if (length(covered) != total_group_size) {
    stop("Group leaves appear to overlap; Q construction assumes disjoint groups.", call. = FALSE)
  }

  p1 <- nrow(valid_set) + (p - length(covered))
  Q <- matrix(0, nrow = p, ncol = p1)

  for (i in seq_len(nrow(valid_set))) {
    ind <- match(valid_set$leaves[[i]], all_leaves)
    Q[ind, i] <- 1
  }

  singletons <- setdiff(all_leaves, covered)
  if (length(singletons)) {
    cols <- seq_len(length(singletons)) + nrow(valid_set)
    Q[cbind(singletons, cols)] <- 1
  }

  X1 <- X %*% Q
  r <- .qr_residual(X1, Y)
  g <- as.numeric(crossprod(X, r))
  g2 <- g^2
  agg <- as.numeric(t(Q) %*% g2)[seq_len(nrow(valid_set))]
  vals <- sqrt(agg) / valid_set$weight

  max(vals) / n
}

.find_coarest_legacy <- function(df, result) {
  all_cover <- FALSE
  max_depth <- max(result$depth)
  coarest_set <- result[result$depth == 0, ]

  if (result[result$depth == 0, ]$weight != 0) {
    return(coarest_set)
  }

  current_depth <- 1
  while (current_depth <= max_depth && all_cover == FALSE) {
    sub_result <- result[result$depth == current_depth, ]
    sub_result <- sub_result[!is.element(df$parent[sub_result$node], coarest_set$node[-1]), ]

    if (nrow(sub_result) == 0) {
      if (coarest_set$weight[1] == 0) coarest_set <- coarest_set[-1, ]
      return(coarest_set)
    }

    coarest_set <- rbind(coarest_set, sub_result[sub_result$weight != 0, ])
    if (length(unique(unlist(coarest_set[-1, ]$leaves))) == length(coarest_set[1, ]$leaves[[1]])) {
      all_cover <- TRUE
    }
    current_depth <- current_depth + 1
  }

  if (coarest_set$weight[1] == 0) coarest_set <- coarest_set[-1, ]
  coarest_set
}

.find_max_param_linear_legacy <- function(Y, X, coarest_set) {
  p <- ncol(X)
  n <- nrow(X)
  leaves <- sort(coarest_set[1, ]$leaves[[1]], decreasing = FALSE)
  valid_set <- coarest_set[coarest_set$weight != 0, ]
  p1 <- p - sum(unlist(lapply(valid_set$leaves, length))) + nrow(valid_set)
  Q <- matrix(0, nrow = p, ncol = p1)
  remain <- rep(1, length(leaves))

  for (i in seq_len(nrow(valid_set))) {
    ind <- match(valid_set[i, ]$leaves[[1]], leaves)
    Q[ind, i] <- 1
    remain[ind] <- 0
  }

  indiv_num <- sum(remain)
  if (indiv_num > 0) {
    Q[which(remain == 1), (p1 - indiv_num + 1):p1] <- diag(indiv_num)
  }

  X1 <- X %*% Q
  beta1 <- MASS::ginv(t(X1) %*% X1) %*% t(X1) %*% Y
  target <- t(X) %*% (Y - X1 %*% beta1)
  vals <- sqrt((t(Q) %*% target^2)[seq_len(nrow(valid_set))]) / valid_set$weight
  max(vals) / n
}

acc_prox_simple_linear <- function(
  Y,
  X,
  tree_result,
  lambda,
  warm_start = FALSE,
  init_beta = NULL,
  intercept = FALSE,
  ridge_param = 0,
  thresh = 1e-6,
  max_iter = 10000,
  algos_cpp = NULL,
  stop_rule = c("objective", "coef")
) {
  stop_rule <- match.arg(stop_rule)
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X

  if (lambda < 0) stop("lambda must be non-negative.", call. = FALSE)
  if (ridge_param < 0) stop("ridge_param must be non-negative.", call. = FALSE)

  if (intercept) {
    Y.mean <- mean(Y)
    X.mean <- colMeans(X)
    Y <- Y - Y.mean
    X <- scale(X, scale = FALSE)
  }

  p <- ncol(X)
  if (warm_start && !is.null(init_beta)) {
    if (length(init_beta) != p) stop("init_beta must have length ncol(X).", call. = FALSE)
    init_beta <- as.numeric(init_beta)
  } else {
    init_beta <- rep(1, p)
  }

  fit <- acc_prox_simple_linear_cpp(
    Y = Y,
    X = X,
    tree_result = tree_result,
    lambda = lambda,
    init_beta = init_beta,
    ridge_param = ridge_param,
    thresh = thresh,
    max_iter = as.integer(max_iter),
    stop_rule = as.integer(match(stop_rule, c("objective", "coef")) - 1L)
  )

  if (!isTRUE(fit$converged)) {
    warning("acc_prox_simple_linear reached max_iter before satisfying the stopping rule.", call. = FALSE)
  }

  beta1 <- as.vector(fit$beta)
  if (!intercept) {
    beta1
  } else {
    list(beta0 = as.numeric(Y.mean - sum(X.mean * beta1)), beta1 = beta1)
  }
}

grid.simple_linear <- function(
  Y,
  X,
  tree_df,
  true_beta = NULL,
  ridge.param = 0,
  seqc = NULL,
  thresh = 1e-5,
  intercept = FALSE,
  max_iter = 10000,
  algos_cpp = NULL,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy")
) {
  stop_rule <- match.arg(stop_rule)
  coarest_rule <- match.arg(coarest_rule)
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X
  tree_df <- as.data.frame(tree_df)
  .validate_tree_df(tree_df, require_weight = TRUE)

  n <- length(Y)
  p <- ncol(X)

  if (is.null(true_beta)) true_beta <- rep(0, p)
  if (length(true_beta) != p) stop("true_beta must have length ncol(X).", call. = FALSE)

  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)
  if (coarest_rule == "legacy") {
    coarest_set <- .find_coarest_legacy(tree_df, tree_result$df_result)
    penalty.max <- .find_max_param_linear_legacy(Y - mean(Y), scale(X, scale = FALSE), coarest_set)
  } else {
    coarest_set <- find_coarest(tree_df, tree_result$df_result)
    penalty.max <- find_max_param.linear(Y - mean(Y), scale(X, scale = FALSE), coarest_set)
  }
  if (!is.finite(penalty.max) || penalty.max <= 0) penalty.max <- 1

  if (is.null(seqc)) {
    seqc <- exp(seq(-4 - log(n), log(penalty.max / 5), length = 50))
  }
  seqc <- as.numeric(seqc)
  if (any(!is.finite(seqc)) || any(seqc < 0) || !length(seqc)) {
    stop("seqc must contain finite non-negative values.", call. = FALSE)
  }

  beta <- matrix(0, nrow = p, ncol = length(seqc))
  beta0 <- rep(0, length(seqc))

  fit <- acc_prox_simple_linear(
    Y,
    X,
    tree_result$list_result,
    lambda = seqc[1],
    ridge_param = ridge.param,
    thresh = thresh,
    intercept = intercept,
    max_iter = max_iter,
    algos_cpp = algos_cpp,
    stop_rule = stop_rule
  )

  if (intercept) {
    beta0[1] <- fit$beta0
    beta[, 1] <- fit$beta1
  } else {
    beta[, 1] <- fit
  }

  if (length(seqc) >= 2L) {
    for (i in 2:length(seqc)) {
      fit <- acc_prox_simple_linear(
        Y,
        X,
        tree_result$list_result,
        lambda = seqc[i],
        warm_start = TRUE,
        init_beta = as.vector(beta[, i - 1]),
        ridge_param = ridge.param,
        thresh = thresh,
        intercept = intercept,
        max_iter = max_iter,
        algos_cpp = algos_cpp,
        stop_rule = stop_rule
      )

      if (intercept) {
        beta0[i] <- fit$beta0
        beta[, i] <- fit$beta1
      } else {
        beta[, i] <- fit
      }
    }
  }

  loss <- apply(matrix(rep(true_beta, length(seqc)), nrow = p, byrow = FALSE) - beta, 2, function(x) sum(x^2)) / p
  output <- list(loss = loss, beta = beta, lambda = seqc)
  if (intercept) output$beta0 <- beta0
  output
}

cv.simple_linear <- function(
  Y,
  X,
  tree_df,
  folds = 5,
  seqc = NULL,
  thresh = 1e-5,
  ridge.param = 0,
  intercept = FALSE,
  Mmratio = 10000,
  max_iter = 10000,
  verbose = FALSE,
  algos_cpp = NULL,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy")
) {
  stop_rule <- match.arg(stop_rule)
  coarest_rule <- match.arg(coarest_rule)
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X
  tree_df <- as.data.frame(tree_df)
  .validate_tree_df(tree_df, require_weight = TRUE)

  n <- length(Y)
  p <- ncol(X)
  stopifnot(n >= 2 * folds)

  random_sequence <- sample(seq_len(n))
  index <- cut(random_sequence, breaks = folds, labels = FALSE)

  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)
  if (coarest_rule == "legacy") {
    coarest_set <- .find_coarest_legacy(tree_df, tree_result$df_result)
    penalty.max <- .find_max_param_linear_legacy(Y, X, coarest_set)
  } else {
    coarest_set <- find_coarest(tree_df, tree_result$df_result)
    penalty.max <- find_max_param.linear(Y, X, coarest_set)
  }
  if (!is.finite(penalty.max) || penalty.max <= 0) penalty.max <- 1

  if (is.null(seqc)) {
    seqc <- exp(seq(log(penalty.max / Mmratio), log(penalty.max), length = 50))
  }
  seqc <- sort(as.numeric(seqc), decreasing = TRUE)
  if (any(!is.finite(seqc)) || any(seqc < 0) || !length(seqc)) {
    stop("seqc must contain finite non-negative values.", call. = FALSE)
  }

  vals.mat <- matrix(0, nrow = folds, ncol = length(seqc))
  colnames(vals.mat) <- seqc

  for (i in seq_len(folds)) {
    if (isTRUE(verbose)) message("fold ", i, "/", folds)

    train_idx <- which(index != i)
    test_idx <- which(index == i)

    X_train <- X[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]
    Y_train <- Y[train_idx]
    Y_test <- Y[test_idx]

    res <- grid.simple_linear(
      Y = Y_train,
      X = X_train,
      tree_df = tree_df,
      true_beta = rep(0, p),
      seqc = seqc,
      thresh = thresh,
      ridge.param = ridge.param,
      intercept = intercept,
      max_iter = max_iter,
      algos_cpp = algos_cpp,
      stop_rule = stop_rule,
      coarest_rule = coarest_rule
    )

    pred <- crossprod(t(X_test), res$beta)
    if (intercept && !is.null(res$beta0)) {
      pred <- sweep(pred, 2, res$beta0, "+")
    }
    vals.mat[i, ] <- colMeans((Y_test - pred)^2)
  }

  vals.vec <- colMeans(vals.mat)
  selected.param <- seqc[which.min(vals.vec)]

  final_fit <- acc_prox_simple_linear(
    Y,
    X,
    tree_result$list_result,
    lambda = selected.param,
    intercept = intercept,
    ridge_param = ridge.param,
    thresh = thresh,
    max_iter = max_iter,
    algos_cpp = algos_cpp,
    stop_rule = stop_rule
  )

  if (intercept) {
    list(
      selected.param = selected.param,
      beta0 = final_fit$beta0,
      beta = final_fit$beta1,
      valid.error = vals.vec,
      lambda = seqc
    )
  } else {
    list(selected.param = selected.param, beta = final_fit, valid.error = vals.vec, lambda = seqc)
  }
}

real_data_one_round.linear <- function(
  y,
  X,
  tree_df,
  split.seed = 123,
  nfolds = 5,
  nval = 31,
  ridge_param = 0,
  Mmratio = 10000,
  max_iter = 10000,
  algos_cpp = NULL,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy")
) {
  stop_rule <- match.arg(stop_rule)
  coarest_rule <- match.arg(coarest_rule)
  checked <- .check_yx(y, X)
  y <- checked$Y
  X <- checked$X
  tree_df <- as.data.frame(tree_df)
  .validate_tree_df(tree_df, require_weight = TRUE)
  .require_namespace("glmnet")

  n <- nrow(X)
  if (nval <= 0 || nval >= n) {
    stop("nval must be positive and smaller than nrow(X).", call. = FALSE)
  }
  if ((n - nval) < nfolds) {
    stop("Training sample size must be at least nfolds.", call. = FALSE)
  }

  set.seed(split.seed)
  test.index <- sample(seq_len(n), nval)
  train.index <- setdiff(seq_len(n), test.index)

  X.train <- X[train.index, , drop = FALSE]
  X.test <- X[test.index, , drop = FALSE]
  y.train <- y[train.index]
  y.test <- y[test.index]

  set.seed(2 * split.seed)
  lasso.mod <- glmnet::cv.glmnet(
    X.train,
    y.train,
    family = "gaussian",
    alpha = 1,
    intercept = FALSE,
    nfolds = nfolds
  )
  lasso.fit <- glmnet::glmnet(
    X.train,
    y.train,
    family = "gaussian",
    alpha = 1,
    intercept = FALSE,
    lambda = lasso.mod$lambda.min
  )
  lasso.beta <- as.numeric(lasso.fit$beta)
  lasso.loss <- mean((y.test - as.numeric(X.test %*% lasso.beta))^2)
  lasso.param <- lasso.mod$lambda.min

  set.seed(2 * split.seed)
  ridge.mod <- glmnet::cv.glmnet(
    X.train,
    y.train,
    family = "gaussian",
    alpha = 0,
    intercept = FALSE,
    nfolds = nfolds
  )
  ridge.fit <- glmnet::glmnet(
    X.train,
    y.train,
    family = "gaussian",
    alpha = 0,
    intercept = FALSE,
    lambda = ridge.mod$lambda.min
  )
  ridge.beta <- as.numeric(ridge.fit$beta)
  ridge.loss <- mean((y.test - as.numeric(X.test %*% ridge.beta))^2)
  ridge.param <- ridge.mod$lambda.min

  set.seed(2 * split.seed)
  our.mod <- cv.simple_linear(
    y.train,
    X.train,
    tree_df,
    folds = nfolds,
    ridge.param = ridge_param,
    Mmratio = Mmratio,
    max_iter = max_iter,
    algos_cpp = algos_cpp,
    stop_rule = stop_rule,
    coarest_rule = coarest_rule
  )
  our.beta <- our.mod$beta
  our.loss <- mean((y.test - as.numeric(X.test %*% our.beta))^2)
  our.param <- our.mod$selected.param

  list(
    split.seed = split.seed,
    train.index = train.index,
    test.index = test.index,
    lasso.beta = lasso.beta,
    lasso.loss = lasso.loss,
    lasso.param = lasso.param,
    ridge.beta = ridge.beta,
    ridge.loss = ridge.loss,
    ridge.param = ridge.param,
    our.beta = our.beta,
    our.loss = our.loss,
    our.param = our.param,
    our.valid.error = our.mod$valid.error,
    our.lambda = our.mod$lambda
  )
}
