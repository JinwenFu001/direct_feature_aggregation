.log1pexp <- function(x) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  pos <- x > 0
  out[pos] <- x[pos] + log1p(exp(-x[pos]))
  out[!pos] <- log1p(exp(x[!pos]))
  out
}

.sigmoid <- function(x) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  exp_x <- exp(x[!pos])
  out[!pos] <- exp_x / (1 + exp_x)
  out
}

.check_binary_y <- function(Y, allow_probabilities = FALSE) {
  if (allow_probabilities) {
    if (any(Y < 0 | Y > 1)) {
      stop("Y must contain values in [0, 1].", call. = FALSE)
    }
  } else if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary and encoded as 0/1.", call. = FALSE)
  }
  invisible(TRUE)
}

.logistic_objective <- function(Y, X, beta, ridge_param = 0) {
  eta <- as.numeric(X %*% beta)
  mean(.log1pexp(eta) - Y * eta) + ridge_param * sum(beta^2) / (2 * nrow(X))
}

negtv_lglkh <- function(Y, X, beta) {
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X
  .check_binary_y(Y, allow_probabilities = TRUE)

  beta <- as.matrix(beta)
  if (nrow(beta) != ncol(X)) {
    stop("beta must have nrow equal to ncol(X).", call. = FALSE)
  }

  eta <- X %*% beta
  log_terms <- matrix(.log1pexp(as.vector(eta)), nrow = nrow(eta))
  unname(colMeans(log_terms - Y * eta))
}

find_max_param.logistic <- function(Y, X, coarest_set) {
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X
  .check_binary_y(Y)

  n <- nrow(X)
  p <- ncol(X)
  all_leaves <- seq_len(p)

  .check_required_columns(coarest_set, c("leaves", "weight"), "coarest_set")
  valid_set <- coarest_set[!is.na(coarest_set$weight) & coarest_set$weight != 0, , drop = FALSE]

  if (!nrow(valid_set)) {
    g <- as.numeric(crossprod(X, Y - 0.5))
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
  fit <- suppressWarnings(stats::glm.fit(
    x = X1,
    y = Y,
    family = stats::binomial(link = "logit"),
    intercept = FALSE
  ))

  beta1 <- as.numeric(fit$coefficients)
  beta1[!is.finite(beta1)] <- 0
  prob <- .sigmoid(as.numeric(X1 %*% beta1))
  target <- as.numeric(crossprod(X, Y - prob))
  vals <- sqrt(as.numeric(t(Q) %*% (target^2))[seq_len(nrow(valid_set))]) / valid_set$weight

  max(vals) / n
}

acc_prox_simple_logistic <- function(
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
  stop_rule = c("objective", "coef")
) {
  stop_rule <- match.arg(stop_rule)
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X
  .check_binary_y(Y)

  if (lambda < 0) stop("lambda must be non-negative.", call. = FALSE)
  if (ridge_param < 0) stop("ridge_param must be non-negative.", call. = FALSE)

  if (intercept) {
    Y.mean <- mean(Y)
    X.mean <- colMeans(X)
    Y <- Y - Y.mean
    X <- scale(X, scale = FALSE)
  }

  p <- ncol(X)
  n <- nrow(X)
  if (warm_start && !is.null(init_beta)) {
    if (length(init_beta) != p) stop("init_beta must have length ncol(X).", call. = FALSE)
    beta0 <- beta1 <- as.numeric(init_beta)
  } else {
    beta0 <- beta1 <- rep(1, p)
  }

  matXX <- crossprod(X)
  eig_vals <- eigen(matXX, symmetric = TRUE, only.values = TRUE)$values
  L0 <- (max(eig_vals, na.rm = TRUE) + ridge_param) / n
  if (!is.finite(L0) || L0 <= 0) L0 <- .Machine$double.eps
  tao <- n / L0

  alpha0 <- 1
  alpha1 <- 0.5
  iter <- 0L
  consecutive_below_thresh <- 0L

  while (consecutive_below_thresh < 5L && iter < max_iter) {
    iter <- iter + 1L
    accept <- FALSE
    Gam <- beta0 + (alpha0 - 1) / alpha1 * (beta1 - beta0)
    eta_gam <- as.numeric(X %*% Gam)
    prob_gam <- .sigmoid(eta_gam)
    grad <- as.numeric(crossprod(X, prob_gam - Y) + ridge_param * Gam) / n
    obj_gam <- mean(.log1pexp(eta_gam) - Y * eta_gam) + ridge_param * sum(Gam^2) / (2 * n)

    while (!accept) {
      eta <- Gam - tao * grad
      eta.new <- prox_tree_cpp(eta, lambda = tao * lambda, tree_list = tree_result)
      diff <- eta.new - Gam

      lhs <- .logistic_objective(Y, X, eta.new, ridge_param = ridge_param)
      rhs <- obj_gam + sum(grad * diff) + sum(diff^2) / (2 * tao)

      if (lhs <= rhs || tao <= 1 / L0) {
        accept <- TRUE
      } else {
        tao <- max(tao / 2, 1 / L0)
      }
    }

    beta0 <- beta1
    beta1 <- eta.new

    if (stop_rule == "objective") {
      old.obj <- .logistic_objective(Y, X, beta0, ridge_param = ridge_param)
      new.obj <- .logistic_objective(Y, X, beta1, ridge_param = ridge_param)
      dis <- abs(new.obj - old.obj) / max(abs(old.obj), .Machine$double.eps)
    } else {
      dis <- sqrt(sum((beta1 - beta0)^2)) / max(sqrt(sum(beta0^2)), .Machine$double.eps)
    }

    if (dis < thresh) {
      consecutive_below_thresh <- consecutive_below_thresh + 1L
    } else {
      consecutive_below_thresh <- 0L
    }

    alpha0 <- alpha1
    alpha1 <- (1 + sqrt(1 + 4 * alpha0^2)) / 2
  }

  if (iter >= max_iter && consecutive_below_thresh < 5L) {
    warning("acc_prox_simple_logistic reached max_iter before satisfying the stopping rule.", call. = FALSE)
  }

  if (!intercept) {
    as.vector(beta1)
  } else {
    list(beta0 = as.numeric(Y.mean - sum(X.mean * beta1)), beta1 = as.vector(beta1))
  }
}

grid.logistic <- function(
  Y,
  X,
  tree_df,
  true_beta = NULL,
  ridge.param = 0,
  seqc = NULL,
  thresh = 1e-5,
  intercept = FALSE,
  max_iter = 10000,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy")
) {
  stop_rule <- match.arg(stop_rule)
  coarest_rule <- match.arg(coarest_rule)
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X
  .check_binary_y(Y)
  tree_df <- as.data.frame(tree_df)
  .validate_tree_df(tree_df, require_weight = TRUE)

  n <- length(Y)
  p <- ncol(X)

  if (is.null(true_beta)) true_beta <- rep(0, p)
  if (length(true_beta) != p) stop("true_beta must have length ncol(X).", call. = FALSE)

  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)
  if (coarest_rule == "legacy") {
    coarest_set <- .find_coarest_legacy(tree_df, tree_result$df_result)
  } else {
    coarest_set <- find_coarest(tree_df, tree_result$df_result)
  }
  penalty.max <- find_max_param.logistic(Y, X, coarest_set)
  if (!is.finite(penalty.max) || penalty.max <= 0) penalty.max <- 1

  if (is.null(seqc)) {
    seqc <- exp(seq(-4 - log(n), log(penalty.max), length = 50))
  }
  seqc <- as.numeric(seqc)
  if (any(!is.finite(seqc)) || any(seqc < 0) || !length(seqc)) {
    stop("seqc must contain finite non-negative values.", call. = FALSE)
  }

  beta <- matrix(0, nrow = p, ncol = length(seqc))
  beta0 <- rep(0, length(seqc))

  fit <- acc_prox_simple_logistic(
    Y,
    X,
    tree_result$list_result,
    lambda = seqc[1],
    ridge_param = ridge.param,
    thresh = thresh,
    intercept = intercept,
    max_iter = max_iter,
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
      fit <- acc_prox_simple_logistic(
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

cv.logistic <- function(
  Y,
  X,
  tree_df,
  folds = 5,
  seqc = NULL,
  neg.ind = NULL,
  thresh = 1e-5,
  ridge.param = 0,
  intercept = FALSE,
  Mmratio = 10000,
  max_iter = 10000,
  verbose = FALSE,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy")
) {
  stop_rule <- match.arg(stop_rule)
  coarest_rule <- match.arg(coarest_rule)
  checked <- .check_yx(Y, X)
  Y <- checked$Y
  X <- checked$X
  .check_binary_y(Y)
  tree_df <- as.data.frame(tree_df)
  .validate_tree_df(tree_df, require_weight = TRUE)

  n <- length(Y)
  p <- ncol(X)
  stopifnot(n >= 2 * folds)

  if (is.null(neg.ind)) {
    random_sequence <- sample(seq_len(n))
    index <- cut(random_sequence, breaks = folds, labels = FALSE)
  } else {
    pos.ind.rand <- sample(setdiff(seq_len(n), neg.ind))
    neg.ind.rand <- sample(neg.ind)
    pos.partition <- cut(pos.ind.rand, breaks = folds, labels = FALSE)
    neg.partition <- cut(neg.ind.rand, breaks = folds, labels = FALSE)
    index <- numeric(n)
    index[setdiff(seq_len(n), neg.ind)] <- pos.partition
    index[neg.ind] <- neg.partition
  }

  tree_result <- gather_leaf_nodes_per_non_leaf(tree_df)
  if (coarest_rule == "legacy") {
    coarest_set <- .find_coarest_legacy(tree_df, tree_result$df_result)
  } else {
    coarest_set <- find_coarest(tree_df, tree_result$df_result)
  }
  penalty.max <- find_max_param.logistic(Y, X, coarest_set)
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

    res <- grid.logistic(
      Y = Y[train_idx],
      X = X[train_idx, , drop = FALSE],
      tree_df = tree_df,
      true_beta = rep(0, p),
      seqc = seqc,
      thresh = thresh,
      ridge.param = ridge.param,
      intercept = intercept,
      max_iter = max_iter,
      stop_rule = stop_rule,
      coarest_rule = coarest_rule
    )

    pred_beta <- res$beta
    if (intercept && !is.null(res$beta0)) {
      eta <- X[test_idx, , drop = FALSE] %*% pred_beta
      eta <- sweep(eta, 2, res$beta0, "+")
      log_terms <- matrix(.log1pexp(as.vector(eta)), nrow = nrow(eta))
      vals.mat[i, ] <- colMeans(log_terms - Y[test_idx] * eta)
    } else {
      vals.mat[i, ] <- negtv_lglkh(Y[test_idx], X[test_idx, , drop = FALSE], beta = pred_beta)
    }
  }

  vals.vec <- colMeans(vals.mat)
  selected.param <- seqc[which.min(vals.vec)]
  final_fit <- acc_prox_simple_logistic(
    Y,
    X,
    tree_result$list_result,
    lambda = selected.param,
    intercept = intercept,
    ridge_param = ridge.param,
    thresh = thresh,
    max_iter = max_iter,
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
