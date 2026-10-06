small_tree <- function(root_weight = 1) {
  data.frame(
    node = 1:5,
    parent = c(4, 4, 5, 5, NA),
    weight = c(0, 0, 0, 1, root_weight)
  )
}

deep_tree <- function() {
  data.frame(
    node = 1:7,
    parent = c(5, 5, 6, 6, 7, 7, NA),
    weight = c(0, 0, 0, 0, 1.0, 0.8, 0.5)
  )
}

reference_one_layer <- function(eta, groups, weights, lambda) {
  beta <- eta

  for (i in seq_along(weights)) {
    ind <- groups[[i]]
    if (!length(ind)) next

    mu <- eta[ind]
    mean_mu <- mean(mu)
    temp <- sqrt(sum((mu - mean_mu)^2))
    d <- 100
    if (temp > 0) d <- weights[i] * lambda / temp

    if (d >= 1) {
      beta[ind] <- mean_mu
    } else {
      beta[ind] <- (1 - d) * mu + d * mean_mu
    }
  }

  as.numeric(beta)
}

reference_prox_tree <- function(eta, lambda, tree_result) {
  new_eta <- as.numeric(eta)

  for (i in seq_along(tree_result$groups)) {
    new_eta <- reference_one_layer(
      eta = new_eta,
      groups = tree_result$groups[[i]],
      weights = tree_result$weights[[i]],
      lambda = lambda
    )
  }

  new_eta
}

reference_acc_core <- function(
  Y,
  X,
  tree_result,
  lambda,
  init_beta,
  ridge_param = 0,
  thresh = 1e-6,
  max_iter = 10000
) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  n <- nrow(X)

  beta0 <- beta1 <- as.numeric(init_beta)
  alpha0 <- 1
  alpha1 <- 0.5

  matXX <- crossprod(X)
  eig_vals <- eigen(matXX, symmetric = TRUE, only.values = TRUE)$values
  L0 <- max(eig_vals, na.rm = TRUE) / n
  if (!is.finite(L0) || L0 <= 0) L0 <- .Machine$double.eps
  tao <- n / L0

  matXY <- as.numeric(crossprod(X, Y))
  iter <- 0L
  consecutive_below_thresh <- 0L

  while (consecutive_below_thresh < 5L && iter < max_iter) {
    iter <- iter + 1L
    accept <- FALSE
    Gam <- beta0 + (alpha0 - 1) / alpha1 * (beta1 - beta0)

    matXXGam <- as.numeric(matXX %*% Gam)
    matXGam <- as.numeric(X %*% Gam)
    base1 <- -matXY + matXXGam + ridge_param * Gam

    while (!accept) {
      eta <- Gam - tao * base1 / n
      eta.new <- reference_prox_tree(eta, lambda = tao * lambda, tree_result = tree_result)

      lhs <- sum((Y - as.numeric(X %*% eta.new))^2) / (2 * n)
      rhs <- sum((Y - matXGam)^2) / (2 * n) +
        sum((matXXGam - matXY) * (eta.new - Gam)) +
        sum((eta.new - Gam)^2) / (2 * tao)

      if (lhs <= rhs || tao <= 1 / L0) {
        accept <- TRUE
      } else {
        tao <- max(tao / 2, 1 / L0)
      }
    }

    beta0 <- beta1
    beta1 <- eta.new

    old.obj <- sum((Y - as.numeric(X %*% beta0))^2) / (2 * n)
    new.obj <- sum((Y - as.numeric(X %*% beta1))^2) / (2 * n)
    dis <- abs(new.obj - old.obj) / max(abs(old.obj), .Machine$double.eps)

    if (dis < thresh) {
      consecutive_below_thresh <- consecutive_below_thresh + 1L
    } else {
      consecutive_below_thresh <- 0L
    }

    alpha0 <- alpha1
    alpha1 <- (1 + sqrt(1 + 4 * alpha0^2)) / 2
  }

  list(
    beta = as.numeric(beta1),
    iter = iter,
    converged = consecutive_below_thresh >= 5L
  )
}

reference_max_param <- function(Y, X, coarest_set) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  valid_set <- coarest_set[!is.na(coarest_set$weight) & coarest_set$weight != 0, , drop = FALSE]

  if (!nrow(valid_set)) {
    r <- qr.resid(qr(X), Y)
    return(max(abs(as.numeric(crossprod(X, r)))) / n)
  }

  covered <- sort(unique(unlist(valid_set$leaves, use.names = FALSE)))
  p1 <- nrow(valid_set) + (p - length(covered))
  Q <- matrix(0, nrow = p, ncol = p1)

  for (i in seq_len(nrow(valid_set))) {
    Q[valid_set$leaves[[i]], i] <- 1
  }

  singletons <- setdiff(seq_len(p), covered)
  if (length(singletons)) {
    cols <- seq_len(length(singletons)) + nrow(valid_set)
    Q[cbind(singletons, cols)] <- 1
  }

  X1 <- X %*% Q
  r <- qr.resid(qr(X1), Y)
  g2 <- as.numeric(crossprod(X, r))^2
  vals <- sqrt(as.numeric(t(Q) %*% g2)[seq_len(nrow(valid_set))]) / valid_set$weight

  max(vals) / n
}
