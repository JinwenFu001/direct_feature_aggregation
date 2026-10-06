# Helper functions for the uu=1 tree-change task.
# These are simulation, tree-construction, and baseline-evaluation helpers.
# The direct-penalty method itself comes from the treeFA package.

ensure_task_packages <- function() {
  required <- c("treeFA", "rare", "mclust", "glmnet")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      ". Install them before running this task.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

simulate_tree <- function(k, p, tao = 0.05) {
  means <- 1 / (1:k)
  group.index <- c(
    cut(1:ceiling(p * 3 / 4), breaks = k / 2, labels = FALSE),
    cut((ceiling(p * 3 / 4) + 1):p, breaks = k / 2, labels = FALSE) + k / 2
  )
  group.size <- table(group.index)
  latent <- c()

  for (i in 1:k) {
    near.index <- ifelse(i < k, i + 1, i - 1)
    latent <- c(latent, stats::rnorm(group.size[i], means[i], tao * abs(means[i] - means[near.index])))
  }

  tree <- stats::hclust(stats::dist(latent))
  list(tree = tree, group = group.index)
}

assign_weight <- function(tree_df, p, weight.order = -1) {
  tree_df$weight[(p + 1):nrow(tree_df)] <- 1
  result <- treeFA::gather_leaf_nodes_per_non_leaf(tree_df)$df_result
  depth <- numeric(nrow(result))

  for (i in seq_len(nrow(result))) {
    depth[i] <- length(result$leaves[[i]])^weight.order
    tree_df$weight[tree_df$node == result$node[i]] <- depth[i]
  }

  tree_df$weight <- tree_df$weight / mean(depth)
  tree_df
}

hclust_to_df <- function(hc, weight.order = -1) {
  tree_df <- data.frame(
    node = integer(),
    parent = integer(),
    name = character(),
    weight = numeric(),
    stringsAsFactors = FALSE
  )

  traverse <- function(node, parent, name_prefix, depth) {
    if (node < 0) {
      node_id <- -node
      name <- paste("leaf", node_id, sep = "")
      weight <- 0
    } else {
      node_id <- nrow(hc$merge) + node + 1
      name <- paste("latent", node_id, sep = "")
      weight <- depth * 10
    }

    tree_df <<- rbind(
      tree_df,
      data.frame(
        node = node_id,
        parent = parent,
        name = name,
        weight = weight,
        stringsAsFactors = FALSE
      )
    )

    if (node > 0) {
      left_child <- hc$merge[node, 1]
      right_child <- hc$merge[node, 2]
      traverse(left_child, node_id, name, depth + 1)
      traverse(right_child, node_id, name, depth + 1)
    }
  }

  traverse(nrow(hc$merge), NA, "root", 0)
  tree_df <- tree_df[order(tree_df$node), , drop = FALSE]
  rownames(tree_df) <- NULL
  assign_weight(tree_df, length(hc$order), weight.order = weight.order)
}

simulate_data <- function(n, group.index, s = 0, beta.pre = NULL, ratio = 5) {
  k <- length(unique(group.index))
  p <- length(group.index)
  A <- matrix(0, nrow = p, ncol = k)

  for (i in 1:k) {
    A[, i] <- group.index == i
  }

  if (is.null(beta.pre)) {
    beta0 <- stats::runif(k, 1.5, 2.5)
    nonzero.ind <- sample(1:k, ceiling(k * (1 - s)))
    nonzero <- rep(0, k)
    nonzero[nonzero.ind] <- rep(c(-1, 1), k %/% 2)[seq_along(nonzero.ind)]
    beta0 <- beta0 * nonzero
  } else {
    beta0 <- beta.pre
  }

  beta <- A %*% beta0
  X <- matrix(stats::rpois(n * p, 0.02 * 5 / (p / k)), nrow = n, ncol = p)
  scale.factor <- sqrt(eigen(crossprod(X))$values[1])
  X <- X / scale.factor
  beta <- beta * scale.factor
  means <- X %*% beta
  sigma <- sqrt(sum(means^2) / ratio / n)
  Y <- means + stats::rnorm(n, 0, sigma)

  list(X = X, Y = as.numeric(Y), true.beta = as.numeric(beta), A = A)
}

make_new_X <- function(X, tree_group) {
  index <- tree_group
  levels <- unique(index)
  Q <- matrix(0, nrow = length(tree_group), ncol = length(levels))

  for (i in seq_along(levels)) {
    Q[index == levels[i], i] <- 1
  }

  list(X1 = X %*% Q, Q = Q)
}

my_ginv <- function(X, tol = 1e-8) {
  s <- svd(X)
  d <- s$d
  d[d < tol * max(d)] <- 0
  d_inv <- 1 / d
  d_inv[!is.finite(d_inv)] <- 0
  s$v %*% diag(d_inv, length(d_inv)) %*% t(s$u)
}

find_latent_nodes_with_exact_leaf_groups <- function(tree_df, leaf_group_vector) {
  root_node <- tree_df$node[is.na(tree_df$parent)]
  all_leaves <- treeFA::gather_leaves(root_node, tree_df)

  leaf_groups_df <- data.frame(
    leaf = all_leaves,
    group = leaf_group_vector
  )
  group_sets <- lapply(split(leaf_groups_df$leaf, leaf_groups_df$group), unique)

  gather_info <- treeFA::gather_leaf_nodes_per_non_leaf(tree_df)
  non_leaf_nodes <- gather_info$df_result

  same_set <- function(a, b) {
    if (length(a) != length(b)) return(FALSE)
    all(sort(a) == sort(b))
  }

  matching_latent_nodes <- vector("list", length(group_sets))
  names(matching_latent_nodes) <- names(group_sets)

  for (g in names(group_sets)) {
    group_leaves <- group_sets[[g]]
    matched_nodes <- non_leaf_nodes$node[
      vapply(non_leaf_nodes$leaves, function(node_leaf_vec) {
        same_set(node_leaf_vec, group_leaves)
      }, logical(1))
    ]
    matching_latent_nodes[[g]] <- matched_nodes
  }

  matching_latent_nodes
}

find_all_parents <- function(tree_df, vec) {
  vec <- vec[!is.na(vec)]
  if (!length(vec)) return(integer())

  parent_vec <- unique(tree_df$parent[match(vec, tree_df$node)])
  old_len <- -1L

  while (length(parent_vec) != old_len) {
    old_len <- length(parent_vec)
    parent_vec <- unique(c(parent_vec, tree_df$parent[match(parent_vec, tree_df$node)]))
    parent_vec <- parent_vec[!is.na(parent_vec)]
  }

  parent_vec
}

simulate_error_given_tree <- function(
  n,
  n0,
  tree_df,
  tree_group,
  tree_hclust,
  n1 = n,
  s = 0,
  reps = 50,
  ridge.param = 0,
  thresh = 1e-5,
  INratio = 5,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy")
) {
  stop_rule <- match.arg(stop_rule)
  coarest_rule <- match.arg(coarest_rule)
  ensure_task_packages()
  cat("Direct penalty settings: stop_rule =", stop_rule, ", coarest_rule =", coarest_rule, "\n")

  rare.minloss <- numeric(reps)
  our.minloss <- numeric(reps)
  rare.rand <- numeric(reps)
  our.rand <- numeric(reps)
  rare.minloss.ideal <- numeric(reps)
  our.minloss.ideal <- numeric(reps)
  rare.rand.ideal <- numeric(reps)
  our.rand.ideal <- numeric(reps)
  oracle.ls.loss <- numeric(reps)
  ls.loss <- numeric(reps)
  ridge.loss <- numeric(reps)
  oracle.ridge.loss <- numeric(reps)

  for (i in seq_len(reps)) {
    cat(i, "/", reps, "\n")
    data <- simulate_data(n + n0 + n1, tree_group, s = s, ratio = INratio)
    p <- ncol(data$X)

    train_index <- sample(seq_len(n + n0 + n1), n)
    valid_index <- sample(setdiff(seq_len(n + n0 + n1), train_index), n1)
    test_index <- setdiff(seq_len(n + n0 + n1), c(train_index, valid_index))

    rare.result <- rare::rarefit(
      data$Y[train_index],
      data$X[train_index, , drop = FALSE],
      hc = tree_hclust,
      alpha = 1,
      intercept = FALSE
    )
    loss1 <- apply(data$X[valid_index, , drop = FALSE] %*% rare.result$beta[[1]], 2, function(x) {
      sum((x - data$Y[valid_index])^2)
    })
    loss3 <- apply(data$X[valid_index, , drop = FALSE] %*% rare.result$beta[[1]], 2, function(x) {
      sum((x - data$X[valid_index, , drop = FALSE] %*% data$true.beta)^2)
    })

    rare.beta <- rare.result$beta[[1]][, which.min(loss1)]
    rare.beta.ideal <- rare.result$beta[[1]][, which.min(loss3)]
    rare.minloss[i] <- sum((data$Y[test_index] - data$X[test_index, , drop = FALSE] %*% rare.beta)^2) / n0
    rare.minloss.ideal[i] <- sum((data$Y[test_index] - data$X[test_index, , drop = FALSE] %*% rare.beta.ideal)^2) / n0

    our.result <- treeFA::grid.simple_linear(
      Y = data$Y[train_index],
      X = data$X[train_index, , drop = FALSE],
      tree_df = tree_df,
      true_beta = data$true.beta,
      ridge.param = ridge.param,
      thresh = thresh,
      stop_rule = stop_rule,
      coarest_rule = coarest_rule
    )
    loss2 <- apply(data$X[valid_index, , drop = FALSE] %*% our.result$beta, 2, function(x) {
      sum((x - data$Y[valid_index])^2)
    })
    loss4 <- apply(data$X[valid_index, , drop = FALSE] %*% our.result$beta, 2, function(x) {
      sum((x - data$X[valid_index, , drop = FALSE] %*% data$true.beta)^2)
    })

    our.beta <- our.result$beta[, which.min(loss2)]
    our.beta.ideal <- our.result$beta[, which.min(loss4)]
    our.minloss[i] <- sum((data$Y[test_index] - data$X[test_index, , drop = FALSE] %*% our.beta)^2) / n0
    our.minloss.ideal[i] <- sum((data$Y[test_index] - data$X[test_index, , drop = FALSE] %*% our.beta.ideal)^2) / n0

    rare.rand[i] <- mclust::adjustedRandIndex(as.numeric(as.factor(data$true.beta)), as.numeric(as.factor(rare.beta)))
    rare.rand.ideal[i] <- mclust::adjustedRandIndex(as.numeric(as.factor(data$true.beta)), as.numeric(as.factor(rare.beta.ideal)))
    our.rand[i] <- mclust::adjustedRandIndex(as.numeric(as.factor(data$true.beta)), as.numeric(as.factor(our.beta)))
    our.rand.ideal[i] <- mclust::adjustedRandIndex(as.numeric(as.factor(data$true.beta)), as.numeric(as.factor(our.beta.ideal)))

    ls.X <- data$X[train_index, , drop = FALSE]
    ls.Y <- data$Y[train_index]
    ls.ginv <- my_ginv(crossprod(ls.X))
    ls.coef <- crossprod(ls.ginv, crossprod(ls.X, ls.Y))
    ls.loss[i] <- sum((data$Y[test_index] - data$X[test_index, , drop = FALSE] %*% ls.coef)^2) / n0

    Q <- make_new_X(data$X, tree_group)
    ols.X <- Q$X1[train_index, , drop = FALSE]
    ols.ginv <- t(my_ginv(crossprod(ols.X)))
    ols.coef <- crossprod(ols.ginv, crossprod(ols.X, ls.Y))
    oracle.ls.loss[i] <- sum((data$Y[test_index] - Q$X1[test_index, , drop = FALSE] %*% ols.coef)^2) / n0

    X_train <- data$X[train_index, , drop = FALSE]
    Y_train <- data$Y[train_index]
    X_valid <- data$X[valid_index, , drop = FALSE]
    Y_valid <- data$Y[valid_index]
    X_test <- data$X[test_index, , drop = FALSE]
    Y_test <- data$Y[test_index]

    fit_ridge <- glmnet::glmnet(X_train, Y_train, alpha = 0, intercept = FALSE)
    lam_seq <- fit_ridge$lambda
    preds_valid_ridge <- stats::predict(fit_ridge, newx = X_valid, s = lam_seq)
    val_mse_ridge <- colMeans((Y_valid - preds_valid_ridge)^2)
    best_lam_std <- lam_seq[which.min(val_mse_ridge)]
    beta_ridge_best_full <- as.numeric(stats::coef(fit_ridge, s = best_lam_std))
    beta_ridge_best <- beta_ridge_best_full[-1]
    preds_ridge_test <- X_test %*% beta_ridge_best
    ridge.loss[i] <- mean((Y_test - preds_ridge_test)^2)

    or.X_train <- Q$X1[train_index, , drop = FALSE]
    or.X_valid <- Q$X1[valid_index, , drop = FALSE]
    or.X_test <- Q$X1[test_index, , drop = FALSE]
    fit_oracle_ridge <- glmnet::glmnet(or.X_train, Y_train, alpha = 0, intercept = FALSE)
    lam_seq_or <- fit_oracle_ridge$lambda
    preds_valid_or <- stats::predict(fit_oracle_ridge, newx = or.X_valid, s = lam_seq_or)
    val_mse_oracle_ridge <- colMeans((Y_valid - preds_valid_or)^2)
    best_lam_oracle <- lam_seq_or[which.min(val_mse_oracle_ridge)]
    beta_group_best_full <- as.numeric(stats::coef(fit_oracle_ridge, s = best_lam_oracle))
    beta_group_best <- beta_group_best_full[-1]
    beta_or_feature <- Q$Q %*% beta_group_best
    preds_test_or <- or.X_test %*% beta_group_best
    oracle.ridge.loss[i] <- mean((Y_test - preds_test_or)^2)

    stopifnot(length(beta_or_feature) == p)
  }

  list(
    our.minloss = our.minloss,
    rare.minloss = rare.minloss,
    our.rand = our.rand,
    rare.rand = rare.rand,
    our.minloss.ideal = our.minloss.ideal,
    rare.minloss.ideal = rare.minloss.ideal,
    our.rand.ideal = our.rand.ideal,
    rare.rand.ideal = rare.rand.ideal,
    oracle.ls.loss = oracle.ls.loss,
    ls.loss = ls.loss,
    ridge.loss = ridge.loss,
    oracle.ridge.loss = oracle.ridge.loss
  )
}
