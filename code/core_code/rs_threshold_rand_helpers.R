# Helpers for computing RS-DL2/RS-CL2 Rand index after penalty-group thresholding.
# These functions do not change selected gamma, fitted beta, or prediction error.

rs_threshold_method_key <- function(method) {
  gsub("-", "_", method, fixed = TRUE)
}

rs_threshold_group_ids_from_beta <- function(beta, tol = 1e-4) {
  x <- round(as.numeric(beta) / tol) * tol
  x[abs(x) < tol] <- 0
  as.integer(factor(
    format(x, scientific = TRUE, digits = 12),
    levels = unique(format(x, scientific = TRUE, digits = 12))
  ))
}

rs_threshold_tree_helpers <- function(tree_df, p) {
  tree_df <- as.data.frame(tree_df)
  tree_df$node <- as.integer(tree_df$node)
  tree_df$parent <- as.integer(tree_df$parent)
  root <- tree_df$node[is.na(tree_df$parent)]
  if (length(root) != 1L) stop("tree_df must have exactly one root.", call. = FALSE)

  children <- split(tree_df$node[!is.na(tree_df$parent)], tree_df$parent[!is.na(tree_df$parent)])

  descendants <- function(node) {
    frontier <- as.integer(children[[as.character(node)]])
    out <- integer()
    while (length(frontier)) {
      out <- c(out, frontier)
      next_frontier <- unlist(children[as.character(frontier)], use.names = FALSE)
      frontier <- as.integer(next_frontier)
    }
    out
  }

  leaves_under <- function(node) {
    sort(intersect(c(node, descendants(node)), seq_len(p)))
  }

  depth <- function(node) {
    d <- 0L
    current <- node
    repeat {
      parent <- tree_df$parent[match(current, tree_df$node)]
      if (is.na(parent)) break
      d <- d + 1L
      current <- as.integer(parent)
    }
    d
  }

  list(
    leaves_under = leaves_under,
    depth = Vectorize(depth, USE.NAMES = FALSE)
  )
}

rs_threshold_penalty_nodes <- function(tree_df) {
  tree_df <- as.data.frame(tree_df)
  tree_df$node <- as.integer(tree_df$node)
  sort(unique(tree_df$parent[!is.na(tree_df$parent)]))
}

rs_threshold_shrink_group <- function(A, g_idx) {
  R <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  for (v in seq_len(nrow(g_idx))) {
    idx <- g_idx[v, 1]:g_idx[v, 2]
    gnorm <- sqrt(sum(A[idx, , drop = FALSE]^2))
    if (gnorm < 1) gnorm <- 1
    R[idx, ] <- A[idx, , drop = FALSE] / gnorm
  }
  R
}

rs_threshold_cal2norm <- function(A, g_idx) {
  total <- 0
  for (v in seq_len(nrow(g_idx))) {
    idx <- g_idx[v, 1]:g_idx[v, 2]
    total <- total + sqrt(sum(A[idx, , drop = FALSE]^2))
  }
  total
}

rs_threshold_spg_group_one_gamma <- function(
  y,
  X,
  gamma,
  C,
  CNorm,
  g_idx,
  mu = 1e-2,
  maxiter = 10000L,
  tol = 1e-7
) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  C <- as.matrix(C)
  g_idx <- round(as.matrix(g_idx))

  J <- ncol(X)
  beta <- rep(0, J)
  w <- beta
  theta <- 1
  C_scaled <- C * gamma
  XX <- crossprod(X)
  XY <- crossprod(X, y)
  L <- as.numeric(max(eigen(XX, symmetric = TRUE, only.values = TRUE)$values)) +
    gamma^2 * CNorm / mu
  if (!is.finite(L) || L <= 0) stop("Invalid Lipschitz constant.", call. = FALSE)

  prev_obj <- Inf
  obj <- Inf
  for (iter in seq_len(maxiter)) {
    A <- rs_threshold_shrink_group(C_scaled %*% w / mu, g_idx)
    grad <- as.numeric(XX %*% w - XY + crossprod(C_scaled, A))
    beta_new <- w - grad / L
    theta_new <- (sqrt(theta^4 + 4 * theta^2) - theta^2) / 2
    w <- beta_new + (1 - theta) / theta * theta_new * (beta_new - beta)

    obj <- sum((y - X %*% beta_new)^2) / 2 +
      rs_threshold_cal2norm(C_scaled %*% beta_new, g_idx)
    rel_change <- abs(obj - prev_obj) / max(abs(prev_obj), .Machine$double.eps)
    beta <- beta_new
    theta <- theta_new
    if (iter > 10L && is.finite(rel_change) && rel_change < tol) break
    prev_obj <- obj
  }

  list(coef = beta, obj = obj, iter = iter)
}

rs_threshold_group_norms <- function(C, coef, g_idx) {
  Ccoef <- as.matrix(C) %*% as.numeric(coef)
  vapply(seq_len(nrow(g_idx)), function(i) {
    idx <- g_idx[i, 1]:g_idx[i, 2]
    sqrt(sum(Ccoef[idx]^2))
  }, numeric(1))
}

rs_threshold_apply_top_down <- function(beta, tree_df, norms, threshold, p) {
  helpers <- rs_threshold_tree_helpers(tree_df, p = p)
  nodes <- rs_threshold_penalty_nodes(tree_df)
  if (length(nodes) != length(norms)) {
    stop("Number of penalty nodes does not match number of group norms.", call. = FALSE)
  }

  leaf_sets <- lapply(nodes, helpers$leaves_under)
  sizes <- vapply(leaf_sets, length, integer(1))
  depths <- helpers$depth(nodes)
  order_idx <- order(depths, -sizes, nodes)

  post_beta <- as.numeric(beta)
  assigned <- rep(FALSE, p)
  triggered <- list()

  for (idx in order_idx) {
    leaves <- leaf_sets[[idx]]
    if (!length(leaves)) next
    if (any(assigned[leaves])) next
    if (is.finite(norms[idx]) && norms[idx] < threshold) {
      post_beta[leaves] <- mean(post_beta[leaves])
      assigned[leaves] <- TRUE
      triggered[[length(triggered) + 1L]] <- data.frame(
        node = nodes[idx],
        depth = depths[idx],
        n_leaves = length(leaves),
        penalty_group_index = idx,
        group_norm = norms[idx],
        leaves = paste(leaves, collapse = " "),
        stringsAsFactors = FALSE
      )
    }
  }

  triggered_df <- if (length(triggered)) {
    do.call(rbind, triggered)
  } else {
    data.frame(
      node = integer(),
      depth = integer(),
      n_leaves = integer(),
      penalty_group_index = integer(),
      group_norm = numeric(),
      leaves = character()
    )
  }

  list(beta = post_beta, triggered = triggered_df)
}

rs_thresholded_beta_for_rand <- function(
  beta,
  tree_df,
  X_train,
  y_train,
  method,
  gamma,
  model = "article_aligned",
  normalize_rows = FALSE,
  mu = NULL,
  threshold = 1e-4,
  norm_mode = "scaled"
) {
  if (!method %in% c("RS-DL2", "RS-CL2")) {
    return(list(beta = beta, diagnostics = NULL, triggered = NULL))
  }
  if (!identical(model, "article_aligned")) {
    stop("RS penalty-group thresholded Rand currently supports model = 'article_aligned' only.", call. = FALSE)
  }
  if (!norm_mode %in% c("scaled", "unscaled")) {
    stop("norm_mode must be either 'scaled' or 'unscaled'.", call. = FALSE)
  }
  if (!is.finite(gamma) || gamma <= 0) {
    stop("gamma must be positive and finite.", call. = FALSE)
  }
  if (is.null(mu)) mu <- 1e-2

  X_use <- as.matrix(X_train)
  if (isTRUE(normalize_rows)) X_use <- rs_normalize_rows(X_use)

  article <- rs_article_penalty_matrices(tree_df = tree_df, p = ncol(X_use), methods = method)
  key <- rs_threshold_method_key(method)
  penalty <- article$penalties[[key]]
  if (is.null(penalty)) stop("Cannot find RS penalty matrix for method: ", method, call. = FALSE)

  coef_fit <- rs_threshold_spg_group_one_gamma(
    y = y_train,
    X = X_use %*% article$A,
    gamma = gamma,
    C = penalty$C,
    CNorm = penalty$CNorm,
    g_idx = penalty$g_idx,
    mu = mu
  )

  unscaled_norms <- rs_threshold_group_norms(penalty$C, coef_fit$coef, penalty$g_idx)
  scaled_norms <- gamma * unscaled_norms
  norms <- if (identical(norm_mode, "scaled")) scaled_norms else unscaled_norms

  post <- rs_threshold_apply_top_down(
    beta = beta,
    tree_df = tree_df,
    norms = norms,
    threshold = threshold,
    p = length(beta)
  )

  diagnostics <- data.frame(
    method = method,
    gamma = gamma,
    threshold = threshold,
    norm_mode = norm_mode,
    spg_iter = coef_fit$iter,
    min_unscaled_group_norm = min(unscaled_norms),
    min_scaled_group_norm = min(scaled_norms),
    n_thresholded_unscaled_penalty_groups = sum(unscaled_norms < threshold),
    n_thresholded_scaled_penalty_groups = sum(scaled_norms < threshold),
    n_thresholded_penalty_groups = sum(norms < threshold),
    n_topdown_aggregated_groups = nrow(post$triggered),
    stringsAsFactors = FALSE
  )

  list(beta = post$beta, diagnostics = diagnostics, triggered = post$triggered)
}
