# Exact misspecified-tree data generator for treeFA / RARE / RS simulations.
# Designed for outer = 10, inner = 40, p around 200, q around 10.
# The generator returns one fixed scenario per outer replicate and one noise draw per inner replicate.

make_balanced_groups <- function(p = 200L, q = 10L) {
  stopifnot(p >= q, q >= 2L)
  sizes <- rep(p %/% q, q)
  if (p %% q > 0L) sizes[seq_len(p %% q)] <- sizes[seq_len(p %% q)] + 1L
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(seq_len(q), function(i) seq.int(starts[i], ends[i]))
}

partition_to_group_index <- function(groups, p) {
  out <- integer(p)
  for (g in seq_along(groups)) out[groups[[g]]] <- g
  if (any(out == 0L)) stop("Some leaves are not covered by the partition.", call. = FALSE)
  out
}

H_from_groups <- function(groups, p) {
  H <- matrix(0, nrow = p, ncol = length(groups))
  for (g in seq_along(groups)) H[groups[[g]], g] <- 1
  H
}

projection_residual_scale <- function(X, H, v, standardized = TRUE, tol = 1e-10) {
  Z <- X %*% H
  Xv <- drop(X %*% v)
  qrz <- qr(Z, tol = tol)
  resid <- qr.resid(qrz, Xv)
  val <- sqrt(sum(resid^2))
  if (standardized) val <- val / sqrt(nrow(X))
  as.numeric(val)
}

gamma_beta_scale <- function(X, H, beta, standardized = TRUE, tol = 1e-10) {
  Z <- X %*% H
  xb <- drop(X %*% beta)
  qrz <- qr(Z, tol = tol)
  resid <- qr.resid(qrz, xb)
  val <- sqrt(sum(resid^2))
  if (standardized) val <- val / sqrt(nrow(X))
  as.numeric(val)
}

# Build a tree_df object directly from a top-level partition.
# Leaves are 1:p. Internal node IDs are contiguous p+1, ..., nrow(tree_df).
# If add_super_root = TRUE, the tree has one root with q children.
# If add_super_root = FALSE, the output is a q-root forest.
make_balanced_tree_df_from_groups <- function(
    groups,
    p,
    weight.order = -1,
    add_super_root = TRUE,
    jitter = TRUE,
    jitter_max = 1L,
    random_leaf_order = TRUE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  q <- length(groups)
  groups <- lapply(groups, function(x) sort(as.integer(x)))

  all_leaves <- sort(unlist(groups, use.names = FALSE))
  if (!identical(all_leaves, seq_len(p))) {
    stop("groups must form a partition of leaves 1:p.", call. = FALSE)
  }

  node_vec <- seq_len(p)
  parent_vec <- rep(NA_integer_, p)
  name_vec <- paste0("leaf", seq_len(p))
  weight_vec <- rep(0, p)

  next_id <- p + 1L
  n_internal_below_top <- sum(vapply(groups, function(g) max(length(g) - 1L, 0L), integer(1)))
  super_id <- if (add_super_root) p + n_internal_below_top + 1L else NA_integer_

  add_internal <- function(parent) {
    id <- next_id
    next_id <<- next_id + 1L
    node_vec <<- c(node_vec, id)
    parent_vec <<- c(parent_vec, as.integer(parent))
    name_vec <<- c(name_vec, paste0("latent", id))
    weight_vec <<- c(weight_vec, NA_real_)
    id
  }

  choose_left <- function(leaves) {
    m <- length(leaves)
    left_n <- floor(m / 2)
    if (jitter && m >= 6L && jitter_max > 0L) {
      shift <- sample(seq.int(-jitter_max, jitter_max), size = 1L)
      left_n <- left_n + shift
    }
    left_n <- max(1L, min(m - 1L, left_n))
    if (random_leaf_order) {
      sort(sample(leaves, size = left_n, replace = FALSE))
    } else {
      leaves[seq_len(left_n)]
    }
  }

  add_subtree <- function(leaves, parent) {
    leaves <- sort(as.integer(leaves))
    if (length(leaves) == 1L) {
      parent_vec[leaves] <<- as.integer(parent)
      return(leaves)
    }
    id <- add_internal(parent)
    left <- choose_left(leaves)
    right <- sort(setdiff(leaves, left))
    add_subtree(left, id)
    add_subtree(right, id)
    id
  }

  top_parent <- if (add_super_root) super_id else NA_integer_
  for (g in seq_along(groups)) add_subtree(groups[[g]], top_parent)

  if (add_super_root) {
    node_vec <- c(node_vec, super_id)
    parent_vec <- c(parent_vec, NA_integer_)
    name_vec <- c(name_vec, paste0("latent", super_id))
    weight_vec <- c(weight_vec, NA_real_)
  }

  tree_df <- data.frame(
    node = as.integer(node_vec),
    parent = as.integer(parent_vec),
    name = as.character(name_vec),
    weight = as.numeric(weight_vec),
    stringsAsFactors = FALSE
  )
  tree_df <- tree_df[order(tree_df$node), , drop = FALSE]
  rownames(tree_df) <- NULL

  # Recompute weights as |descendant leaves|^weight.order, normalized over internal nodes.
  children <- split(tree_df$node[!is.na(tree_df$parent)], tree_df$parent[!is.na(tree_df$parent)])
  memo <- new.env(parent = emptyenv())

  leaf_count <- function(node) {
    key <- as.character(node)
    if (exists(key, envir = memo, inherits = FALSE)) return(get(key, envir = memo))
    if (node <= p) {
      val <- 1L
    } else {
      ch <- children[[key]]
      if (is.null(ch)) {
        val <- 0L
      } else {
        val <- sum(vapply(ch, leaf_count, numeric(1)))
      }
    }
    assign(key, val, envir = memo)
    val
  }

  internal <- tree_df$node[tree_df$node > p]
  if (length(internal)) {
    counts <- vapply(internal, leaf_count, numeric(1))
    w <- counts^weight.order
    w <- w / mean(w)
    tree_df$weight[match(internal, tree_df$node)] <- w
  }
  tree_df$weight[tree_df$node <= p] <- 0

  tree_df
}

swap_between_two_groups <- function(groups, baseline, carrier, m) {
  G0 <- groups[[baseline]]
  Gi <- groups[[carrier]]
  if (m < 1L || m >= min(length(G0), length(Gi))) {
    stop("m must satisfy 1 <= m < min(length(baseline group), length(carrier group)).", call. = FALSE)
  }

  S0 <- sample(G0, size = m, replace = FALSE)
  Si <- sample(Gi, size = m, replace = FALSE)

  new_groups <- groups
  new_groups[[baseline]] <- sort(c(setdiff(G0, S0), Si))
  new_groups[[carrier]] <- sort(c(setdiff(Gi, Si), S0))

  list(groups = new_groups, S_baseline = sort(S0), S_carrier = sort(Si))
}

make_exact_misspec_scenario <- function(
    X_cal,
    c_levels,
    p = ncol(X_cal),
    q = 10L,
    q_groups = NULL,
    m_grid = NULL,
    baseline = NULL,
    carriers = NULL,
    weight.order = -1,
    add_super_root = TRUE,
    tree_seed = 1L,
    swap_seed = 1L,
    max_swap_tries = 200L,
    tol = 1e-8
) {
  stopifnot(length(c_levels) >= 1L, abs(c_levels[1]) < tol)
  M <- length(c_levels) - 1L
  if (M > q - 1L) stop("This construction needs length(c_levels) - 1 <= q - 1.", call. = FALSE)

  if (is.null(q_groups)) q_groups <- make_balanced_groups(p = p, q = q)
  q_groups <- lapply(q_groups, function(x) sort(as.integer(x)))

  if (is.null(baseline) || is.null(carriers)) {
    set.seed(swap_seed)
    baseline <- sample(seq_len(q), size = 1L)
    carriers <- sample(setdiff(seq_len(q), baseline), size = M, replace = FALSE)
  }
  if (length(carriers) < M) stop("Need at least M carrier branches.", call. = FALSE)

  if (is.null(m_grid)) {
    min_size <- min(vapply(q_groups, length, integer(1)))
    max_m <- max(1L, floor(min_size / 2L))
    m_grid <- ceiling(seq(1L, max_m, length.out = M))
  }
  m_grid <- as.integer(m_grid)
  if (length(m_grid) < M) stop("m_grid must have length at least length(c_levels)-1.", call. = FALSE)

  partitions <- vector("list", M + 1L)
  partitions[[1L]] <- q_groups

  tree_df_list <- vector("list", M + 1L)
  names(tree_df_list) <- paste0("level_", seq_len(M + 1L) - 1L)
  tree_df_list[[1L]] <- make_balanced_tree_df_from_groups(
    groups = partitions[[1L]], p = p, weight.order = weight.order,
    add_super_root = add_super_root, seed = tree_seed
  )

  beta_star <- numeric(p)
  gamma_info <- data.frame(
    level = seq_len(M + 1L) - 1L,
    target_c = as.numeric(c_levels),
    gamma_hat = NA_real_,
    scale_s = NA_real_,
    alpha = NA_real_,
    baseline = as.integer(baseline),
    carrier = NA_integer_,
    m = NA_integer_
  )

  swap_info <- vector("list", M)

  for (i in seq_len(M)) {
    c_i <- c_levels[i + 1L]
    carrier <- carriers[i]
    m_i <- m_grid[i]

    found <- FALSE
    for (try_id in seq_len(max_swap_tries)) {
      set.seed(swap_seed + 1000L * i + try_id)
      sw <- swap_between_two_groups(q_groups, baseline = baseline, carrier = carrier, m = m_i)
      part_i <- sw$groups
      H_i <- H_from_groups(part_i, p = p)

      u_i <- numeric(p)
      u_i[q_groups[[carrier]]] <- 1
      s_i <- projection_residual_scale(X_cal, H_i, u_i, standardized = TRUE)

      if (is.finite(s_i) && s_i > tol) {
        found <- TRUE
        break
      }
    }

    if (!found) {
      stop("Could not find a nonzero residual scale for level ", i,
           ". Try another X, branch, or m_grid.", call. = FALSE)
    }

    alpha_i <- c_i / s_i
    beta_star <- beta_star + alpha_i * u_i

    partitions[[i + 1L]] <- part_i
    tree_df_list[[i + 1L]] <- make_balanced_tree_df_from_groups(
      groups = part_i, p = p, weight.order = weight.order,
      add_super_root = add_super_root, seed = tree_seed + i
    )

    gamma_info$scale_s[i + 1L] <- s_i
    gamma_info$alpha[i + 1L] <- alpha_i
    gamma_info$carrier[i + 1L] <- carrier
    gamma_info$m[i + 1L] <- m_i
    swap_info[[i]] <- sw
  }

  # Verify exact standardized approximation error under the unique q-cut partition.
  for (i in seq_len(M + 1L)) {
    H_i <- H_from_groups(partitions[[i]], p = p)
    gamma_info$gamma_hat[i] <- gamma_beta_scale(X_cal, H_i, beta_star, standardized = TRUE)
  }

  err <- max(abs(gamma_info$gamma_hat - gamma_info$target_c))
  if (!is.finite(err) || err > 1e-6) {
    stop("Gamma verification failed. Max absolute error = ", signif(err, 4), call. = FALSE)
  }

  names(partitions) <- names(tree_df_list)
  list(
    beta_star = as.numeric(beta_star),
    true_group = partition_to_group_index(q_groups, p = p),
    true_groups = q_groups,
    partitions = partitions,
    tree_df_list = tree_df_list,
    gamma_info = gamma_info,
    swap_info = swap_info,
    baseline = baseline,
    carriers = carriers,
    m_grid = m_grid
  )
}

misspec_outer_inner_id <- function(rep_id, inner_reps = 40L) {
  outer <- ((rep_id - 1L) %/% inner_reps) + 1L
  inner <- ((rep_id - 1L) %% inner_reps) + 1L
  list(outer = outer, inner = inner)
}

simulate_data_misspecified_tree <- function(
    outer_id,
    inner_id,
    c_levels,
    n = 50L,
    n0 = 500L,
    n1 = 50L,
    p = 200L,
    q = 10L,
    ratio = 5,
    poisson_rate = 0.02,
    weight.order = -1,
    q_groups = NULL,
    m_grid = NULL,
    calibrate_on = c("train", "train_valid", "all"),
    add_super_root = TRUE,
    base_seed = 20260706L
) {
  calibrate_on <- match.arg(calibrate_on)
  M <- length(c_levels) - 1L
  if (M > q - 1L) stop("Need at most q-1 positive c-levels.", call. = FALSE)

  n_total <- n + n0 + n1

  scenario_seed <- as.integer(base_seed + 10000L * outer_id)
  response_seed <- as.integer(base_seed + 1000000L + 10000L * outer_id + inner_id)

  set.seed(scenario_seed)
  X <- matrix(stats::rpois(n_total * p, poisson_rate), nrow = n_total, ncol = p)

  train_index <- sample(seq_len(n_total), n)
  valid_index <- sample(setdiff(seq_len(n_total), train_index), n1)
  test_index <- setdiff(seq_len(n_total), c(train_index, valid_index))

  X_cal <- switch(
    calibrate_on,
    train = X[train_index, , drop = FALSE],
    train_valid = X[c(train_index, valid_index), , drop = FALSE],
    all = X
  )

  scenario <- make_exact_misspec_scenario(
    X_cal = X_cal,
    c_levels = c_levels,
    p = p,
    q = q,
    q_groups = q_groups,
    m_grid = m_grid,
    weight.order = weight.order,
    add_super_root = add_super_root,
    tree_seed = scenario_seed + 10L,
    swap_seed = scenario_seed + 20L
  )

  means <- drop(X %*% scenario$beta_star)
  sigma <- sqrt(sum(means^2) / ratio / n_total)

  set.seed(response_seed)
  Y <- means + stats::rnorm(n_total, mean = 0, sd = sigma)

  out <- list(
    outer_id = outer_id,
    inner_id = inner_id,
    c_levels = c_levels,
    p = p,
    q = q,
    n = n,
    n0 = n0,
    n1 = n1,
    X = X,
    Y = as.numeric(Y),
    true.beta = scenario$beta_star,
    A = H_from_groups(scenario$true_groups, p = p),
    true_group = scenario$true_group,
    means = as.numeric(means),
    sigma = sigma,
    train_index = train_index,
    valid_index = valid_index,
    test_index = test_index,
    X_train = X[train_index, , drop = FALSE],
    Y_train = as.numeric(Y[train_index]),
    X_valid = X[valid_index, , drop = FALSE],
    Y_valid = as.numeric(Y[valid_index]),
    X_test = X[test_index, , drop = FALSE],
    Y_test = as.numeric(Y[test_index]),
    tree_df_list = scenario$tree_df_list,
    partitions = scenario$partitions,
    true_groups = scenario$true_groups,
    gamma_info = scenario$gamma_info,
    baseline = scenario$baseline,
    carriers = scenario$carriers,
    m_grid = scenario$m_grid,
    calibrate_on = calibrate_on,
    scenario_seed = scenario_seed,
    response_seed = response_seed
  )

  out
}
