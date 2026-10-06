# Joint aggregation/selection is intentionally independent of the legacy solvers.
.selection_check_scalar <- function(x, name, lower = 0, integer = FALSE,
                                    strictly_positive = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || (strictly_positive && x <= 0) ||
      (integer && (x != floor(x) || x > .Machine$integer.max))) {
    stop(name, " must be a finite ", if (integer) "integer " else "numeric ",
         "scalar ", if (strictly_positive) "> 0." else paste0(">= ", lower, "."),
         call. = FALSE)
  }
  invisible(TRUE)
}

.selection_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  invisible(TRUE)
}

.selection_prepare_tree <- function(tree_result, p) {
  if (is.data.frame(tree_result)) {
    tree <- tree_result
    if (!all(c("node", "parent", "weight") %in% names(tree)) || !nrow(tree))
      stop("tree_df must have nonempty node, parent and weight columns.", call. = FALSE)
    ids <- tree$node
    parents <- tree$parent
    if (is.logical(parents) && all(is.na(parents))) parents <- as.numeric(parents)
    weights <- tree$weight
    valid_ids <- function(x) is.numeric(x) && !is.complex(x) &&
      all(is.finite(x)) && all(x >= 1 & x <= .Machine$integer.max & x == floor(x))
    if (!valid_ids(ids) || anyDuplicated(ids))
      stop("Tree node IDs must be unique positive integers.", call. = FALSE)
    if (!is.numeric(parents) || is.complex(parents) ||
        !valid_ids(parents[!is.na(parents)]))
      stop("Tree parent IDs must be positive integers or NA.", call. = FALSE)
    if (!is.numeric(weights) || is.complex(weights) ||
        any(!is.finite(weights)) || any(weights < 0))
      stop("Tree weights must be finite and non-negative.", call. = FALSE)
    root <- which(is.na(parents))
    parent_index <- match(parents, ids)
    if (length(root) != 1L || any(is.na(parent_index[-root])))
      stop("tree_df must have one root and valid parent IDs.", call. = FALSE)
    children <- split(which(!is.na(parent_index)), parent_index[!is.na(parent_index)])
    # Traverse once, then build leaf sets in reverse (children before parents).
    order <- integer(length(ids))
    order[1L] <- root
    head <- tail <- 1L
    while (head <= tail) {
      kids <- children[[as.character(order[head])]]
      if (length(kids)) {
        if (tail + length(kids) > length(ids))
          stop("Cycle detected in tree_df.", call. = FALSE)
        order[tail + seq_along(kids)] <- kids
        tail <- tail + length(kids)
      }
      head <- head + 1L
    }
    if (tail != length(ids))
      stop("tree_df must be connected and contain no cycles.", call. = FALSE)
    leaf_rows <- which(!seq_along(ids) %in% parent_index)
    if (!identical(sort(as.integer(ids[leaf_rows])), seq_len(p)))
      stop("Tree leaves must be indexed as 1:p to match the columns of X.", call. = FALSE)
    descendants <- vector("list", length(ids))
    groups <- list()
    group_weights <- numeric()
    for (i in rev(order)) {
      kids <- children[[as.character(i)]]
      if (!length(kids)) {
        descendants[[i]] <- as.integer(ids[i])
      } else {
        descendants[[i]] <- unlist(descendants[kids], use.names = FALSE)
        groups[[length(groups) + 1L]] <- descendants[[i]]
        group_weights <- c(group_weights, weights[i])
      }
    }
    return(list(groups = groups, weights = group_weights, p = p))
  }

  if (is.list(tree_result) && !is.null(tree_result$list_result))
    tree_result <- tree_result$list_result
  if (!is.list(tree_result) || !is.list(tree_result$groups) ||
      !is.list(tree_result$weights) ||
      length(tree_result$groups) != length(tree_result$weights))
    stop("tree_result must be a tree data frame or a gathered tree list.", call. = FALSE)

  layers <- tree_result$groups
  weights <- tree_result$weights
  # Assign each coordinate its nearest coarser ancestor. This validates laminar
  # structure without checking all pairs of groups or building dense operators.
  owner <- integer(p)
  seen_layer <- integer(p)
  group_id <- 0L
  for (level in rev(seq_along(layers))) {
    layer <- layers[[level]]
    w <- weights[[level]]
    if (!is.list(layer) || !is.numeric(w) || is.complex(w) ||
        length(layer) != length(w) || any(!is.finite(w)) || any(w < 0))
      stop("Tree groups and finite non-negative weights must have matching lengths.", call. = FALSE)
    for (idx in layer) {
      if (!is.numeric(idx) || is.complex(idx) || !length(idx) ||
          any(!is.finite(idx)) || any(idx != floor(idx)) ||
          any(idx < 1 | idx > p) || anyDuplicated(idx))
        stop("Tree group indices must be unique integers in 1:p.", call. = FALSE)
      idx <- as.integer(idx)
      if (any(seen_layer[idx] == level))
        stop("Tree groups in the same layer must not overlap.", call. = FALSE)
      if (length(unique(owner[idx])) != 1L)
        stop("Tree groups must be nested and ordered from children to parents.", call. = FALSE)
      group_id <- group_id + 1L
      owner[idx] <- group_id
      seen_layer[idx] <- level
    }
  }
  groups <- do.call(c, unname(layers))
  if (is.null(groups)) groups <- list()
  list(groups = lapply(groups, as.integer),
       weights = as.numeric(unlist(weights, use.names = FALSE)), p = p)
}

.selection_prepare_model <- function(Y, X, tree_result, family = "gaussian",
                                     intercept = FALSE, ridge_param = 0,
                                     step_size = NULL) {
  if (!family %in% c("gaussian", "binomial"))
    stop("Unknown selection family.", call. = FALSE)
  .selection_check_flag(intercept, "intercept")
  .selection_check_scalar(ridge_param, "ridge_param")
  X <- as.matrix(X)
  if (!is.numeric(X) || is.complex(X) || !is.numeric(Y) || is.complex(Y) ||
      nrow(X) < 1L || ncol(X) < 1L || length(Y) != nrow(X) ||
      any(!is.finite(X)) || any(!is.finite(Y)))
    stop("X and Y must be finite numeric data with nrow(X) = length(Y) and positive dimensions.", call. = FALSE)
  storage.mode(X) <- "double"
  Y <- as.numeric(Y)
  if (family == "binomial") {
    if (!all(Y %in% c(0, 1)))
      stop("Y must be binary and encoded as 0/1.", call. = FALSE)
    if (intercept && length(unique(Y)) < 2L)
      stop("A single-class binomial response has no finite unpenalized intercept.", call. = FALSE)
  }
  tree <- .selection_prepare_tree(tree_result, ncol(X))
  X_mean <- rep(0, ncol(X))
  Y_mean <- 0
  if (family == "gaussian" && intercept) {
    X_mean <- colMeans(X)
    Y_mean <- mean(Y)
    X <- sweep(X, 2L, X_mean, "-")
    Y <- Y - Y_mean
  }
  solver_intercept <- family == "binomial" && intercept
  A <- if (solver_intercept) cbind(1, X) else X
  # Base R's spectral norm uses an SVD; unlike a short power iteration this is
  # not an unchecked underestimate. A small margin protects the strict bound.
  spectral <- norm(A, type = "2")
  L <- ((if (family == "gaussian") 1 else 0.25) * spectral^2 + ridge_param) / nrow(X)
  L <- L * (1 + 1e-12)
  if (!is.finite(L))
    stop("The Lipschitz bound overflowed; rescale the supplied data.", call. = FALSE)
  if (is.null(step_size)) step_size <- if (L > 0) 1 / L else 1
  .selection_check_scalar(step_size, "step_size", strictly_positive = TRUE)
  if (L > 0 && step_size * L >= 2)
    stop("step_size must satisfy 0 < step_size < 2 / L.", call. = FALSE)
  list(Y = Y, X = X, tree = tree, family = family,
       family_code = if (family == "gaussian") 0L else 1L,
       intercept = intercept, solver_intercept = solver_intercept,
       ridge_param = ridge_param, step_size = step_size, lipschitz = L,
       p = ncol(X), X_mean = X_mean, Y_mean = Y_mean,
       feature_names = colnames(X))
}

.selection_fit_prepared <- function(prepared, lambda, tau, abs_tol = 1e-7,
                                    rel_tol = 1e-6, max_iter = 10000L,
                                    init_z = NULL, keep_history = FALSE,
                                    warn = TRUE) {
  .selection_check_scalar(lambda, "lambda")
  .selection_check_scalar(tau, "tau")
  .selection_check_scalar(abs_tol, "abs_tol")
  .selection_check_scalar(rel_tol, "rel_tol")
  if (abs_tol == 0 && rel_tol == 0)
    stop("At least one of abs_tol and rel_tol must be positive.", call. = FALSE)
  .selection_check_scalar(max_iter, "max_iter", integer = TRUE, strictly_positive = TRUE)
  .selection_check_flag(keep_history, "keep_history")
  dimension <- prepared$p + as.integer(prepared$solver_intercept)
  if (is.null(init_z)) init_z <- numeric(dimension)
  if (!is.numeric(init_z) || is.complex(init_z) || length(init_z) != dimension ||
      any(!is.finite(init_z)))
    stop("init_z must be a finite numeric vector of length ", dimension, ".", call. = FALSE)
  fit <- selection_fit_cpp(
    Y = prepared$Y, X = prepared$X,
    groups = prepared$tree$groups, weights = prepared$tree$weights,
    lambda = lambda, tau = tau, family = prepared$family_code,
    intercept = prepared$solver_intercept, ridge_param = prepared$ridge_param,
    step_size = prepared$step_size, abs_tol = abs_tol, rel_tol = rel_tol,
    max_iter = as.integer(max_iter), init_z = as.numeric(init_z),
    keep_history = keep_history
  )
  fit$beta <- as.numeric(fit$beta)
  names(fit$beta) <- prepared$feature_names
  if (prepared$family == "gaussian" && prepared$intercept)
    fit$beta0 <- as.numeric(prepared$Y_mean - sum(prepared$X_mean * fit$beta))
  if (!is.finite(fit$beta0) || any(!is.finite(fit$beta)) || !is.finite(fit$objective)) {
    fit$converged <- FALSE
    fit$status <- "nonfinite"
  }
  fit$state <- list(z = as.numeric(fit$z), step_size = prepared$step_size)
  fit$z <- NULL
  fit$lambda <- lambda
  fit$tau <- tau
  fit$ridge_param <- prepared$ridge_param
  fit$family <- prepared$family
  fit$intercept <- prepared$intercept
  fit$step_size <- prepared$step_size
  fit$lipschitz <- prepared$lipschitz
  fit$n_nonzero <- sum(fit$beta != 0)
  if (warn && !isTRUE(fit$converged))
    warning("Davis-Yin selection fit did not converge (", fit$status,
            "); inspect residual_scaled and consider increasing max_iter.", call. = FALSE)
  fit
}

dy_prox_simple_linear <- function(Y, X, tree_result, lambda, tau,
                                  intercept = FALSE, ridge_param = 0,
                                  step_size = NULL, abs_tol = 1e-7,
                                  rel_tol = 1e-6, max_iter = 10000L,
                                  init_z = NULL, keep_history = FALSE) {
  prepared <- .selection_prepare_model(Y, X, tree_result, "gaussian", intercept,
                                       ridge_param, step_size)
  .selection_fit_prepared(prepared, lambda, tau, abs_tol, rel_tol, max_iter,
                          init_z, keep_history)
}

dy_prox_simple_logistic <- function(Y, X, tree_result, lambda, tau,
                                    intercept = FALSE, ridge_param = 0,
                                    step_size = NULL, abs_tol = 1e-7,
                                    rel_tol = 1e-6, max_iter = 10000L,
                                    init_z = NULL, keep_history = FALSE) {
  prepared <- .selection_prepare_model(Y, X, tree_result, "binomial", intercept,
                                       ridge_param, step_size)
  .selection_fit_prepared(prepared, lambda, tau, abs_tol, rel_tol, max_iter,
                          init_z, keep_history)
}
