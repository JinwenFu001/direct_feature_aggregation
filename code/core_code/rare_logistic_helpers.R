# Logistic RARE helpers for tree_df-based simulations.
# These mirror the original reproduction code, but use explicit namespaces and
# stable treeFA logistic loss helpers so cluster scripts can source them safely.

check_rare_logistic_packages <- function() {
  required <- c("treeFA", "glmnet", "Matrix")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      ". Install them before running logistic RARE.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rare_log1pexp <- function(x) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  pos <- x > 0
  out[pos] <- x[pos] + log1p(exp(-x[pos]))
  out[!pos] <- log1p(exp(x[!pos]))
  out
}

rare_negtv_lglkh <- function(Y, X, beta) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  if (length(Y) != nrow(X)) stop("length(Y) must equal nrow(X).", call. = FALSE)
  if (any(Y < 0 | Y > 1)) stop("Y must contain values in [0, 1].", call. = FALSE)

  beta <- as.matrix(beta)
  if (nrow(beta) != ncol(X)) stop("beta must have nrow equal to ncol(X).", call. = FALSE)

  eta <- X %*% beta
  log_terms <- matrix(.rare_log1pexp(as.vector(eta)), nrow = nrow(eta))
  unname(colMeans(log_terms - Y * eta))
}

df_to_A <- function(tree_df, p = NULL, sparse = FALSE) {
  tree_df <- as.data.frame(tree_df)
  if (is.null(p)) {
    p <- treeFA::find_p(tree_df)
  }
  p <- as.integer(p)
  if (!is.finite(p) || p <= 0) stop("p must be a positive integer.", call. = FALSE)

  tree_result <- treeFA::gather_leaf_nodes_per_non_leaf(tree_df)$df_result
  tree_result <- tree_result[order(tree_result$depth, decreasing = TRUE), , drop = FALSE]

  A <- diag(p)
  for (i in seq_len(nrow(tree_result))) {
    if (!is.na(tree_result$weight[i]) && tree_result$weight[i] != 0) {
      new_col <- rep(0, p)
      leaves <- tree_result$leaves[[i]]
      if (any(leaves < 1 | leaves > p)) {
        stop("tree_df contains leaves outside 1:p.", call. = FALSE)
      }
      new_col[leaves] <- 1
      A <- cbind(A, new_col)
    }
  }

  if (isTRUE(sparse)) {
    Matrix::Matrix(A, sparse = TRUE)
  } else {
    A
  }
}

rare_logistic_penalty_factor <- function(tree_df, A, p = nrow(A)) {
  tree_result <- treeFA::gather_leaf_nodes_per_non_leaf(tree_df)$df_result
  tree_result <- tree_result[order(tree_result$depth, decreasing = TRUE), , drop = FALSE]
  nonzero_rows <- which(!is.na(tree_result$weight) & tree_result$weight != 0)
  nonzero_nodes <- tree_result$node[nonzero_rows]

  penalty.factor <- rep(1, ncol(A))
  root_nodes <- tree_result$node[tree_result$depth == 0 & !is.na(tree_result$weight) & tree_result$weight != 0]
  if (length(root_nodes)) {
    root_cols <- p + match(root_nodes, nonzero_nodes)
    root_cols <- root_cols[root_cols <= ncol(A)]
    penalty.factor[root_cols] <- 0
  }

  if ((ncol(A) - p) != length(nonzero_nodes)) {
    # This should not happen when A is produced by df_to_A(), but the check helps
    # catch externally supplied matrices with a different tree-column convention.
    warning(
      "A has ", ncol(A) - p, " internal columns, but tree_df has ",
      length(nonzero_nodes), " nonzero-weight internal nodes.",
      call. = FALSE
    )
  }

  penalty.factor
}

rarefit.logistic <- function(
  y,
  X,
  tree_df = NULL,
  A = NULL,
  intercept = FALSE,
  lambda = NULL,
  nlam = 50,
  lam.min.ratio = 1e-4,
  eps = 1e-5,
  maxite = 1e6,
  sparse_A = TRUE
) {
  check_rare_logistic_packages()

  X <- as.matrix(X)
  y <- as.numeric(y)
  if (length(y) != nrow(X)) {
    stop("length(y) must equal nrow(X).", call. = FALSE)
  }
  if (!all(y %in% c(0, 1))) {
    stop("y must be binary and encoded as 0/1.", call. = FALSE)
  }

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(A)) {
    if (is.null(tree_df)) stop("Either tree_df or A must be supplied.", call. = FALSE)
    A <- df_to_A(tree_df, p = p, sparse = sparse_A)
  } else if (isTRUE(sparse_A) && !inherits(A, "sparseMatrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  } else if (!isTRUE(sparse_A)) {
    A <- as.matrix(A)
  }

  if (nrow(A) != p) {
    stop("A must have nrow equal to ncol(X).", call. = FALSE)
  }

  if (!is.null(tree_df)) {
    penalty.factor <- rare_logistic_penalty_factor(tree_df, A = A, p = p)
  } else {
    penalty.factor <- rep(1, ncol(A))
  }

  XA <- X %*% A
  if (is.null(lambda)) {
    lambda.max <- max(abs(as.numeric(Matrix::crossprod(XA, y - 0.5)))) / n
    if (!is.finite(lambda.max) || lambda.max <= 0) lambda.max <- 1
    lambda <- lambda.max * exp(seq(0, log(lam.min.ratio), length.out = nlam))
  } else {
    lambda <- as.numeric(lambda)
    if (any(!is.finite(lambda)) || any(lambda < 0)) {
      stop("lambda must contain finite non-negative values.", call. = FALSE)
    }
    nlam <- length(lambda)
  }

  fit <- glmnet::glmnet(
    x = XA,
    y = y,
    family = "binomial",
    lambda = lambda,
    standardize = FALSE,
    intercept = intercept,
    penalty.factor = penalty.factor,
    thresh = eps,
    maxit = maxite
  )

  gamma <- as.matrix(fit$beta)
  beta <- as.matrix(A %*% gamma)
  beta0 <- if (isTRUE(intercept)) as.numeric(fit$a0) else numeric(ncol(beta))

  list(
    beta0 = beta0,
    beta = beta,
    gamma = gamma,
    lambda = fit$lambda,
    A = A,
    penalty.factor = penalty.factor,
    intercept = intercept,
    glmnet.fit = fit
  )
}

cv.rare.logistic <- function(
  Y,
  X,
  tree_df,
  folds = 5,
  seqc = NULL,
  neg.ind = NULL,
  thresh = 1e-5,
  intercept = FALSE,
  Mmratio = 1e4,
  verbose = FALSE
) {
  check_rare_logistic_packages()
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  if (length(Y) != nrow(X)) stop("length(Y) must equal nrow(X).", call. = FALSE)
  if (!all(Y %in% c(0, 1))) stop("Y must be binary and encoded as 0/1.", call. = FALSE)

  # These checks are local to CV. The shared rarefit.logistic() and its
  # iteration limit are unchanged, including for non-CV simulations.
  usable_columns <- function(fit) {
    jerr <- fit$glmnet.fit$jerr
    if (length(jerr) != 1L || !is.finite(jerr)) {
      stop("RARE CV: glmnet returned an invalid convergence status.", call. = FALSE)
    }
    if (jerr != 0 && !(jerr < 0 && jerr > -10000)) {
      stop("RARE CV: glmnet failed with jerr = ", jerr,
           "; this is not an iteration-limit failure.", call. = FALSE)
    }
    columns <- seq_along(fit$lambda)
    # At jerr = -k, only the points before k converged. A placeholder at
    # jerr = -1 must never be treated as a fitted model.
    if (jerr < 0) columns <- columns[columns < -jerr]
    if (!length(columns)) return(integer())
    if (ncol(fit$beta) != length(fit$lambda) ||
        ncol(fit$gamma) != length(fit$lambda) ||
        length(fit$beta0) != length(fit$lambda)) {
      stop("RARE CV: coefficient columns do not match returned lambda values.",
           call. = FALSE)
    }
    finite <- is.finite(fit$lambda[columns]) & fit$lambda[columns] >= 0 &
      is.finite(fit$beta0[columns]) &
      colSums(!is.finite(fit$beta[, columns, drop = FALSE])) == 0L &
      colSums(!is.finite(fit$gamma[, columns, drop = FALSE])) == 0L
    columns[finite]
  }
  fit_diagnostics <- function(fit, columns) {
    list(jerr = fit$glmnet.fit$jerr, returned.lambda = fit$lambda,
         usable.lambda = fit$lambda[columns])
  }

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

  base.diagnostics <- NULL
  if (is.null(seqc)) {
    base_fit <- rarefit.logistic(
      y = Y,
      X = X,
      tree_df = tree_df,
      intercept = intercept,
      nlam = 50,
      lam.min.ratio = 1 / Mmratio,
      eps = thresh
    )
    base_columns <- usable_columns(base_fit)
    if (!length(base_columns)) {
      stop("RARE CV: the initial full-data path has no usable lambda.", call. = FALSE)
    }
    # Preserve the original rule: obtain the default grid from the full-data
    # path, then compare its usable lambda values across all CV folds.
    seqc <- base_fit$lambda[base_columns]
    base.diagnostics <- fit_diagnostics(base_fit, base_columns)
  }
  seqc <- sort(as.numeric(seqc), decreasing = TRUE)
  if (any(!is.finite(seqc)) || any(seqc < 0) || !length(seqc)) {
    stop("seqc must contain finite non-negative values.", call. = FALSE)
  }

  vals.mat <- matrix(NA_real_, nrow = folds, ncol = length(seqc))
  colnames(vals.mat) <- seqc
  fold.diagnostics <- vector("list", folds)

  for (i in seq_len(folds)) {
    if (isTRUE(verbose)) message("fold ", i, "/", folds)
    train_idx <- which(index != i)
    test_idx <- which(index == i)

    fit <- rarefit.logistic(
      y = Y[train_idx],
      X = X[train_idx, , drop = FALSE],
      tree_df = tree_df,
      intercept = intercept,
      lambda = seqc,
      eps = thresh
    )

    # Match by lambda value, never by the number of returned columns. Missing
    # points remain NA and are ineligible even if other folds fitted them.
    usable <- usable_columns(fit)
    fold.diagnostics[[i]] <- fit_diagnostics(fit, usable)
    if (anyNA(match(fit$lambda[usable], seqc))) {
      stop("RARE CV: glmnet returned a lambda outside the requested grid.", call. = FALSE)
    }
    fit_columns <- match(seqc, fit$lambda[usable])
    present <- which(!is.na(fit_columns))
    if (!length(present)) next
    losses <- rare_negtv_lglkh(
      Y[test_idx],
      X[test_idx, , drop = FALSE],
      beta = fit$beta[, usable[fit_columns[present]], drop = FALSE]
    )
    if (length(losses) != length(present)) {
      stop("RARE CV: loss length does not match the returned lambda values.", call. = FALSE)
    }
    good <- is.finite(losses)
    vals.mat[i, present[good]] <- losses[good]
  }

  eligible <- colSums(is.finite(vals.mat)) == folds
  if (!any(eligible)) {
    stop("RARE CV: no lambda has a usable validation loss in every fold. ",
         "No model was selected; do not average over missing folds.", call. = FALSE)
  }
  if (any(!eligible)) {
    warning("RARE CV: excluded ", sum(!eligible), "/", length(seqc),
            " lambda values without a usable loss in every fold.", call. = FALSE)
  }
  vals.vec <- rep(Inf, length(seqc))
  names(vals.vec) <- colnames(vals.mat)
  vals.vec[eligible] <- colMeans(vals.mat[, eligible, drop = FALSE])
  selected.param <- seqc[which.min(vals.vec)]
  final_fit <- rarefit.logistic(
    y = Y,
    X = X,
    tree_df = tree_df,
    intercept = intercept,
    lambda = selected.param,
    eps = thresh
  )
  final_columns <- usable_columns(final_fit)
  if (length(final_columns) != 1L || final_fit$lambda[final_columns] != selected.param) {
    stop("RARE CV: the selected lambda did not yield a valid full-data refit.", call. = FALSE)
  }

  list(
    selected.param = selected.param,
    beta = final_fit,
    valid.error = vals.vec,
    lambda = seqc,
    eligible = eligible,
    excluded.lambda = seqc[!eligible],
    fold.error = vals.mat,
    fold.diagnostics = fold.diagnostics,
    base.diagnostics = base.diagnostics
  )
}


log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}

# Use the new treeFA tree API to assign normalized group-size weights and
# construct the RARE expansion. The full tree has one unpenalized root, last in A.
make_tree_objects <- function(tree_df, p, weight.order) {
  tree_df <- as.data.frame(tree_df)
  if (!identical(as.integer(treeFA::tree_leaf_ids(tree_df)), seq_len(p))) {
    stop("The tree leaves must be numbered 1:ncol(X).", call. = FALSE)
  }
  tree_df$weight <- 0
  groups <- treeFA::gather_leaf_nodes_per_non_leaf(tree_df)$df_result
  if (!nrow(groups) || sum(groups$depth == 0) != 1L) {
    stop("This analysis requires a full tree with one internal root.", call. = FALSE)
  }
  raw_weights <- lengths(groups$leaves)^weight.order
  tree_df$weight[match(groups$node, tree_df$node)] <- raw_weights / mean(raw_weights)
  groups <- groups[order(groups$depth, decreasing = TRUE), , drop = FALSE]
  leaf_sets <- c(lapply(seq_len(p), identity), groups$leaves)
  A <- Matrix::sparseMatrix(
    i = unlist(leaf_sets, use.names = FALSE),
    j = rep(seq_along(leaf_sets), lengths(leaf_sets)),
    x = 1, dims = c(p, length(leaf_sets))
  )
  A <- methods::as(A, "dgCMatrix")
  if (!all(as.numeric(A[, ncol(A)]) == 1)) {
    stop("The last RARE expansion column must be the root.", call. = FALSE)
  }
  list(tree_df = tree_df, A = A)
}

# Pure Gaussian RARE: use the official package, alpha=1, without an intercept.
# At alpha=1 its glmnet path is in descending lambda order. rarefit reports
# the supplied grid even when glmnet returns only a prefix, so align explicitly.
fit_rare_path <- function(y, X, A, lambda) {
  lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  if (!length(lambda) || any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("Invalid RARE lambda grid.", call. = FALSE)
  }
  warnings <- character()
  fit <- withCallingHandlers(
    rare::rarefit(y = y, X = X, A = A, intercept = FALSE, alpha = 1,
                  lambda = lambda, rho = 1e-2, eps1 = 1e-6, eps2 = 1e-5,
                  maxite = 1000000L),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
    }
  )
  B <- as.matrix(fit$beta[[1L]])
  G <- as.matrix(fit$gamma[[1L]])
  if (nrow(B) != ncol(X) || nrow(G) != ncol(A) || ncol(G) != ncol(B) ||
      ncol(B) > length(lambda)) {
    stop("RARE returned inconsistent coefficient dimensions.", call. = FALSE)
  }
  count <- ncol(B)
  matches <- regmatches(warnings, regexec("error code[[:space:]]+(-?[0-9]+)", warnings))
  codes <- vapply(matches[lengths(matches) == 2L], function(x) as.integer(x[2L]), integer(1))
  if (any(codes > 0L | codes <= -10000L)) {
    stop("RARE/glmnet failed for a reason other than its iteration limit: ",
         paste(codes, collapse = ", "), call. = FALSE)
  }
  failed_at <- -codes[codes < 0L & codes > -10000L]
  if (length(failed_at)) count <- min(count, min(failed_at) - 1L)
  # A first-lambda failure can return a dummy zero column; it is never eligible.
  B_aligned <- matrix(NA_real_, nrow = ncol(X), ncol = length(lambda))
  available <- rep(FALSE, length(lambda))
  if (count > 0L) {
    columns <- seq_len(count)
    good <- colSums(!is.finite(B[, columns, drop = FALSE])) == 0L &
      colSums(!is.finite(G[, columns, drop = FALSE])) == 0L
    columns <- columns[good]
    B_aligned[, columns] <- B[, columns, drop = FALSE]
    available[columns] <- TRUE
  }
  list(beta = B_aligned, lambda = lambda, available = available,
       diagnostics = list(warnings = warnings, glmnet_codes = codes,
                          available = available, maxite = 1000000L))
}

cv_rare <- function(y, X, A, foldid) {
  # This is rare::rarefit's default 50-point grid, computed on training data.
  lambda_max <- max(abs(as.numeric(crossprod(X, y)))) / nrow(X)
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("RARE lambda grid has no positive upper endpoint.", call. = FALSE)
  }
  lambda <- lambda_max * exp(seq(0, log(1e-4), length.out = 50L))
  kfold <- max(foldid)
  errors <- matrix(NA_real_, kfold, length(lambda), dimnames = list(NULL, lambda))
  diagnostics <- vector("list", kfold)
  for (k in seq_len(kfold)) {
    log_msg("RARE CV fold ", k, "/", kfold)
    train <- foldid != k
    fit <- fit_rare_path(y[train], X[train, , drop = FALSE], A, lambda)
    columns <- which(fit$available)
    if (length(columns)) {
      predictions <- X[!train, , drop = FALSE] %*% fit$beta[, columns, drop = FALSE]
      losses <- colMeans((y[!train] - predictions)^2)
      finite <- is.finite(losses)
      errors[k, columns[finite]] <- losses[finite]
    }
    diagnostics[[k]] <- fit$diagnostics
  }
  eligible <- colSums(is.finite(errors)) == kfold
  if (!any(eligible)) stop("No RARE lambda is usable in every CV fold.", call. = FALSE)
  if (any(!eligible)) {
    warning("RARE CV excluded ", sum(!eligible), " of ", length(lambda),
            " lambda values; no iteration limit was increased.", call. = FALSE)
  }
  scores <- rep(Inf, length(lambda))
  scores[eligible] <- colMeans(errors[, eligible, drop = FALSE])
  selected <- which.min(scores)
  final <- fit_rare_path(y, X, A, lambda[selected])
  if (!isTRUE(final$available[1L])) {
    stop("Selected RARE lambda did not produce a usable full-training-data fit.", call. = FALSE)
  }
  list(beta = as.numeric(final$beta[, 1L]), selected.param = lambda[selected],
       lambda = lambda, valid.error = scores, eligible = eligible,
       fold.error = errors, fold.diagnostics = diagnostics,
       refit.diagnostics = final$diagnostics)
}