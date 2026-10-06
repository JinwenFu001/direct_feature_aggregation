script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

setup_libpaths <- function(code_dir) {
  code_dir <- normalizePath(code_dir, mustWork = FALSE)
  candidates <- c(
    Sys.getenv("STUDY2_R_LIB", unset = ""),
    Sys.getenv("TREEFA_R_LIB", unset = ""),
    file.path(code_dir, "rlib"),
    file.path(dirname(code_dir), "rlib")
  )
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (length(candidates)) .libPaths(unique(c(candidates, .libPaths())))
  invisible(.libPaths())
}

resolve_cluster_code_dir <- function(code_dir) {
  code_dir <- normalizePath(code_dir, mustWork = TRUE)
  direct_required <- c("rs_matlab_functions_cluster.R", "variant_tree_helpers.R")
  if (all(file.exists(file.path(code_dir, direct_required)))) {
    return(code_dir)
  }

  nested <- file.path(code_dir, "cluster_treeFA_RARE_RS")
  if (dir.exists(nested) && all(file.exists(file.path(nested, direct_required)))) {
    return(normalizePath(nested, mustWork = TRUE))
  }

  stop(
    "Cannot find original cluster RS helper files under STUDY2_CODE_DIR. Expected either ",
    file.path(code_dir, "rs_matlab_functions_cluster.R"),
    " or ",
    file.path(nested, "rs_matlab_functions_cluster.R"),
    ".",
    call. = FALSE
  )
}

source_cluster_code <- function(code_dir) {
  code_dir <- resolve_cluster_code_dir(code_dir)
  Sys.setenv(TREEFA_CLUSTER_CODE_DIR = code_dir)
  source(file.path(code_dir, "rs_matlab_functions_cluster.R"))
  source(file.path(code_dir, "variant_tree_helpers.R"))
  code_dir
}

study2_rs_code_dir <- function(code_dir) {
  code_dir <- resolve_cluster_code_dir(code_dir)
  candidates <- c(file.path(code_dir, "RS-Code"), file.path(code_dir, "RS-code"))
  candidates <- candidates[dir.exists(candidates)]
  if (!length(candidates)) {
    stop("Cannot find RS-Code/ under: ", code_dir, call. = FALSE)
  }
  normalizePath(candidates[1], mustWork = TRUE)
}

check_packages <- function(run_rare = TRUE, run_rs = TRUE) {
  required <- character()
  if (run_rare) required <- c(required, "rare", "Matrix")
  if (run_rs) required <- c(required, "glmnet")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R package(s): ", paste(unique(missing), collapse = ", "),
      ". Install them before running Study II.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

study2_taxonomy <- function(p = 100) {
  if (!identical(as.integer(p), 100L)) {
    stop("Study II taxonomy is hard-coded for p = 100.", call. = FALSE)
  }
  cbind(
    seq_len(p),
    rep(seq_len(10), each = 10),
    rep(seq_len(5), each = 20),
    c(rep(1L, 40), rep(2L, 40), rep(3L, 20))
  )
}

study2_tree_df <- function(p = 100) {
  if (!identical(as.integer(p), 100L)) {
    stop("Study II tree is hard-coded for p = 100.", call. = FALSE)
  }

  level10 <- p + seq_len(10)
  level20 <- p + 10L + seq_len(5)
  level40 <- p + 15L + seq_len(2)
  root <- p + 18L

  leaves <- data.frame(
    node = seq_len(p),
    parent = rep(level10, each = 10),
    name = paste0("leaf", seq_len(p)),
    weight = 0,
    stringsAsFactors = FALSE
  )
  nodes10 <- data.frame(
    node = level10,
    parent = rep(level20, each = 2),
    name = paste0("block10_", seq_len(10)),
    weight = 1,
    stringsAsFactors = FALSE
  )
  nodes20 <- data.frame(
    node = level20,
    parent = c(level40[1], level40[1], level40[2], level40[2], root),
    name = paste0("block20_", seq_len(5)),
    weight = 1,
    stringsAsFactors = FALSE
  )
  nodes40 <- data.frame(
    node = level40,
    parent = c(root, root),
    name = paste0("block40_", seq_len(2)),
    weight = 1,
    stringsAsFactors = FALSE
  )
  root_df <- data.frame(
    node = root,
    parent = NA_integer_,
    name = "root",
    weight = 1,
    stringsAsFactors = FALSE
  )

  out <- rbind(leaves, nodes10, nodes20, nodes40, root_df)
  out <- out[order(out$node), , drop = FALSE]
  rownames(out) <- NULL
  out
}

study2_threshold <- function(X, threshold = 0.005) {
  X <- as.matrix(X)
  if (any(X < 0)) stop("X must be non-negative.", call. = FALSE)
  row_sums <- rowSums(X)
  if (any(abs(row_sums - 1) > 1e-9)) {
    X <- X / row_sums
  }
  if (threshold > min(apply(X, 1, max))) {
    stop("Threshold is too large and would create at least one all-zero row.", call. = FALSE)
  }
  X_out <- X * (X >= threshold)
  X_out / rowSums(X_out)
}

study2_gamma_grid <- function() {
  seq(0.5, 10, by = 0.5) * 1e-2
}

study2_design <- function(
  seed = 20200630,
  n = 500,
  p = 100,
  threshold = 0.005,
  snr = 1
) {
  set.seed(seed)
  temp <- matrix(stats::rnorm(n * (p - 1L)), nrow = n, ncol = p - 1L)
  exp_temp <- exp(temp)
  X_true <- cbind(exp_temp, 1) / (rowSums(exp_temp) + 1)
  X_obs <- study2_threshold(X_true, threshold = threshold)

  beta <- c(
    rep(1, 20),
    rep(-2, 10),
    rep(0.5, 10),
    rep(2, 40),
    stats::rnorm(20)
  )
  sigma <- sqrt(stats::var(as.numeric(X_true %*% beta)) / snr)

  list(
    n = n,
    p = p,
    threshold = threshold,
    snr = snr,
    X_true = X_true,
    X_obs = X_obs,
    beta = beta,
    beta_group = c(20, 10, 10, 40, rep(1, 20)),
    sigma = sigma,
    taxonomy = study2_taxonomy(p),
    tree_df = study2_tree_df(p),
    zero_prop = mean(rowMeans(X_obs == 0))
  )
}

study2_split <- function(design, rep_id, response_seed_base = 20300630) {
  set.seed(response_seed_base + as.integer(rep_id))
  n <- design$n
  E <- stats::rnorm(n, mean = 0, sd = design$sigma)
  y <- as.numeric(design$X_true %*% design$beta + E)
  train_index <- seq_len(100)
  test_index <- 101:n

  list(
    rep_id = as.integer(rep_id),
    y = y,
    X_train = design$X_obs[train_index, , drop = FALSE],
    y_train = y[train_index],
    X_test = design$X_obs[test_index, , drop = FALSE],
    y_test = y[test_index],
    train_index = train_index,
    test_index = test_index
  )
}

population_var <- function(x) {
  mean((x - mean(x))^2)
}

study2_equisparsity <- function(beta) {
  beta <- as.numeric(beta)
  if (anyNA(beta) || length(beta) < 80L) return(NA_real_)
  20 * population_var(beta[1:20]) +
    10 * population_var(beta[21:30]) +
    10 * population_var(beta[31:40]) +
    40 * population_var(beta[41:80])
}

study2_eval_beta <- function(method, beta, split, selected_param = NA_real_,
                             valid_mse = NA_real_, intercept = 0,
                             elapsed_seconds = NA_real_, status = "ok") {
  beta <- as.numeric(beta)
  train_resid <- as.numeric(split$y_train - intercept - split$X_train %*% beta)
  test_resid <- as.numeric(split$y_test - intercept - split$X_test %*% beta)
  data.frame(
    rep_id = split$rep_id,
    method = method,
    selected_param = selected_param,
    valid_mse = valid_mse,
    train_sse = sum(train_resid^2),
    test_sse = sum(test_resid^2),
    test_mse = mean(test_resid^2),
    equi_sparsity = study2_equisparsity(beta),
    elapsed_seconds = elapsed_seconds,
    status = status,
    stringsAsFactors = FALSE
  )
}

study2_cv_folds <- function(n, nfold = 5, seed = 1) {
  set.seed(seed)
  sample(rep(seq_len(nfold), length.out = n))
}

study2_fit_rare_cv <- function(X_train, y_train, X_test, y_test, tree_df,
                               gamma_range, nfold = 5, seed = 1) {
  if (!requireNamespace("rare", quietly = TRUE)) {
    stop("Package 'rare' is required for RARE.", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for sparse RARE expansion.", call. = FALSE)
  }

  start <- proc.time()[["elapsed"]]
  n_train <- nrow(X_train)
  lambda <- as.numeric(gamma_range) / n_train
  A <- methods::as(
    Matrix::Matrix(rs_tree_df_to_expansion_matrix(tree_df, p = ncol(X_train)), sparse = TRUE),
    "dgCMatrix"
  )
  foldid <- study2_cv_folds(n_train, nfold = nfold, seed = seed)
  cv_mse_by_lambda <- numeric(length(lambda))

  for (fold in seq_len(nfold)) {
    valid <- foldid == fold
    fit <- rare::rarefit(
      y_train[!valid],
      X_train[!valid, , drop = FALSE],
      A = A,
      alpha = 1,
      intercept = FALSE,
      lambda = lambda
    )
    beta_path <- as.matrix(fit$beta[[1]])
    if (ncol(beta_path) != length(lambda)) {
      stop("RARE returned an unexpected number of lambda values.", call. = FALSE)
    }
    pred <- X_train[valid, , drop = FALSE] %*% beta_path
    cv_mse_by_lambda <- cv_mse_by_lambda +
      colMeans(sweep(pred, 1, y_train[valid], "-")^2) / nfold
  }

  best_idx <- which.min(cv_mse_by_lambda)
  final_fit <- rare::rarefit(
    y_train,
    X_train,
    A = A,
    alpha = 1,
    intercept = FALSE,
    lambda = lambda
  )
  beta_path <- as.matrix(final_fit$beta[[1]])
  beta <- beta_path[, best_idx]
  elapsed <- proc.time()[["elapsed"]] - start

  list(
    fit = final_fit,
    beta = beta,
    A = A,
    lambda = lambda,
    gamma = gamma_range,
    best_idx = best_idx,
    selected_gamma = gamma_range[best_idx],
    selected_lambda = lambda[best_idx],
    cv_mse = cv_mse_by_lambda[best_idx],
    cv_mse_by_lambda = cv_mse_by_lambda,
    test_mse = mean((y_test - X_test %*% beta)^2),
    elapsed_seconds = elapsed
  )
}

clean_rs_method_name <- function(x) {
  sub("^RS-L1$", "RS-L", x)
}
