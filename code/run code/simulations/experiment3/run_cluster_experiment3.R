#!/usr/bin/env Rscript

cat("\n=== Experiment3: Bien tree k-increasing simulation with treeFA + RARE + RS + oLS + oRidge ===\n\n")
flush.console()

log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = ""), "\n", sep = "")
  flush.console()
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    file_path <- sub("^--file=", "", file_arg[1])
    # Rscript may encode spaces as ~+~ in commandArgs(), even for quoted paths.
    if (!file.exists(file_path)) file_path <- gsub("~+~", " ", file_path, fixed = TRUE)
    return(dirname(normalizePath(file_path, mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

parse_numeric_vector <- function(x, default) {
  if (!nzchar(x)) return(default)
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

parse_character_vector <- function(x, default) {
  if (!nzchar(x)) return(default)
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

experiment_dir <- Sys.getenv("TREEFA_EXPERIMENT3_DIR", unset = script_dir())
experiment_dir <- normalizePath(experiment_dir, mustWork = TRUE)

sim_dir <- Sys.getenv("TREEFA_SIM_DIR", unset = experiment_dir)
sim_dir <- normalizePath(sim_dir, mustWork = TRUE)

# This script lives in code/run code/simulations/experiment3/.
project_root <- Sys.getenv(
  "TREEFA_PROJECT_ROOT",
  unset = file.path(script_dir(), "..", "..", "..", "..")
)
project_root <- normalizePath(project_root, mustWork = TRUE)
default_code_dir <- file.path(project_root, "code", "core_code")
code_dir <- Sys.getenv("TREEFA_CLUSTER_CODE_DIR", unset = default_code_dir)
code_dir <- normalizePath(code_dir, mustWork = TRUE)
Sys.setenv(TREEFA_CLUSTER_CODE_DIR = code_dir)

r_lib <- Sys.getenv("TREEFA_R_LIB", unset = file.path(code_dir, "rlib"))
.libPaths(c(r_lib, .libPaths()))

source(file.path(code_dir, "rs_matlab_functions_cluster.R"))
source(file.path(code_dir, "other_functions_treeFA.R"))
source(file.path(code_dir, "variant_tree_helpers.R"))
source(file.path(code_dir, "adaptive_rs_helpers.R"))
source(file.path(code_dir, "rs_threshold_rand_helpers.R"))

required <- c("treeFA", "rare", "mclust", "glmnet", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing required R package(s): ", paste(missing, collapse = ", "),
    ". Install treeFA and dependencies in TREEFA_R_LIB before submitting.",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(treeFA))

simulate_data_experiment3 <- function(n, group.index, s = 0, beta.pre = NULL, ratio = 5) {
  k <- length(unique(group.index))
  p <- length(group.index)
  A <- matrix(0, nrow = p, ncol = k)

  for (i in seq_len(k)) {
    A[, i] <- group.index == i
  }

  if (is.null(beta.pre)) {
    beta0 <- stats::runif(k, 1.5, 2.5)
  } else {
    beta0 <- beta.pre
  }

  nonzero.ind <- sample(seq_len(k), ceiling(k * (1 - s)))
  nonzero <- rep(0, k)
  nonzero[nonzero.ind] <- rep(c(-1, 1), k %/% 2)[seq_along(nonzero.ind)]
  beta0 <- beta0 * nonzero

  beta <- A %*% beta0
  X <- matrix(stats::rpois(n * p, 0.02), nrow = n, ncol = p)
  means <- X %*% beta
  sigma <- sqrt(sum(means^2) / ratio / n)
  Y <- means + stats::rnorm(n, 0, sigma)

  list(X = X, Y = as.numeric(Y), true.beta = as.numeric(beta), A = A)
}

make_eval_row <- function(
  k_value,
  p,
  method,
  selected_param,
  beta,
  X_valid,
  Y_valid,
  X_test,
  Y_test,
  true_beta,
  beta_ideal = NULL,
  intercept = 0,
  rand_beta = beta
) {
  true_test_signal <- as.numeric(X_test %*% true_beta)

  out <- data.frame(
    variant = paste0("k_", k_value),
    k = k_value,
    p = p,
    method = method,
    selected_param = selected_param,
    valid_mse = mean((as.numeric(Y_valid) - intercept - X_valid %*% beta)^2),
    test_mse = mean((true_test_signal - intercept - X_test %*% beta)^2),
    test_mse_ideal = NA_real_,
    rand = mclust::adjustedRandIndex(
      as.numeric(as.factor(true_beta)),
      as.numeric(as.factor(rand_beta))
    ),
    rand_ideal = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(beta_ideal)) {
    out$test_mse_ideal <- mean((true_test_signal - X_test %*% beta_ideal)^2)
    out$rand_ideal <- mclust::adjustedRandIndex(
      as.numeric(as.factor(true_beta)),
      as.numeric(as.factor(beta_ideal))
    )
  }

  out
}

uu_text <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = Sys.getenv("TREEFA_UU", unset = "1"))
uu <- as.numeric(uu_text)
if (!is.finite(uu) || uu < 1) stop("uu must be a positive number. Got: ", uu_text, call. = FALSE)

nreps <- as.integer(Sys.getenv("EXPERIMENT3_NREPS", unset = "200"))
n <- as.integer(Sys.getenv("EXPERIMENT3_N", unset = "50"))
n0 <- as.integer(Sys.getenv("EXPERIMENT3_N0", unset = as.character(n * 10L)))
n1 <- as.integer(Sys.getenv("EXPERIMENT3_N1", unset = as.character(n)))
p <- as.integer(Sys.getenv("EXPERIMENT3_P", unset = "100"))
ks <- parse_numeric_vector(Sys.getenv("EXPERIMENT3_KS", unset = ""), seq(10, 50, 10))
ks <- as.integer(ks)
s <- as.numeric(Sys.getenv("EXPERIMENT3_S", unset = "0"))
INratio <- as.numeric(Sys.getenv("EXPERIMENT3_INRATIO", unset = "5"))
thresh <- as.numeric(Sys.getenv("TREEFA_THRESH", unset = "1e-5"))

weight_grid <- -1 / 2  # Only inverse-square-root leaf-count weights.
weight_index <- ((uu - 1) %/% nreps) + 1L
if (weight_index > length(weight_grid)) {
  stop("uu=", uu, " is outside the configured array range 1:", nreps * length(weight_grid), call. = FALSE)
}
weight.order <- weight_grid[weight_index]
uu_rep <- (uu - 1) %% nreps + 1L

task_stop_rule <- Sys.getenv("TREEFA_STOP_RULE", unset = "coef")
task_coarest_rule <- Sys.getenv("TREEFA_COAREST_RULE", unset = "revised")

rs_methods <- c("RS-DL2", "RS-CL2")
rs_selection <- "validation"
rs_model <- Sys.getenv("TREE_VARIANT_RS_MODEL", unset = "article_aligned")
rs_normalize_rows <- identical(tolower(Sys.getenv("TREE_VARIANT_RS_NORMALIZE_ROWS", unset = "false")), "true")
rs_mu_text <- Sys.getenv("TREE_VARIANT_RS_MU", unset = "")
rs_mu <- if (nzchar(rs_mu_text)) as.numeric(rs_mu_text) else NULL
rs_gamma_nlam <- as.integer(Sys.getenv("TREE_VARIANT_RS_GAMMA_NLAM", unset = "50"))
rs_gamma_min_ratio <- as.numeric(Sys.getenv("TREE_VARIANT_RS_GAMMA_MIN_RATIO", unset = "1e-4"))
rs_keep_files <- identical(tolower(Sys.getenv("TREEFA_KEEP_RS_FILES", unset = "false")), "true")

rs_rand_threshold_methods <- parse_character_vector(
  Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD_METHODS", unset = ""),
  c("RS-DL2", "RS-CL2")
)
rs_rand_threshold <- as.numeric(Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD", unset = "1e-4"))
rs_rand_threshold_norm <- Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD_NORM", unset = "scaled")
if (!rs_rand_threshold_norm %in% c("scaled", "unscaled")) {
  stop("TREE_VARIANT_RS_RAND_THRESHOLD_NORM must be 'scaled' or 'unscaled'.", call. = FALSE)
}

experiment_name <- Sys.getenv("TREEFA_EXPERIMENT_NAME", unset = "experiment3")
out_dir <- Sys.getenv("TREEFA_RESULT_DIR", unset = file.path(project_root, "output", experiment_name))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rs_work_root <- Sys.getenv("TREEFA_RS_WORKDIR", unset = file.path(out_dir, "rs_runs", paste0("array_", uu)))
dir.create(rs_work_root, recursive = TRUE, showWarnings = FALSE)

cat("Experiment dir:", experiment_dir, "\n")
cat("Simulation dir:", sim_dir, "\n")
cat("Code dir:", code_dir, "\n")
cat("R lib:", r_lib, "\n")
cat("Result dir:", out_dir, "\n")
cat("RS work dir:", rs_work_root, "\n")
cat("uu =", uu, "; replicate =", uu_rep, "; weight.order =", weight.order, "\n")
cat("n =", n, "; n0 =", n0, "; n1 =", n1, "; p =", p, "; ks =", paste(ks, collapse = ","), "\n")
cat("Experiment3 data: X ~ Pois(0.02), no spectral scaling, sigma from INratio =", INratio, "\n")
cat("Tuning uses observed Y_valid; test_mse uses the noiseless signal X_test beta_true.\n")
cat("treeFA stop_rule =", task_stop_rule, "; coarest_rule =", task_coarest_rule, "\n")
cat("RS methods:", paste(rs_methods, collapse = ", "), "\n")
cat("RS model =", rs_model, "; selection =", rs_selection, "; normalize_rows =", rs_normalize_rows, "\n")
cat("Fixed RS/RARE gamma grid: nlam =", rs_gamma_nlam, "; min ratio =", rs_gamma_min_ratio, "\n")
cat(
  "RS Rand threshold methods:", paste(rs_rand_threshold_methods, collapse = ", "),
  "; threshold =", rs_rand_threshold,
  "; norm =", rs_rand_threshold_norm, "\n\n"
)
flush.console()

set.seed(uu)

combined_rows <- list()
treefa_rows <- list()
rare_rows <- list()
rs_rows <- list()
oracle_rows <- list()
selected_betas <- list()
true_betas <- list()
tree_groups <- list()
tree_list <- list()
collapsed_tree_list <- list()
split_indices <- list()
rs_gamma_ranges <- list()

for (idx in seq_along(ks)) {
  k_value <- ks[idx]
  label <- paste0("k_", k_value)
  start_time <- proc.time()[["elapsed"]]

  cat("\n--- k =", k_value, "; p =", p, "(", idx, "/", length(ks), ") ---\n")
  flush.console()

  trees <- simulate_tree(k = k_value, p = p)
  tree_df <- hclust_to_df(trees$tree, weight.order = weight.order)
  collapsed_tree_df <- collapse_zero_weight_tree(tree_df)

  data <- simulate_data_experiment3(n + n0 + n1, trees$group, s = s, ratio = INratio)
  train_index <- sample(seq_len(n + n0 + n1), n)
  valid_index <- sample(setdiff(seq_len(n + n0 + n1), train_index), n1)
  test_index <- setdiff(seq_len(n + n0 + n1), c(train_index, valid_index))

  X_train <- data$X[train_index, , drop = FALSE]
  Y_train <- data$Y[train_index]
  X_valid <- data$X[valid_index, , drop = FALSE]
  Y_valid <- data$Y[valid_index]
  X_test <- data$X[test_index, , drop = FALSE]
  Y_test <- data$Y[test_index]

  rs_gamma_range <- rs_rare_scaled_gamma_range(
    X_train,
    Y_train,
    intercept = FALSE,
    nlam = rs_gamma_nlam,
    lam.min.ratio = rs_gamma_min_ratio
  )
  rs_gamma_ranges[[label]] <- rs_gamma_range

  log_msg("Running treeFA for ", label, "...")
  treefa_fit <- treeFA::grid.simple_linear(
    Y = Y_train,
    X = X_train,
    tree_df = tree_df,
    true_beta = data$true.beta,
    ridge.param = 0,
    thresh = thresh,
    stop_rule = task_stop_rule,
    coarest_rule = task_coarest_rule
  )
  treefa_valid_loss <- colSums((as.vector(Y_valid) - X_valid %*% treefa_fit$beta)^2)
  treefa_ideal_loss <- colSums((as.vector(X_valid %*% data$true.beta) - X_valid %*% treefa_fit$beta)^2)
  treefa_best <- which.min(treefa_valid_loss)
  treefa_ideal_best <- which.min(treefa_ideal_loss)
  treefa_beta <- treefa_fit$beta[, treefa_best]
  treefa_beta_ideal <- treefa_fit$beta[, treefa_ideal_best]
  treefa_rows[[label]] <- make_eval_row(
    k_value = k_value,
    p = p,
    method = "treeFA",
    selected_param = treefa_fit$lambda[treefa_best],
    beta = treefa_beta,
    X_valid = X_valid,
    Y_valid = Y_valid,
    X_test = X_test,
    Y_test = Y_test,
    true_beta = data$true.beta,
    beta_ideal = treefa_beta_ideal
  )

  log_msg("Running RARE for ", label, "...")
  A_variant <- make_sparse_expansion_from_tree_df(collapsed_tree_df, p = p)
  rare_fit <- rare::rarefit(
    Y_train,
    X_train,
    A = A_variant,
    alpha = 1,
    intercept = FALSE,
    lambda = rs_gamma_range / nrow(X_train)
  )
  rare_valid_loss <- colSums((as.vector(Y_valid) - X_valid %*% rare_fit$beta[[1]])^2)
  rare_ideal_loss <- colSums((as.vector(X_valid %*% data$true.beta) - X_valid %*% rare_fit$beta[[1]])^2)
  rare_best <- which.min(rare_valid_loss)
  rare_ideal_best <- which.min(rare_ideal_loss)
  rare_beta <- rare_fit$beta[[1]][, rare_best]
  rare_beta_ideal <- rare_fit$beta[[1]][, rare_ideal_best]
  rare_rows[[label]] <- make_eval_row(
    k_value = k_value,
    p = p,
    method = "RARE",
    selected_param = rare_fit$lambda[rare_best] * nrow(X_train),
    beta = rare_beta,
    X_valid = X_valid,
    Y_valid = Y_valid,
    X_test = X_test,
    Y_test = Y_test,
    true_beta = data$true.beta,
    beta_ideal = rare_beta_ideal
  )

  log_msg("Running MATLAB RS for ", label, "...")
  rs_fit <- run_rs_matlab_once(
    X_train = X_train,
    y_train = Y_train,
    X_valid = X_valid,
    y_valid = Y_valid,
    X_test = X_test,
    y_test = Y_test,
    tree_df = collapsed_tree_df,
    work_dir = file.path(rs_work_root, label),
    seed = uu,
    gamma_range = rs_gamma_range,
    selection = rs_selection,
    model = rs_model,
    normalize_rows = rs_normalize_rows,
    mu = rs_mu,
    keep_files = rs_keep_files,
    methods = rs_methods
  )

  present_rs_methods <- intersect(rs_methods, rs_fit$summary$method)
  if (!length(present_rs_methods)) {
    stop("None of EXPERIMENT3_RS_METHODS were returned by MATLAB RS.", call. = FALSE)
  }

  rs_rows[[label]] <- do.call(rbind, lapply(present_rs_methods, function(method) {
    j <- match(method, rs_fit$summary$method)
    beta <- rs_fit$beta[, method]
    intercept <- rs_fit$summary$intercept[j]
    rand_beta <- beta

    if (method %in% rs_rand_threshold_methods) {
      thresholded <- rs_thresholded_beta_for_rand(
        beta = beta,
        tree_df = collapsed_tree_df,
        X_train = X_train,
        y_train = Y_train,
        method = method,
        gamma = rs_fit$summary$gamma[j],
        model = rs_model,
        normalize_rows = rs_normalize_rows,
        mu = rs_mu,
        threshold = rs_rand_threshold,
        norm_mode = rs_rand_threshold_norm
      )
      rand_beta <- thresholded$beta
      log_msg(
        "RS Rand threshold for ", label, " ", method,
        ": top-down groups = ", thresholded$diagnostics$n_topdown_aggregated_groups,
        ", penalty groups below threshold = ", thresholded$diagnostics$n_thresholded_penalty_groups
      )
    }

    make_eval_row(
      k_value = k_value,
      p = p,
      method = method,
      selected_param = rs_fit$summary$gamma[j],
      beta = beta,
      X_valid = X_valid,
      Y_valid = Y_valid,
      X_test = X_test,
      Y_test = Y_test,
      true_beta = data$true.beta,
      intercept = intercept,
      rand_beta = rand_beta
    )
  }))

  ## Oracle baselines use the true grouping, as in experiments 1 and 2.
  oracle_Q <- data$A
  oracle_X_train <- X_train %*% oracle_Q
  oracle_X_valid <- X_valid %*% oracle_Q

  log_msg("Running oLS for ", label, "...")
  ## Match my_ginv(crossprod(ols.X)), including the original 1e-8 tolerance.
  oracle_svd <- svd(crossprod(oracle_X_train))
  oracle_d <- oracle_svd$d
  oracle_d[oracle_d < 1e-8 * max(oracle_d)] <- 0
  oracle_d_inv <- 1 / oracle_d
  oracle_d_inv[!is.finite(oracle_d_inv)] <- 0
  oracle_ginv <- oracle_svd$v %*% diag(oracle_d_inv, length(oracle_d_inv)) %*% t(oracle_svd$u)
  oracle_ls_coef <- oracle_ginv %*% crossprod(oracle_X_train, Y_train)
  oracle_ls_beta <- as.numeric(oracle_Q %*% oracle_ls_coef)

  log_msg("Running oRidge for ", label, "...")
  ## Keep the original glmnet defaults: ridge, no intercept, and default standardization.
  oracle_ridge_fit <- glmnet::glmnet(oracle_X_train, Y_train, alpha = 0, intercept = FALSE)
  oracle_ridge_lambda <- oracle_ridge_fit$lambda
  oracle_ridge_valid_pred <- stats::predict(
    oracle_ridge_fit, newx = oracle_X_valid, s = oracle_ridge_lambda
  )
  oracle_ridge_valid_mse <- colMeans((Y_valid - oracle_ridge_valid_pred)^2)
  oracle_ridge_best <- which.min(oracle_ridge_valid_mse)
  oracle_ridge_selected_lambda <- unname(oracle_ridge_lambda[oracle_ridge_best])
  oracle_ridge_coef <- as.numeric(stats::coef(oracle_ridge_fit, s = oracle_ridge_selected_lambda))[-1]
  oracle_ridge_beta <- as.numeric(oracle_Q %*% oracle_ridge_coef)

  oracle_rows[[label]] <- rbind(
    make_eval_row(
      k_value = k_value,
      p = p,
      method = "oLS",
      selected_param = NA_real_,
      beta = oracle_ls_beta,
      X_valid = X_valid,
      Y_valid = Y_valid,
      X_test = X_test,
      Y_test = Y_test,
      true_beta = data$true.beta
    ),
    make_eval_row(
      k_value = k_value,
      p = p,
      method = "oRidge",
      selected_param = oracle_ridge_selected_lambda,
      beta = oracle_ridge_beta,
      X_valid = X_valid,
      Y_valid = Y_valid,
      X_test = X_test,
      Y_test = Y_test,
      true_beta = data$true.beta
    )
  )

  selected_betas[[label]] <- cbind(
    treeFA = treefa_beta,
    RARE = rare_beta,
    rs_fit$beta[, present_rs_methods, drop = FALSE],
    oLS = oracle_ls_beta,
    oRidge = oracle_ridge_beta
  )
  true_betas[[label]] <- data$true.beta
  tree_groups[[label]] <- trees$group
  tree_list[[label]] <- tree_df
  collapsed_tree_list[[label]] <- collapsed_tree_df
  split_indices[[label]] <- list(train = train_index, valid = valid_index, test = test_index)

  combined_rows[[label]] <- rbind(treefa_rows[[label]], rare_rows[[label]], rs_rows[[label]], oracle_rows[[label]])

  log_msg("Finished ", label, ". Elapsed seconds: ", round(proc.time()[["elapsed"]] - start_time, 2))
  print(combined_rows[[label]])
  flush.console()
}

treefa_result <- do.call(rbind, treefa_rows)
rare_result <- do.call(rbind, rare_rows)
rs_result <- do.call(rbind, rs_rows)
oracle_result <- do.call(rbind, oracle_rows)
combined_result <- do.call(rbind, combined_rows)
result <- combined_result

test_mse_table <- stats::reshape(
  combined_result[, c("variant", "method", "test_mse")],
  idvar = "variant",
  timevar = "method",
  direction = "wide"
)
rand_table <- stats::reshape(
  combined_result[, c("variant", "method", "rand")],
  idvar = "variant",
  timevar = "method",
  direction = "wide"
)

out_file <- file.path(
  out_dir,
  paste0(
    "Valid_BienTree_n50_p100_kIncre_weight",
    weight.order,
    "_uu",
    uu_rep,
    "_treeFA_RARE_RS.RData"
  )
)

save(
  result,
  combined_result,
  treefa_result,
  rare_result,
  rs_result,
  oracle_result,
  selected_betas,
  true_betas,
  tree_groups,
  tree_list,
  collapsed_tree_list,
  split_indices,
  rs_gamma_ranges,
  rs_methods,
  rs_gamma_nlam,
  rs_gamma_min_ratio,
  rs_keep_files,
  rs_rand_threshold_methods,
  rs_rand_threshold,
  rs_rand_threshold_norm,
  test_mse_table,
  rand_table,
  uu,
  uu_rep,
  ks,
  n,
  n0,
  n1,
  p,
  INratio,
  weight.order,
  file = out_file
)

cat("\n=== FINISHED ===\n")
cat("Saved result to:\n  ", out_file, "\n\n", sep = "")
cat("Combined result:\n")
print(combined_result)
