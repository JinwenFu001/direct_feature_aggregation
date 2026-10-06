#!/usr/bin/env Rscript

cat("\n=== Experiment6: misspecified-tree paired simulation, O-LS + O-Ridge + treeFA + RARE + RS ===\n\n")
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

source_misspec_dgp <- function(experiment_dir, sim_dir, code_dir) {
  candidates <- c(
    Sys.getenv("MISSPEC_DGP_FILE", unset = ""),
    file.path(code_dir, "misspecified_tree_dgp.R"),
    file.path(experiment_dir, "misspecified_tree_dgp.R"),
    file.path(sim_dir, "misspecified_tree_dgp.R")
  )
  candidates <- candidates[nzchar(candidates)]
  dgp_file <- candidates[file.exists(candidates)][1]
  if (is.na(dgp_file)) {
    stop(
      "Cannot find misspecified_tree_dgp.R. Put it in TREEFA_CLUSTER_CODE_DIR or set MISSPEC_DGP_FILE.",
      call. = FALSE
    )
  }
  source(dgp_file)
  normalizePath(dgp_file, mustWork = TRUE)
}

coef_group_id <- function(beta) {
  as.integer(as.factor(as.numeric(beta)))
}

safe_ari <- function(x, y) {
  mclust::adjustedRandIndex(as.integer(x), as.integer(y))
}

noise_free_pe <- function(X, beta, true_beta, intercept = 0) {
  mean((as.numeric(X %*% true_beta) - intercept - as.numeric(X %*% beta))^2)
}

observed_mse <- function(X, y, beta, intercept = 0) {
  mean((as.numeric(y) - intercept - as.numeric(X %*% beta))^2)
}

least_squares_coef <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  decomp <- svd(X)
  if (!length(decomp$d)) return(numeric(ncol(X)))

  tol <- max(dim(X)) * max(decomp$d) * .Machine$double.eps
  keep <- decomp$d > tol
  if (!any(keep)) return(numeric(ncol(X)))

  as.numeric(
    decomp$v[, keep, drop = FALSE] %*%
      (crossprod(decomp$u[, keep, drop = FALSE], y) / decomp$d[keep])
  )
}

make_eval_row_exp6 <- function(
  level_label,
  level,
  target_c,
  gamma_hat,
  method,
  selected_param,
  beta,
  X_valid,
  Y_valid,
  X_test,
  Y_test,
  true_beta,
  true_group,
  oracle_group,
  beta_ideal = NULL,
  selected_param_ideal = NA_real_,
  intercept = 0,
  rand_beta = beta
) {
  beta <- as.numeric(beta)
  rand_beta <- as.numeric(rand_beta)
  beta_groups <- coef_group_id(rand_beta)
  true_beta_groups <- coef_group_id(true_beta)

  oracle_ari <- safe_ari(true_group, oracle_group)
  rand_true_group <- safe_ari(true_group, beta_groups)
  rand_true_beta <- safe_ari(true_beta_groups, beta_groups)

  out <- data.frame(
    variant = level_label,
    level = level,
    target_c = target_c,
    gamma_hat = gamma_hat,
    method = method,
    selected_param = selected_param,
    selected_param_ideal = selected_param_ideal,
    valid_mse = observed_mse(X_valid, Y_valid, beta, intercept = intercept),
    test_mse = noise_free_pe(X_test, beta, true_beta, intercept = intercept),
    test_y_mse = observed_mse(X_test, Y_test, beta, intercept = intercept),
    test_mse_ideal = NA_real_,
    rand = rand_true_beta,
    rand_true_group = rand_true_group,
    rand_true_beta = rand_true_beta,
    oracle_rand_true_group = oracle_ari,
    normalized_rand_true_group = if (is.finite(oracle_ari) && abs(oracle_ari) > .Machine$double.eps) {
      rand_true_group / oracle_ari
    } else {
      NA_real_
    },
    rand_ideal = NA_real_,
    rand_true_group_ideal = NA_real_,
    rand_true_beta_ideal = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(beta_ideal)) {
    beta_ideal <- as.numeric(beta_ideal)
    ideal_groups <- coef_group_id(beta_ideal)
    out$test_mse_ideal <- noise_free_pe(X_test, beta_ideal, true_beta, intercept = 0)
    out$rand_ideal <- safe_ari(true_beta_groups, ideal_groups)
    out$rand_true_group_ideal <- safe_ari(true_group, ideal_groups)
    out$rand_true_beta_ideal <- safe_ari(true_beta_groups, ideal_groups)
  }

  out
}

experiment_dir <- Sys.getenv("TREEFA_EXPERIMENT6_DIR", unset = script_dir())
experiment_dir <- normalizePath(experiment_dir, mustWork = TRUE)

sim_dir <- Sys.getenv("TREEFA_SIM_DIR", unset = experiment_dir)
sim_dir <- normalizePath(sim_dir, mustWork = TRUE)

# This script lives in code/run code/simulations/experiment6/.
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
dgp_file <- source_misspec_dgp(experiment_dir = experiment_dir, sim_dir = sim_dir, code_dir = code_dir)

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

uu_text <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = Sys.getenv("TREEFA_UU", unset = "1"))
uu <- as.numeric(uu_text)
if (!is.finite(uu) || uu < 1) stop("uu must be a positive number. Got: ", uu_text, call. = FALSE)

outer_reps <- as.integer(Sys.getenv("EXPERIMENT6_OUTER_REPS", unset = "10"))
inner_reps <- as.integer(Sys.getenv("EXPERIMENT6_INNER_REPS", unset = "40"))
nreps <- outer_reps * inner_reps

weight_grid <- -0.5
weight_index <- ((uu - 1) %/% nreps) + 1L
if (weight_index > length(weight_grid)) {
  stop("uu=", uu, " is outside the configured array range 1:", nreps * length(weight_grid), call. = FALSE)
}
weight.order <- weight_grid[weight_index]
uu_rep <- (uu - 1) %% nreps + 1L
ids <- misspec_outer_inner_id(rep_id = uu_rep, inner_reps = inner_reps)

n <- as.integer(Sys.getenv("EXPERIMENT6_N", unset = "50"))
n0 <- as.integer(Sys.getenv("EXPERIMENT6_N0", unset = "500"))
n1 <- as.integer(Sys.getenv("EXPERIMENT6_N1", unset = "50"))
p <- as.integer(Sys.getenv("EXPERIMENT6_P", unset = "200"))
q <- as.integer(Sys.getenv("EXPERIMENT6_Q", unset = "10"))
INratio <- as.numeric(Sys.getenv("EXPERIMENT6_INRATIO", unset = "5"))
poisson_rate <- as.numeric(Sys.getenv("EXPERIMENT6_POISSON_RATE", unset = "0.02"))
c_levels <- parse_numeric_vector(
  Sys.getenv("EXPERIMENT6_C_LEVELS", unset = "0,0.10,0.20,0.35,0.55,0.80,1.10"),
  c(0, 0.10, 0.20, 0.35, 0.55, 0.80, 1.10)
)
calibrate_on <- Sys.getenv("EXPERIMENT6_CALIBRATE_ON", unset = "train")
add_super_root <- identical(tolower(Sys.getenv("EXPERIMENT6_ADD_SUPER_ROOT", unset = "true")), "true")
base_seed <- as.integer(Sys.getenv("EXPERIMENT6_BASE_SEED", unset = "20260706"))

thresh <- as.numeric(Sys.getenv("TREEFA_THRESH", unset = "1e-5"))
max_iter <- as.integer(Sys.getenv("TREEFA_MAX_ITER", unset = "1000000"))
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
  Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD_METHODS", unset = "RS-DL2,RS-CL2"),
  c("RS-DL2", "RS-CL2")
)
rs_rand_threshold <- as.numeric(Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD", unset = "1e-4"))
rs_rand_threshold_norm <- Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD_NORM", unset = "scaled")
if (!rs_rand_threshold_norm %in% c("scaled", "unscaled")) {
  stop("TREE_VARIANT_RS_RAND_THRESHOLD_NORM must be 'scaled' or 'unscaled'.", call. = FALSE)
}

experiment_name <- Sys.getenv("TREEFA_EXPERIMENT_NAME", unset = "experiment6")
out_dir <- Sys.getenv("TREEFA_RESULT_DIR", unset = file.path(project_root, "output", experiment_name))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
experiment6_result_version <- "experiment6_with_oracle_v1"
includes_oracle_estimators <- TRUE

rs_work_root <- Sys.getenv("TREEFA_RS_WORKDIR", unset = file.path(out_dir, "rs_runs", paste0("array_", uu)))
dir.create(rs_work_root, recursive = TRUE, showWarnings = FALSE)

cat("Experiment dir:", experiment_dir, "\n")
cat("Simulation dir:", sim_dir, "\n")
cat("Code dir:", code_dir, "\n")
cat("DGP file:", dgp_file, "\n")
cat("R lib:", r_lib, "\n")
cat("Result dir:", out_dir, "\n")
cat("RS work dir:", rs_work_root, "\n")
cat("uu =", uu, "; uu_rep =", uu_rep, "; outer =", ids$outer, "; inner =", ids$inner, "\n")
cat("outer_reps =", outer_reps, "; inner_reps =", inner_reps, "; weight.order =", weight.order, "\n")
cat("n =", n, "; n0 =", n0, "; n1 =", n1, "; p =", p, "; q =", q, "\n")
cat("c_levels =", paste(c_levels, collapse = ","), "; INratio =", INratio, "; calibrate_on =", calibrate_on, "\n")
cat("treeFA stop_rule =", task_stop_rule, "; coarest_rule =", task_coarest_rule, "; thresh =", thresh, "; max_iter =", max_iter, "\n")
cat("RS methods:", paste(rs_methods, collapse = ", "), "; gamma nlam =", rs_gamma_nlam, "; min ratio =", rs_gamma_min_ratio, "\n")
cat("RS rand threshold:", paste(rs_rand_threshold_methods, collapse = ", "), "; threshold =", rs_rand_threshold, "; norm =", rs_rand_threshold_norm, "\n\n")
flush.console()

dat <- simulate_data_misspecified_tree(
  outer_id = ids$outer,
  inner_id = ids$inner,
  c_levels = c_levels,
  n = n,
  n0 = n0,
  n1 = n1,
  p = p,
  q = q,
  ratio = INratio,
  poisson_rate = poisson_rate,
  weight.order = weight.order,
  calibrate_on = calibrate_on,
  add_super_root = add_super_root,
  base_seed = base_seed
)

gamma_error <- max(abs(dat$gamma_info$gamma_hat - dat$gamma_info$target_c))
if (!is.finite(gamma_error) || gamma_error > 1e-6) {
  stop("DGP gamma verification failed. Max error = ", signif(gamma_error, 4), call. = FALSE)
}

X_train <- dat$X_train
Y_train <- dat$Y_train
X_valid <- dat$X_valid
Y_valid <- dat$Y_valid
X_test <- dat$X_test
Y_test <- dat$Y_test
signal_valid <- as.numeric(X_valid %*% dat$true.beta)
signal_test <- as.numeric(X_test %*% dat$true.beta)

## Fit the two oracle estimators once per replication. They use the true
## grouping and therefore do not depend on the misspecification level c.
true_group <- as.integer(dat$true_group)
true_group_levels <- sort(unique(true_group))
true_group_expansion <- Matrix::sparseMatrix(
  i = seq_len(p),
  j = match(true_group, true_group_levels),
  x = 1,
  dims = c(p, length(true_group_levels)),
  dimnames = list(NULL, paste0("G", true_group_levels))
)

oracle_X_train <- as.matrix(X_train %*% true_group_expansion)
oracle_X_valid <- as.matrix(X_valid %*% true_group_expansion)
oracle_X_test <- as.matrix(X_test %*% true_group_expansion)

log_msg("Running O-LS once for this replication...")
oracle_ls_group_coef <- least_squares_coef(oracle_X_train, Y_train)
oracle_ls_beta <- as.numeric(true_group_expansion %*% oracle_ls_group_coef)

log_msg("Running O-Ridge once for this replication...")
oracle_ridge_fit <- glmnet::glmnet(
  x = oracle_X_train,
  y = as.numeric(Y_train),
  alpha = 0,
  intercept = FALSE
)
oracle_ridge_lambda <- oracle_ridge_fit$lambda
oracle_ridge_valid_pred <- stats::predict(
  oracle_ridge_fit,
  newx = oracle_X_valid,
  s = oracle_ridge_lambda
)
oracle_ridge_valid_loss <- colMeans(
  (as.numeric(Y_valid) - oracle_ridge_valid_pred)^2
)
oracle_ridge_ideal_loss <- colMeans(
  (signal_valid - oracle_ridge_valid_pred)^2
)
oracle_ridge_best <- which.min(oracle_ridge_valid_loss)
oracle_ridge_ideal_best <- which.min(oracle_ridge_ideal_loss)

oracle_ridge_group_coef <- as.numeric(
  stats::coef(
    oracle_ridge_fit,
    s = oracle_ridge_lambda[oracle_ridge_best]
  )[-1, , drop = FALSE]
)
oracle_ridge_group_coef_ideal <- as.numeric(
  stats::coef(
    oracle_ridge_fit,
    s = oracle_ridge_lambda[oracle_ridge_ideal_best]
  )[-1, , drop = FALSE]
)
oracle_ridge_beta <- as.numeric(
  true_group_expansion %*% oracle_ridge_group_coef
)
oracle_ridge_beta_ideal <- as.numeric(
  true_group_expansion %*% oracle_ridge_group_coef_ideal
)

oracle_betas <- cbind(
  `O-LS` = oracle_ls_beta,
  `O-Ridge` = oracle_ridge_beta
)
oracle_betas_ideal <- cbind(
  `O-LS` = oracle_ls_beta,
  `O-Ridge` = oracle_ridge_beta_ideal
)
oracle_group_coefficients <- cbind(
  `O-LS` = oracle_ls_group_coef,
  `O-Ridge` = oracle_ridge_group_coef
)
oracle_group_coefficients_ideal <- cbind(
  `O-LS` = oracle_ls_group_coef,
  `O-Ridge` = oracle_ridge_group_coef_ideal
)
oracle_selected_params <- c(
  `O-LS` = NA_real_,
  `O-Ridge` = oracle_ridge_lambda[oracle_ridge_best]
)
oracle_selected_params_ideal <- c(
  `O-LS` = NA_real_,
  `O-Ridge` = oracle_ridge_lambda[oracle_ridge_ideal_best]
)

rs_gamma_range <- rs_rare_scaled_gamma_range(
  X_train,
  Y_train,
  intercept = FALSE,
  nlam = rs_gamma_nlam,
  lam.min.ratio = rs_gamma_min_ratio
)

cat("gamma_info:\n")
print(dat$gamma_info)
cat("RS/RARE gamma range:", signif(min(rs_gamma_range), 4), "to", signif(max(rs_gamma_range), 4), "\n\n")
flush.console()

combined_rows <- list()
oracle_rows <- list()
treefa_rows <- list()
rare_rows <- list()
rs_rows <- list()
selected_betas <- list()
rs_fits <- list()
rs_threshold_diagnostics <- list()
tree_list <- dat$tree_df_list
collapsed_tree_list <- list()
oracle_groups <- list()

for (level_idx in seq_along(c_levels)) {
  start_time <- proc.time()[["elapsed"]]
  level <- level_idx - 1L
  level_label <- paste0("c_", formatC(c_levels[level_idx], format = "f", digits = 2))
  tree_df <- dat$tree_df_list[[level_idx]]
  collapsed_tree_df <- collapse_zero_weight_tree(tree_df)
  collapsed_tree_list[[level_label]] <- collapsed_tree_df
  oracle_group <- partition_to_group_index(dat$partitions[[level_idx]], p = p)
  oracle_groups[[level_label]] <- oracle_group

  cat("\n--- level =", level, "; target c =", c_levels[level_idx], "(", level_idx, "/", length(c_levels), ") ---\n")
  flush.console()

  oracle_rows[[level_label]] <- rbind(
    make_eval_row_exp6(
      level_label = level_label,
      level = level,
      target_c = c_levels[level_idx],
      gamma_hat = dat$gamma_info$gamma_hat[level_idx],
      method = "O-LS",
      selected_param = NA_real_,
      beta = oracle_ls_beta,
      X_valid = X_valid,
      Y_valid = Y_valid,
      X_test = X_test,
      Y_test = Y_test,
      true_beta = dat$true.beta,
      true_group = dat$true_group,
      oracle_group = oracle_group,
      beta_ideal = oracle_ls_beta
    ),
    make_eval_row_exp6(
      level_label = level_label,
      level = level,
      target_c = c_levels[level_idx],
      gamma_hat = dat$gamma_info$gamma_hat[level_idx],
      method = "O-Ridge",
      selected_param = oracle_selected_params[["O-Ridge"]],
      selected_param_ideal = oracle_selected_params_ideal[["O-Ridge"]],
      beta = oracle_ridge_beta,
      X_valid = X_valid,
      Y_valid = Y_valid,
      X_test = X_test,
      Y_test = Y_test,
      true_beta = dat$true.beta,
      true_group = dat$true_group,
      oracle_group = oracle_group,
      beta_ideal = oracle_ridge_beta_ideal
    )
  )
  ## The normalized tree-oracle ARI is not meaningful for estimators that
  ## are given the true grouping and are not restricted to the input tree.
  oracle_rows[[level_label]]$normalized_rand_true_group <- NA_real_

  log_msg("Running treeFA for ", level_label, "...")
  treefa_fit <- treeFA::grid.simple_linear(
    Y = Y_train,
    X = X_train,
    tree_df = tree_df,
    true_beta = dat$true.beta,
    ridge.param = 0,
    thresh = thresh,
    max_iter = max_iter,
    stop_rule = task_stop_rule,
    coarest_rule = task_coarest_rule
  )
  treefa_valid_loss <- colMeans((as.numeric(Y_valid) - X_valid %*% treefa_fit$beta)^2)
  treefa_ideal_loss <- colMeans((signal_valid - X_valid %*% treefa_fit$beta)^2)
  treefa_best <- which.min(treefa_valid_loss)
  treefa_ideal_best <- which.min(treefa_ideal_loss)
  treefa_beta <- treefa_fit$beta[, treefa_best]
  treefa_beta_ideal <- treefa_fit$beta[, treefa_ideal_best]
  treefa_rows[[level_label]] <- make_eval_row_exp6(
    level_label = level_label,
    level = level,
    target_c = c_levels[level_idx],
    gamma_hat = dat$gamma_info$gamma_hat[level_idx],
    method = "treeFA",
    selected_param = treefa_fit$lambda[treefa_best],
    selected_param_ideal = treefa_fit$lambda[treefa_ideal_best],
    beta = treefa_beta,
    X_valid = X_valid,
    Y_valid = Y_valid,
    X_test = X_test,
    Y_test = Y_test,
    true_beta = dat$true.beta,
    true_group = dat$true_group,
    oracle_group = oracle_group,
    beta_ideal = treefa_beta_ideal
  )

  log_msg("Running RARE for ", level_label, "...")
  A_variant <- make_sparse_expansion_from_tree_df(collapsed_tree_df, p = p)
  rare_fit <- rare::rarefit(
    Y_train,
    X_train,
    A = A_variant,
    alpha = 1,
    intercept = FALSE,
    lambda = rs_gamma_range / nrow(X_train)
  )
  rare_valid_loss <- colMeans((as.numeric(Y_valid) - X_valid %*% rare_fit$beta[[1]])^2)
  rare_ideal_loss <- colMeans((signal_valid - X_valid %*% rare_fit$beta[[1]])^2)
  rare_best <- which.min(rare_valid_loss)
  rare_ideal_best <- which.min(rare_ideal_loss)
  rare_beta <- rare_fit$beta[[1]][, rare_best]
  rare_beta_ideal <- rare_fit$beta[[1]][, rare_ideal_best]
  rare_rows[[level_label]] <- make_eval_row_exp6(
    level_label = level_label,
    level = level,
    target_c = c_levels[level_idx],
    gamma_hat = dat$gamma_info$gamma_hat[level_idx],
    method = "RARE",
    selected_param = rare_fit$lambda[rare_best] * nrow(X_train),
    selected_param_ideal = rare_fit$lambda[rare_ideal_best] * nrow(X_train),
    beta = rare_beta,
    X_valid = X_valid,
    Y_valid = Y_valid,
    X_test = X_test,
    Y_test = Y_test,
    true_beta = dat$true.beta,
    true_group = dat$true_group,
    oracle_group = oracle_group,
    beta_ideal = rare_beta_ideal
  )

  log_msg("Running MATLAB RS for ", level_label, "...")
  rs_fit <- run_rs_matlab_once(
    X_train = X_train,
    y_train = Y_train,
    X_valid = X_valid,
    y_valid = Y_valid,
    X_test = X_test,
    y_test = Y_test,
    tree_df = collapsed_tree_df,
    work_dir = file.path(rs_work_root, level_label),
    seed = uu,
    gamma_range = rs_gamma_range,
    selection = rs_selection,
    model = rs_model,
    normalize_rows = rs_normalize_rows,
    mu = rs_mu,
    methods = rs_methods,
    keep_files = rs_keep_files
  )
  rs_fits[[level_label]] <- rs_fit

  present_rs_methods <- intersect(rs_methods, rs_fit$summary$method)
  if (!length(present_rs_methods)) {
    stop("None of EXPERIMENT6_RS_METHODS were returned by MATLAB RS.", call. = FALSE)
  }

  rs_rows[[level_label]] <- do.call(rbind, lapply(present_rs_methods, function(method) {
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
      rs_threshold_diagnostics[[paste(level_label, method, sep = "_")]] <<- thresholded$diagnostics
      log_msg(
        "RS Rand threshold for ", level_label, " ", method,
        ": top-down groups = ", thresholded$diagnostics$n_topdown_aggregated_groups,
        ", penalty groups below threshold = ", thresholded$diagnostics$n_thresholded_penalty_groups
      )
    }

    make_eval_row_exp6(
      level_label = level_label,
      level = level,
      target_c = c_levels[level_idx],
      gamma_hat = dat$gamma_info$gamma_hat[level_idx],
      method = method,
      selected_param = rs_fit$summary$gamma[j],
      beta = beta,
      X_valid = X_valid,
      Y_valid = Y_valid,
      X_test = X_test,
      Y_test = Y_test,
      true_beta = dat$true.beta,
      true_group = dat$true_group,
      oracle_group = oracle_group,
      intercept = intercept,
      rand_beta = rand_beta
    )
  }))

  selected_betas[[level_label]] <- cbind(
    oracle_betas,
    treeFA = treefa_beta,
    RARE = rare_beta,
    rs_fit$beta[, present_rs_methods, drop = FALSE]
  )

  combined_rows[[level_label]] <- rbind(
    oracle_rows[[level_label]],
    treefa_rows[[level_label]],
    rare_rows[[level_label]],
    rs_rows[[level_label]]
  )
  log_msg("Finished ", level_label, ". elapsed seconds: ", round(proc.time()[["elapsed"]] - start_time, 2))
  print(combined_rows[[level_label]])
  flush.console()
}

oracle_result <- do.call(rbind, oracle_rows)
treefa_result <- do.call(rbind, treefa_rows)
rare_result <- do.call(rbind, rare_rows)
rs_result <- do.call(rbind, rs_rows)
combined_result <- do.call(rbind, combined_rows)
result <- combined_result

test_mse_table <- stats::reshape(
  combined_result[, c("variant", "method", "test_mse")],
  idvar = "variant",
  timevar = "method",
  direction = "wide"
)
test_y_mse_table <- stats::reshape(
  combined_result[, c("variant", "method", "test_y_mse")],
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
rand_true_group_table <- stats::reshape(
  combined_result[, c("variant", "method", "rand_true_group")],
  idvar = "variant",
  timevar = "method",
  direction = "wide"
)

out_file <- file.path(
  out_dir,
  paste0(
    "MisspecifiedTree_n50_p200_q10_outer10_inner40_weight",
    weight.order,
    "_uu",
    uu_rep,
    "_treeFA_RARE_RS.RData"
  )
)

save(
  result,
  experiment6_result_version,
  includes_oracle_estimators,
  combined_result,
  oracle_result,
  treefa_result,
  rare_result,
  rs_result,
  selected_betas,
  oracle_betas,
  oracle_betas_ideal,
  oracle_group_coefficients,
  oracle_group_coefficients_ideal,
  oracle_selected_params,
  oracle_selected_params_ideal,
  true_group_expansion,
  tree_list,
  collapsed_tree_list,
  oracle_groups,
  rs_fits,
  rs_threshold_diagnostics,
  rs_gamma_range,
  rs_methods,
  dat,
  test_mse_table,
  test_y_mse_table,
  rand_table,
  rand_true_group_table,
  uu,
  uu_rep,
  ids,
  outer_reps,
  inner_reps,
  n,
  n0,
  n1,
  p,
  q,
  c_levels,
  INratio,
  weight.order,
  file = out_file
)

cat("\n=== Experiment6 finished ===\n")
cat("Saved result to:\n  ", out_file, "\n\n", sep = "")
cat("Combined result:\n")
print(combined_result)
