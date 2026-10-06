#!/usr/bin/env Rscript

cat("\n=== treeFA + RARE + RS-DL2/RS-CL2 on variant trees ===\n\n")
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

sim_dir <- Sys.getenv("TREEFA_SIM_DIR", unset = script_dir())
sim_dir <- normalizePath(sim_dir, mustWork = TRUE)

# This script lives in code/run code/simulations/experiment1/.
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
    ". Install treeFA with install_treeFA_cluster.R and install CRAN/GitHub dependencies before submitting the array.",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(treeFA))

uu_text <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = Sys.getenv("TREEFA_UU", unset = "1"))
uu <- as.numeric(uu_text)
if (!is.finite(uu) || uu < 1) stop("uu must be a positive number. Got: ", uu_text, call. = FALSE)

task_stop_rule <- Sys.getenv("TREEFA_STOP_RULE", unset = "coef")
task_coarest_rule <- Sys.getenv("TREEFA_COAREST_RULE", unset = "revised")
rs_selection <- Sys.getenv("TREE_VARIANT_RS_SELECTION", unset = "validation")
rs_model <- Sys.getenv("TREE_VARIANT_RS_MODEL", unset = "article_aligned")
rs_normalize_rows <- identical(tolower(Sys.getenv("TREE_VARIANT_RS_NORMALIZE_ROWS", unset = "false")), "true")
rs_mu_text <- Sys.getenv("TREE_VARIANT_RS_MU", unset = "")
rs_mu <- if (nzchar(rs_mu_text)) as.numeric(rs_mu_text) else NULL

rs_adaptive_methods_text <- Sys.getenv("TREE_VARIANT_RS_ADAPTIVE_METHODS", unset = "")
rs_adaptive_methods <- trimws(strsplit(rs_adaptive_methods_text, ",", fixed = TRUE)[[1]])
rs_adaptive_methods <- rs_adaptive_methods[nzchar(rs_adaptive_methods)]
rs_adaptive_edge <- as.integer(Sys.getenv("TREE_VARIANT_RS_ADAPTIVE_EDGE", unset = "2"))
rs_adaptive_expand_factor <- as.numeric(Sys.getenv("TREE_VARIANT_RS_ADAPTIVE_EXPAND_FACTOR", unset = "10"))
rs_adaptive_max_expansions <- as.integer(Sys.getenv("TREE_VARIANT_RS_ADAPTIVE_MAX_EXPANSIONS", unset = "3"))
rs_gamma_nlam <- as.integer(Sys.getenv("TREE_VARIANT_RS_GAMMA_NLAM", unset = "50"))
rs_gamma_min_ratio <- as.numeric(Sys.getenv("TREE_VARIANT_RS_GAMMA_MIN_RATIO", unset = "1e-4"))
rs_keep_files <- identical(tolower(Sys.getenv("TREEFA_KEEP_RS_FILES", unset = "false")), "true")
rs_rand_threshold_methods_text <- Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD_METHODS", unset = "RS-DL2,RS-CL2")
rs_rand_threshold_methods <- trimws(strsplit(rs_rand_threshold_methods_text, ",", fixed = TRUE)[[1]])
rs_rand_threshold_methods <- rs_rand_threshold_methods[nzchar(rs_rand_threshold_methods)]
rs_rand_threshold <- as.numeric(Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD", unset = "1e-4"))
rs_rand_threshold_norm <- Sys.getenv("TREE_VARIANT_RS_RAND_THRESHOLD_NORM", unset = "scaled")
if (!is.finite(rs_rand_threshold) || rs_rand_threshold < 0) {
  stop("TREE_VARIANT_RS_RAND_THRESHOLD must be finite and non-negative.", call. = FALSE)
}
if (!rs_rand_threshold_norm %in% c("scaled", "unscaled")) {
  stop("TREE_VARIANT_RS_RAND_THRESHOLD_NORM must be 'scaled' or 'unscaled'.", call. = FALSE)
}

set.seed(20)
n <- 50
n0 <- 500
n1 <- n
k <- 10
p <- 60
nreps <- 200
weight_index <- ((uu - 1) %/% nreps) + 1
weight_grid <- c(-1 / 2, -1, -2)
if (weight_index > length(weight_grid)) {
  stop("uu=", uu, " is outside the configured array range 1:", nreps * length(weight_grid), call. = FALSE)
}
weight.order <- weight_grid[weight_index]
uu_rep <- (uu - 1) %% nreps + 1

experiment_name <- Sys.getenv("TREEFA_EXPERIMENT_NAME", unset = "experiment1")
default_out_dir <- file.path(project_root, "output", experiment_name)
out_dir <- Sys.getenv("TREEFA_RESULT_DIR", unset = default_out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rs_work_root <- Sys.getenv("TREEFA_RS_WORKDIR", unset = file.path(out_dir, "rs_runs"))
dir.create(rs_work_root, recursive = TRUE, showWarnings = FALSE)

cat("Simulation dir:", sim_dir, "\n")
cat("Code dir:", code_dir, "\n")
cat("R lib:", r_lib, "\n")
cat("Result dir:", out_dir, "\n")
cat("RS work dir:", rs_work_root, "\n")
cat("uu =", uu, "; replicate =", uu_rep, "; weight.order =", weight.order, "\n")
cat("treeFA stop_rule =", task_stop_rule, "; coarest_rule =", task_coarest_rule, "\n")
cat("RS model =", rs_model, "; selection =", rs_selection, "; normalize_rows =", rs_normalize_rows, "\n")
cat("Adaptive RS methods:", if (length(rs_adaptive_methods)) paste(rs_adaptive_methods, collapse = ", ") else "none", "\n")
cat("Fixed RS/RARE gamma grid: nlam =", rs_gamma_nlam, "; min ratio =", rs_gamma_min_ratio, "\n\n")
cat(
  "RS Rand threshold methods:",
  if (length(rs_rand_threshold_methods)) paste(rs_rand_threshold_methods, collapse = ", ") else "none",
  "; threshold =", rs_rand_threshold,
  "; norm =", rs_rand_threshold_norm, "\n"
)
cat("Keep RS temporary files:", rs_keep_files, "\n\n")
flush.console()

trees <- simulate_tree(k = k, p = p)
tree_df0 <- hclust_to_df(trees$tree, weight.order)
merge.nodes <- unlist(find_latent_nodes_with_exact_leaf_groups(tree_df0, trees$group), use.names = FALSE)
L1 <- find_all_parents(tree_df0, merge.nodes)
L0 <- seq_len(nrow(tree_df0))[-c(seq_len(p), L1, merge.nodes)]

l0_seq <- c(10, 20, 30)
l1_seq <- c(3, 6)

tree1_l0 <- tree_df0
tree1_l0$weight[sample(L0, l0_seq[1], replace = FALSE)] <- 0
tree2_l0 <- tree_df0
tree2_l0$weight[sample(L0, l0_seq[2], replace = FALSE)] <- 0
tree3_l0 <- tree_df0
tree3_l0$weight[sample(L0, l0_seq[3], replace = FALSE)] <- 0
tree1_l1 <- tree_df0
tree1_l1$weight[sample(L1, l1_seq[1], replace = FALSE)] <- 0
tree2_l1 <- tree_df0
tree2_l1$weight[sample(L1, l1_seq[2], replace = FALSE)] <- 0

df_list <- list(tree_df0, tree1_l0, tree2_l0, tree3_l0, tree1_l1, tree2_l1)
row_labels <- c("tree_df0", "tree1_l0", "tree2_l0", "tree3_l0", "tree1_l1", "tree2_l1")

collapsed_df_list <- lapply(df_list, collapse_zero_weight_tree)
variant_tree_summary <- do.call(rbind, Map(function(label, original, collapsed) {
  out <- tree_variant_summary(original, collapsed)
  out$variant <- label
  out[, c("variant", setdiff(names(out), "variant")), drop = FALSE]
}, row_labels, df_list, collapsed_df_list))

set.seed(uu)
data <- simulate_data(n + n0 + n1, trees$group, s = 0, ratio = 5)
train_index <- sample(seq_len(n + n0 + n1), n)
valid_index <- sample(setdiff(seq_len(n + n0 + n1), train_index), n1)
test_index <- setdiff(seq_len(n + n0 + n1), c(train_index, valid_index))

X_train <- data$X[train_index, , drop = FALSE]
Y_train <- data$Y[train_index]
X_valid <- data$X[valid_index, , drop = FALSE]
Y_valid <- data$Y[valid_index]
X_test <- data$X[test_index, , drop = FALSE]
Y_test <- data$Y[test_index]
# Keep noisy validation responses for tuning; evaluate against the noiseless test mean.
mu_test <- as.vector(X_test %*% data$true.beta)

rs_gamma_range <- rs_rare_scaled_gamma_range(
  X_train,
  Y_train,
  intercept = FALSE,
  nlam = rs_gamma_nlam,
  lam.min.ratio = rs_gamma_min_ratio
)
cat(
  "RS/RARE gamma range: ",
  signif(min(rs_gamma_range), 4), " to ", signif(max(rs_gamma_range), 4),
  " (", length(rs_gamma_range), " values)\n\n",
  sep = ""
)
flush.console()

make_eval_row <- function(variant, method, selected_param, beta, beta_ideal = NULL, valid_mse = NA_real_) {
  out <- data.frame(
    variant = variant,
    method = method,
    selected_param = selected_param,
    valid_mse = valid_mse,
    test_mse = mean((mu_test - X_test %*% beta)^2),
    test_mse_ideal = NA_real_,
    rand = mclust::adjustedRandIndex(
      as.numeric(as.factor(data$true.beta)),
      as.numeric(as.factor(beta))
    ),
    rand_ideal = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(beta_ideal)) {
    out$test_mse_ideal <- mean((mu_test - X_test %*% beta_ideal)^2)
    out$rand_ideal <- mclust::adjustedRandIndex(
      as.numeric(as.factor(data$true.beta)),
      as.numeric(as.factor(beta_ideal))
    )
  }
  out
}

# Oracle fits use the true groups and the same training/validation split as all methods.
# This is the make_new_X construction used in the original simulation code.
oracle_groups <- unique(trees$group)
oracle_Q <- matrix(0, nrow = p, ncol = length(oracle_groups))
for (j in seq_along(oracle_groups)) {
  oracle_Q[trees$group == oracle_groups[j], j] <- 1
}
oracle_X <- data$X %*% oracle_Q
oracle_X_train <- oracle_X[train_index, , drop = FALSE]
oracle_X_valid <- oracle_X[valid_index, , drop = FALSE]

# Match the original my_ginv(crossprod(ols.X)), including its tolerance.
oracle_svd <- svd(crossprod(oracle_X_train))
oracle_d <- oracle_svd$d
oracle_d[oracle_d < 1e-8 * max(oracle_d)] <- 0
oracle_d_inv <- 1 / oracle_d
oracle_d_inv[!is.finite(oracle_d_inv)] <- 0
oracle_ginv <- oracle_svd$v %*% diag(oracle_d_inv, length(oracle_d_inv)) %*% t(oracle_svd$u)
oracle_ls_coef <- oracle_ginv %*% crossprod(oracle_X_train, Y_train)
oracle_ls_beta <- as.numeric(oracle_Q %*% oracle_ls_coef)
oracle_ls_valid_mse <- mean((Y_valid - oracle_X_valid %*% oracle_ls_coef)^2)

# Match the original glmnet defaults: alpha = 0, no intercept, and validation tuning.
oracle_ridge_fit <- glmnet::glmnet(oracle_X_train, Y_train, alpha = 0, intercept = FALSE)
oracle_ridge_lambda <- oracle_ridge_fit$lambda
oracle_ridge_valid_pred <- stats::predict(
  oracle_ridge_fit, newx = oracle_X_valid, s = oracle_ridge_lambda
)
oracle_ridge_valid_mse <- colMeans((Y_valid - oracle_ridge_valid_pred)^2)
oracle_ridge_best <- which.min(oracle_ridge_valid_mse)
oracle_ridge_selected_lambda <- oracle_ridge_lambda[oracle_ridge_best]
oracle_ridge_coef <- as.numeric(stats::coef(oracle_ridge_fit, s = oracle_ridge_selected_lambda))[-1]
oracle_ridge_beta <- as.numeric(oracle_Q %*% oracle_ridge_coef)

true_group_ids <- as.numeric(as.factor(data$true.beta))

treefa_rows <- list()
rare_rows <- list()
rs_rows <- list()
oracle_rows <- list()
selected_betas <- list()

log_msg("Entering variant loop with ", length(row_labels), " tree variants.")

for (i in seq_along(row_labels)) {
  label <- row_labels[i]
  variant_start <- proc.time()[["elapsed"]]
  cat("\n--- Variant:", label, "---\n")
  cat(
    "Collapsed nodes:",
    nrow(df_list[[i]]), "->", nrow(collapsed_df_list[[i]]),
    "; removed zero internal nodes =",
    length(attr(collapsed_df_list[[i]], "removed_zero_nodes")), "\n"
  )
  flush.console()

  log_msg("Running treeFA for ", label, "...")
  treefa_fit <- treeFA::grid.simple_linear(
    Y = Y_train,
    X = X_train,
    tree_df = df_list[[i]],
    true_beta = data$true.beta,
    ridge.param = 0,
    thresh = 1e-5,
    stop_rule = task_stop_rule,
    coarest_rule = task_coarest_rule
  )
  treefa_valid_loss <- colSums((as.vector(Y_valid) - X_valid %*% treefa_fit$beta)^2)
  treefa_ideal_loss <- colSums((as.vector(X_valid %*% data$true.beta) - X_valid %*% treefa_fit$beta)^2)
  treefa_beta <- treefa_fit$beta[, which.min(treefa_valid_loss)]
  treefa_beta_ideal <- treefa_fit$beta[, which.min(treefa_ideal_loss)]
  treefa_rows[[i]] <- make_eval_row(
    variant = label,
    method = "treeFA",
    selected_param = treefa_fit$lambda[which.min(treefa_valid_loss)],
    beta = treefa_beta,
    beta_ideal = treefa_beta_ideal,
    valid_mse = min(treefa_valid_loss) / length(Y_valid)
  )
  log_msg("Finished treeFA for ", label, ".")

  log_msg("Running RARE with collapsed modified tree for ", label, "...")
  A_variant <- make_sparse_expansion_from_tree_df(collapsed_df_list[[i]], p = p)
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
  rare_rows[[i]] <- make_eval_row(
    variant = label,
    method = "RARE",
    selected_param = rare_fit$lambda[rare_best] * nrow(X_train),
    beta = rare_beta,
    beta_ideal = rare_beta_ideal,
    valid_mse = min(rare_valid_loss) / length(Y_valid)
  )
  log_msg("Finished RARE for ", label, ".")

  log_msg("Running MATLAB RS-DL2/RS-CL2 with collapsed modified tree for ", label, "...")
  rs_work_dir <- file.path(rs_work_root, paste0("array", uu, "_", label))
  log_msg("RS work directory: ", rs_work_dir)
  if (length(rs_adaptive_methods)) {
    rs_fit <- run_rs_matlab_adaptive_once(
      X_train = X_train,
      y_train = Y_train,
      X_valid = X_valid,
      y_valid = Y_valid,
      X_test = X_test,
      y_test = mu_test,
      tree_df = collapsed_df_list[[i]],
      work_dir = rs_work_dir,
      seed = uu,
      gamma_range = rs_gamma_range,
      selection = rs_selection,
      model = rs_model,
      normalize_rows = rs_normalize_rows,
      mu = rs_mu,
      adaptive_methods = rs_adaptive_methods,
      edge = rs_adaptive_edge,
      expand_factor = rs_adaptive_expand_factor,
      max_expansions = rs_adaptive_max_expansions,
      grid_length = length(rs_gamma_range),
      verbose = TRUE
    )
  } else {
    rs_fit <- run_rs_matlab_once(
      X_train = X_train,
      y_train = Y_train,
      X_valid = X_valid,
      y_valid = Y_valid,
      X_test = X_test,
      y_test = mu_test,
      tree_df = collapsed_df_list[[i]],
      work_dir = rs_work_dir,
      seed = uu,
      gamma_range = rs_gamma_range,
      selection = rs_selection,
      model = rs_model,
      normalize_rows = rs_normalize_rows,
      mu = rs_mu,
      methods = c("RS-DL2", "RS-CL2"),
      keep_files = rs_keep_files
    )
  }

  rs_rows[[i]] <- do.call(rbind, lapply(seq_len(nrow(rs_fit$summary)), function(j) {
    method <- rs_fit$summary$method[j]
    beta <- rs_fit$beta[, method]
    intercept <- rs_fit$summary$intercept[j]
    rand_beta <- beta
    if (method %in% rs_rand_threshold_methods) {
      thresholded <- rs_thresholded_beta_for_rand(
        beta = beta,
        tree_df = collapsed_df_list[[i]],
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
    data.frame(
      variant = label,
      method = method,
      selected_param = rs_fit$summary$gamma[j],
      valid_mse = rs_fit$summary$valid_mse[j],
      test_mse = mean((mu_test - intercept - X_test %*% beta)^2),
      test_mse_ideal = NA_real_,
      rand = mclust::adjustedRandIndex(
        true_group_ids,
        as.numeric(as.factor(rand_beta))
      ),
      rand_ideal = NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  log_msg("Finished MATLAB RS for ", label, ".")

  oracle_rows[[i]] <- rbind(
    make_eval_row(
      variant = label,
      method = "oLS",
      selected_param = NA_real_,
      beta = oracle_ls_beta,
      valid_mse = oracle_ls_valid_mse
    ),
    make_eval_row(
      variant = label,
      method = "oRidge",
      selected_param = oracle_ridge_selected_lambda,
      beta = oracle_ridge_beta,
      valid_mse = oracle_ridge_valid_mse[oracle_ridge_best]
    )
  )

  selected_betas[[label]] <- cbind(
    treeFA = treefa_beta,
    RARE = rare_beta,
    `RS-DL2` = rs_fit$beta[, "RS-DL2"],
    `RS-CL2` = rs_fit$beta[, "RS-CL2"],
    oLS = oracle_ls_beta,
    oRidge = oracle_ridge_beta
  )

  variant_result <- rbind(treefa_rows[[i]], rare_rows[[i]], rs_rows[[i]], oracle_rows[[i]])
  log_msg("Finished all methods for ", label, ". Elapsed seconds: ", round(proc.time()[["elapsed"]] - variant_start, 2))
  print(variant_result)
  flush.console()
}

treefa_result <- do.call(rbind, treefa_rows)
rare_result <- do.call(rbind, rare_rows)
rs_result <- do.call(rbind, rs_rows)
oracle_result <- do.call(rbind, oracle_rows)
combined_result <- rbind(treefa_result, rare_result, rs_result, oracle_result)
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
    "Valid_BienTree_n50_p60_k10_treeChange_weight",
    weight.order,
    "_uu",
    uu_rep,
    "_treeFA_RARE_RS_modifiedTree.RData"
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
  test_mse_table,
  rand_table,
  rs_gamma_range,
  rs_gamma_nlam,
  rs_gamma_min_ratio,
  rs_adaptive_methods,
  rs_adaptive_edge,
  rs_adaptive_expand_factor,
  rs_adaptive_max_expansions,
  rs_keep_files,
  rs_rand_threshold_methods,
  rs_rand_threshold,
  rs_rand_threshold_norm,
  variant_tree_summary,
  df_list,
  collapsed_df_list,
  data,
  train_index,
  valid_index,
  test_index,
  file = out_file
)

cat("\n=== FINISHED ===\n")
cat("Saved result to:\n  ", out_file, "\n\n", sep = "")
cat("Variant tree summary:\n")
print(variant_tree_summary)
cat("\nCombined result:\n")
print(combined_result)
