#!/usr/bin/env Rscript

cat("\n=== Experiment5: logistic Bien tree p-increasing simulation, treeFA + RARE ===\n\n")
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

elapsed_now <- function() {
  unname(proc.time()[["elapsed"]])
}

make_runtime_row <- function(
  p_value,
  k_value,
  rep_id,
  method,
  method_order,
  runtime_seconds,
  n_tuning_values,
  status = "ok"
) {
  data.frame(
    variant = paste0("p_", p_value),
    p = as.integer(p_value),
    k = as.integer(k_value),
    within_call_rep = as.integer(rep_id),
    method = as.character(method),
    method_order = as.integer(method_order),
    runtime_seconds = as.numeric(runtime_seconds),
    n_tuning_values = as.integer(n_tuning_values),
    timing_source = "R wall time: fit path + observed-validation selection",
    status = as.character(status),
    stringsAsFactors = FALSE
  )
}

experiment_dir <- Sys.getenv("TREEFA_EXPERIMENT5_DIR", unset = script_dir())
experiment_dir <- normalizePath(experiment_dir, mustWork = TRUE)

sim_dir <- Sys.getenv("TREEFA_SIM_DIR", unset = experiment_dir)
sim_dir <- normalizePath(sim_dir, mustWork = TRUE)

# This script lives in code/run code/simulations/experiment5/.
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

source(file.path(code_dir, "other_functions_treeFA.R"))
source(file.path(code_dir, "rare_logistic_helpers.R"))

required <- c("treeFA", "mclust", "glmnet", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing required R package(s): ", paste(missing, collapse = ", "),
    ". Install treeFA and dependencies in TREEFA_R_LIB before submitting.",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(treeFA))
if (!"grid.logistic" %in% getNamespaceExports("treeFA")) {
  stop(
    "Loaded treeFA does not export grid.logistic. Reinstall the updated treeFA package in TREEFA_R_LIB before running experiment5.",
    call. = FALSE
  )
}
if (!exists("rarefit.logistic", mode = "function")) {
  stop("rarefit.logistic was not found. Check that rare_logistic_helpers.R exists in code_dir.", call. = FALSE)
}

simulate_data_binary_experiment5 <- function(n, group.index, s = 0, beta.pre = NULL) {
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

  beta <- as.numeric(A %*% beta0)
  X <- matrix(stats::rpois(n * p, 0.02), nrow = n, ncol = p)
  means <- as.numeric(X %*% beta)
  prob <- stats::plogis(means)
  Y <- stats::rbinom(n = n, size = 1, prob = prob)

  list(X = X, Y = Y, true.beta = beta, A = A, prob = prob)
}

simulate_error_valid_logistic_experiment5 <- function(
  n,
  n0,
  p,
  k,
  n1 = n,
  s = 0,
  reps = 1,
  ridge.param = 0,
  thresh = 1e-5,
  weight.order = -1 / 2,
  max_iter = 1e6,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy")
) {
  n1 <- n
  stop_rule <- match.arg(stop_rule)
  coarest_rule <- match.arg(coarest_rule)

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
  selected_betas <- vector("list", reps)
  true_betas <- vector("list", reps)
  tree_groups <- vector("list", reps)
  runtime_rows <- vector("list", reps)

  for (rep_id in seq_len(reps)) {
    cat(rep_id, "/", reps, "\n")
    flush.console()

    trees <- simulate_tree(k, p)
    tree_df <- hclust_to_df(trees$tree, weight.order = weight.order)
    data <- simulate_data_binary_experiment5(n = n + n0 + n1, group.index = trees$group, s = s)

    train_index <- sample(seq_len(n + n0 + n1), n)
    valid_index <- sample(setdiff(seq_len(n + n0 + n1), train_index), n1)
    test_index <- setdiff(seq_len(n + n0 + n1), c(train_index, valid_index))

    X_train <- data$X[train_index, , drop = FALSE]
    Y_train <- data$Y[train_index]
    X_valid <- data$X[valid_index, , drop = FALSE]
    Y_valid <- data$Y[valid_index]
    X_test <- data$X[test_index, , drop = FALSE]
    Y_test <- data$Y[test_index]
    # Evaluate the original logistic loss against the noiseless response probability.
    prob_test <- data$prob[test_index]

    rare_time_start <- elapsed_now()
    rare.result <- rarefit.logistic(
      y = Y_train,
      X = X_train,
      tree_df = tree_df,
      intercept = FALSE
    )
    rare.valid.loss <- rare_negtv_lglkh(Y_valid, X_valid, rare.result$beta)
    rare.best <- which.min(rare.valid.loss)
    rare.beta <- rare.result$beta[, rare.best]
    rare_runtime <- elapsed_now() - rare_time_start

    rare.ideal.loss <- rare_negtv_lglkh(
      stats::plogis(as.numeric(X_valid %*% data$true.beta)),
      X_valid,
      rare.result$beta
    )
    rare.beta.ideal <- rare.result$beta[, which.min(rare.ideal.loss)]
    rare.minloss[rep_id] <- rare_negtv_lglkh(prob_test, X_test, rare.beta)
    rare.minloss.ideal[rep_id] <- rare_negtv_lglkh(prob_test, X_test, rare.beta.ideal)

    our_time_start <- elapsed_now()
    our.result <- treeFA::grid.logistic(
      Y = Y_train,
      X = X_train,
      tree_df = tree_df,
      true_beta = data$true.beta,
      ridge.param = ridge.param,
      thresh = thresh,
      max_iter = max_iter,
      stop_rule = stop_rule,
      coarest_rule = coarest_rule
    )
    our.valid.loss <- rare_negtv_lglkh(Y_valid, X_valid, our.result$beta)
    our.best <- which.min(our.valid.loss)
    our.beta <- our.result$beta[, our.best]
    our_runtime <- elapsed_now() - our_time_start

    our.ideal.loss <- rare_negtv_lglkh(
      stats::plogis(as.numeric(X_valid %*% data$true.beta)),
      X_valid,
      our.result$beta
    )
    our.beta.ideal <- our.result$beta[, which.min(our.ideal.loss)]
    our.minloss[rep_id] <- rare_negtv_lglkh(prob_test, X_test, our.beta)
    our.minloss.ideal[rep_id] <- rare_negtv_lglkh(prob_test, X_test, our.beta.ideal)

    runtime_rows[[rep_id]] <- rbind(
      make_runtime_row(
        p_value = p,
        k_value = k,
        rep_id = rep_id,
        method = "RARE",
        method_order = 1L,
        runtime_seconds = rare_runtime,
        n_tuning_values = ncol(rare.result$beta)
      ),
      make_runtime_row(
        p_value = p,
        k_value = k,
        rep_id = rep_id,
        method = "treeFA",
        method_order = 2L,
        runtime_seconds = our_runtime,
        n_tuning_values = ncol(our.result$beta)
      )
    )

    selected_betas[[rep_id]] <- cbind(
      treeFA = as.numeric(our.beta),
      RARE = as.numeric(rare.beta)
    )
    true_betas[[rep_id]] <- as.numeric(data$true.beta)
    tree_groups[[rep_id]] <- as.integer(trees$group)

    rare.rand[rep_id] <- mclust::adjustedRandIndex(
      as.numeric(as.factor(data$true.beta)),
      as.numeric(as.factor(rare.beta))
    )
    rare.rand.ideal[rep_id] <- mclust::adjustedRandIndex(
      as.numeric(as.factor(data$true.beta)),
      as.numeric(as.factor(rare.beta.ideal))
    )
    our.rand[rep_id] <- mclust::adjustedRandIndex(
      as.numeric(as.factor(data$true.beta)),
      as.numeric(as.factor(our.beta))
    )
    our.rand.ideal[rep_id] <- mclust::adjustedRandIndex(
      as.numeric(as.factor(data$true.beta)),
      as.numeric(as.factor(our.beta.ideal))
    )

    ls.res <- glmnet::glmnet(
      y = Y_train,
      x = X_train,
      lambda = 0,
      intercept = FALSE,
      nlambda = 50
    )
    ls.coef <- as.numeric(ls.res$beta)
    ls.loss[rep_id] <- rare_negtv_lglkh(prob_test, X_test, ls.coef)

    Q <- make_new_X(data$X, trees$group)
    ols.X <- Q$X1[train_index, , drop = FALSE]
    ols.res <- glmnet::glmnet(
      y = Y_train,
      x = ols.X,
      lambda = 0,
      intercept = FALSE,
      nlambda = 50
    )
    ols.coef <- as.numeric(ols.res$beta)
    oracle.ls.loss[rep_id] <- rare_negtv_lglkh(prob_test, Q$X1[test_index, , drop = FALSE], ols.coef)

    fit_ridge <- glmnet::glmnet(
      X_train,
      Y_train,
      alpha = 0,
      intercept = FALSE,
      family = "binomial",
      nlambda = 50
    )
    val_mse_ridge <- rare_negtv_lglkh(Y_valid, X_valid, as.matrix(fit_ridge$beta))
    ridge.beta <- fit_ridge$beta[, which.min(val_mse_ridge)]
    ridge.loss[rep_id] <- rare_negtv_lglkh(prob_test, X_test, ridge.beta)

    or.X_train <- Q$X1[train_index, , drop = FALSE]
    or.X_valid <- Q$X1[valid_index, , drop = FALSE]
    or.X_test <- Q$X1[test_index, , drop = FALSE]

    fit_oracle_ridge <- glmnet::glmnet(
      or.X_train,
      Y_train,
      alpha = 0,
      intercept = FALSE,
      family = "binomial",
      nlambda = 50
    )
    val_mse_oracle_ridge <- rare_negtv_lglkh(Y_valid, or.X_valid, as.matrix(fit_oracle_ridge$beta))
    or.ridge.beta <- fit_oracle_ridge$beta[, which.min(val_mse_oracle_ridge)]
    oracle.ridge.loss[rep_id] <- rare_negtv_lglkh(prob_test, or.X_test, or.ridge.beta)

    # Keep the existing oracle fits and save their coefficients in feature coordinates.
    selected_betas[[rep_id]] <- cbind(
      selected_betas[[rep_id]],
      oLS = as.numeric(Q$Q %*% ols.coef),
      oRidge = as.numeric(Q$Q %*% or.ridge.beta)
    )
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
    oracle.ridge.loss = oracle.ridge.loss,
    selected_betas = selected_betas,
    true_betas = true_betas,
    tree_groups = tree_groups,
    runtime_result = do.call(rbind, runtime_rows)
  )
}

uu_text <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = Sys.getenv("TREEFA_UU", unset = "1"))
uu <- as.numeric(uu_text)
if (!is.finite(uu) || uu < 1) stop("uu must be a positive number. Got: ", uu_text, call. = FALSE)

nreps <- as.integer(Sys.getenv("EXPERIMENT5_NREPS", unset = "200"))
n <- as.integer(Sys.getenv("EXPERIMENT5_N", unset = "50"))
n0 <- as.integer(Sys.getenv("EXPERIMENT5_N0", unset = as.character(n * 10L)))
n1 <- as.integer(Sys.getenv("EXPERIMENT5_N1", unset = as.character(n)))
k <- as.integer(Sys.getenv("EXPERIMENT5_K", unset = "20"))
ps <- parse_numeric_vector(Sys.getenv("EXPERIMENT5_PS", unset = ""), c(50, 100, 200, 400, 600, 800, 1000))
ps <- as.integer(ps)
s <- as.numeric(Sys.getenv("EXPERIMENT5_S", unset = "0"))
ridge.param <- as.numeric(Sys.getenv("EXPERIMENT5_RIDGE_PARAM", unset = "0"))
thresh <- as.numeric(Sys.getenv("TREEFA_THRESH", unset = "1e-5"))
max_iter <- as.integer(Sys.getenv("TREEFA_MAX_ITER", unset = "1000000"))

weight_grid <- -1 / 2
weight_index <- ((uu - 1) %/% nreps) + 1L
if (weight_index > length(weight_grid)) {
  stop("uu=", uu, " is outside the configured array range 1:", nreps * length(weight_grid), call. = FALSE)
}
weight.order <- weight_grid[weight_index]
uu_rep <- (uu - 1) %% nreps + 1L

task_stop_rule <- Sys.getenv("TREEFA_STOP_RULE", unset = "objective")
task_coarest_rule <- Sys.getenv("TREEFA_COAREST_RULE", unset = "revised")

experiment_name <- Sys.getenv("TREEFA_EXPERIMENT_NAME", unset = "experiment5")
out_dir <- Sys.getenv("TREEFA_RESULT_DIR", unset = file.path(project_root, "output", experiment_name))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Experiment dir:", experiment_dir, "\n")
cat("Simulation dir:", sim_dir, "\n")
cat("Code dir:", code_dir, "\n")
cat("R lib:", r_lib, "\n")
cat("Result dir:", out_dir, "\n")
cat("uu =", uu, "; replicate =", uu_rep, "; weight.order =", weight.order, "\n")
cat("n =", n, "; n0 =", n0, "; n1 =", n1, "; k =", k, "; ps =", paste(ps, collapse = ","), "\n")
cat("treeFA stop_rule =", task_stop_rule, "; coarest_rule =", task_coarest_rule, "; thresh =", thresh, "; max_iter =", max_iter, "\n\n")
flush.console()

set.seed(uu)

metric_names <- c(
  "our.minloss", "rare.minloss", "our.rand", "rare.rand",
  "our.minloss.ideal", "rare.minloss.ideal",
  "our.rand.ideal", "rare.rand.ideal",
  "oracle.ls.loss", "ls.loss", "ridge.loss", "oracle.ridge.loss"
)

result <- data.frame()
selected_betas <- list()
true_betas <- list()
tree_groups <- list()
runtime_rows <- list()
for (idx in seq_along(ps)) {
  p_value <- ps[idx]
  label <- paste0("p_", p_value)
  start_time <- proc.time()[["elapsed"]]
  cat("\n--- p =", p_value, "(", idx, "/", length(ps), ") ---\n")
  flush.console()

  single.output <- simulate_error_valid_logistic_experiment5(
    n = n,
    n0 = n0,
    p = p_value,
    k = k,
    n1 = n1,
    s = s,
    reps = 1,
    ridge.param = ridge.param,
    thresh = thresh,
    max_iter = max_iter,
    weight.order = weight.order,
    stop_rule = task_stop_rule,
    coarest_rule = task_coarest_rule
  )

  single.result <- as.data.frame(single.output[metric_names])
  selected_betas[[label]] <- single.output$selected_betas[[1L]]
  true_betas[[label]] <- single.output$true_betas[[1L]]
  tree_groups[[label]] <- single.output$tree_groups[[1L]]
  runtime_rows[[label]] <- single.output$runtime_result

  result <- rbind(result, single.result)
  elapsed <- proc.time()[["elapsed"]] - start_time
  cat(idx, ",", p_value, "is finished; elapsed seconds =", round(elapsed, 2), "\n")
  print(single.result)
  cat("Runtime by method (seconds):\n")
  print(single.output$runtime_result[, c("method", "runtime_seconds", "n_tuning_values", "method_order")])
  flush.console()
}

rownames(result) <- ps

runtime_result <- do.call(rbind, runtime_rows)
rownames(runtime_result) <- NULL
runtime_result$replication <- uu_rep
runtime_result$weight.order <- weight.order
runtime_result <- runtime_result[, c(
  "replication", "weight.order", "variant", "p", "k", "within_call_rep",
  "method", "runtime_seconds", "n_tuning_values", "timing_source",
  "method_order", "status"
)]

runtime_table <- stats::reshape(
  runtime_result[, c(
    "replication", "weight.order", "variant", "p", "k", "within_call_rep",
    "method", "runtime_seconds"
  )],
  idvar = c("replication", "weight.order", "variant", "p", "k", "within_call_rep"),
  timevar = "method",
  direction = "wide"
)
rownames(runtime_table) <- NULL

out_file <- file.path(
  out_dir,
  paste0(
    "Valid_BienTree_Logistic_n50_pIncre_k20_weightChange_weight",
    weight.order,
    "_uu",
    uu_rep,
    ".RData"
  )
)

save(
  result,
  selected_betas,
  true_betas,
  tree_groups,
  runtime_result,
  runtime_table,
  uu,
  uu_rep,
  ps,
  n,
  n0,
  n1,
  k,
  weight.order,
  file = out_file
)

cat("\n=== Experiment5 finished ===\n")
cat("Saved result to:\n  ", out_file, "\n", sep = "")
cat("Result dimensions:", paste(dim(result), collapse = " x "), "\n")
print(result)
cat("\nRuntime by method (seconds):\n")
print(runtime_result)
