#!/usr/bin/env Rscript
# Wang 2020, continuous response BMI. All model fits use installed package APIs.
# Reuse make_tree_objects(), fit_rare_path() and cv_rare() from the supplied helper.

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

# This script lives in code/run code/real data study/wang_2020_continuous/.
project_root <- Sys.getenv(
  "TREEFA_PROJECT_ROOT",
  unset = file.path(script_dir(), "..", "..", "..", "..")
)
project_root <- normalizePath(project_root, mustWork = TRUE)
code_dir <- Sys.getenv(
  "TREEFA_CLUSTER_CODE_DIR",
  unset = file.path(project_root, "code", "core_code")
)
code_dir <- normalizePath(code_dir, mustWork = TRUE)
work_dir <- Sys.getenv("TREEFA_STUDY_DIR", unset = script_dir())
work_dir <- normalizePath(work_dir, mustWork = TRUE)
data_dir <- Sys.getenv(
  "TREEFA_REAL_DATA_DIR",
  unset = file.path(project_root, "data")
)
out_dir <- Sys.getenv(
  "TREEFA_RESULT_DIR",
  unset = file.path(project_root, "output", "real data study", "wang_2020_continuous")
)
helper_file <- file.path(code_dir, "rare_logistic_helpers.R")
data_file <- file.path(data_dir, "wang_2020_continuous_data.RData")
r_lib <- Sys.getenv("TREEFA_R_LIB", unset = file.path(code_dir, "rlib"))
.libPaths(c(r_lib, .libPaths()))

uu <- suppressWarnings(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
nreps <- 200L
weight.order <- -0.5
nfolds <- 5L
nval <- 31L
treeFA_Mmratio <- 10000
rare_Mmratio <- 10000  # cv_rare() uses 50 lambda values with min/max ratio 1e-4.
if (length(uu) != 1L || !is.finite(uu) || uu != floor(uu) || uu < 1 || uu > nreps) {
  stop("SLURM_ARRAY_TASK_ID must be an integer in 1:200.", call. = FALSE)
}

required_packages <- c("treeFA", "rare", "Matrix", "glmnet")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
required_tree_api <- c("real_data_one_round.linear", "tree_leaf_ids",
                       "gather_leaf_nodes_per_non_leaf")
if (!all(required_tree_api %in% getNamespaceExports("treeFA")) ||
    !all(c("stop_rule", "coarest_rule") %in%
         names(formals(treeFA::real_data_one_round.linear)))) {
  stop("The installed treeFA package lacks the revised linear API. ",
       "Install the supplied treeFA package before submitting.", call. = FALSE)
}
if (!all(c("y", "X", "A", "alpha", "lambda", "intercept", "maxite") %in%
         names(formals(rare::rarefit)))) {
  stop("The installed rare::rarefit has an incompatible interface.", call. = FALSE)
}
required_files <- c(helper_file, data_file)
missing_inputs <- required_files[!file.exists(required_files)]
if (length(missing_inputs)) {
  stop("Missing input file(s): ", paste(missing_inputs, collapse = "; "), call. = FALSE)
}
rare_helpers <- new.env(parent = globalenv())
source(helper_file, local = rare_helpers)
required_helpers <- c("log_msg", "make_tree_objects", "fit_rare_path", "cv_rare")
missing_helpers <- required_helpers[
  !vapply(required_helpers, exists, logical(1), envir = rare_helpers,
          mode = "function", inherits = FALSE)
]
if (length(missing_helpers)) {
  stop("Missing function(s) in ", helper_file, ": ",
       paste(missing_helpers, collapse = ", "), call. = FALSE)
}
log_msg <- rare_helpers$log_msg

# Read the original response directly from data$y; keep the supplied row order.
data_env <- new.env(parent = emptyenv())
load(data_file, envir = data_env)
if (!exists("data", envir = data_env, inherits = FALSE) ||
    !all(c("X", "y", "tree_df") %in% names(data_env$data))) {
  stop("wang_2020_continuous_data.RData must contain data$X, data$y and data$tree_df.", call. = FALSE)
}
X <- as.matrix(data_env$data$X)
y <- data_env$data$y
if (!is.numeric(X) || !is.numeric(y) || length(y) != nrow(X) ||
    any(!is.finite(X)) || any(!is.finite(y)) || nrow(X) - nval < 2L * nfolds) {
  stop("X and BMI must be finite numeric data with matching row counts and enough training rows.",
       call. = FALSE)
}
# Reuse the supplied helper for the same full-tree weights and RARE expansion.
tree_objects <- rare_helpers$make_tree_objects(
  tree_df = data_env$data$tree_df, p = ncol(X), weight.order = weight.order
)
tree_df <- tree_objects$tree_df
A <- tree_objects$A

log_msg("Task=", uu, "; weight=", weight.order, "; outcome=BMI")
log_msg("RARE helpers: ", helper_file)
log_msg("treeFA ", as.character(utils::packageVersion("treeFA")), " at ", find.package("treeFA"))
log_msg("rare ", as.character(utils::packageVersion("rare")), " at ", find.package("rare"))
log_msg("Fitting treeFA, LASSO and Ridge through the new treeFA package")
base_result <- treeFA::real_data_one_round.linear(
  y = y, X = X, tree_df = tree_df, split.seed = uu, nfolds = nfolds, nval = nval,
  ridge_param = 0, Mmratio = treeFA_Mmratio, max_iter = 10000L,
  stop_rule = "coef", coarest_rule = "revised"
)
train.index <- base_result$train.index
test.index <- base_result$test.index
X_train <- X[train.index, , drop = FALSE]
X_test <- X[test.index, , drop = FALSE]
y_train <- y[train.index]
y_test <- y[test.index]
set.seed(2 * uu)
foldid <- as.integer(cut(sample(seq_along(y_train)), breaks = nfolds, labels = FALSE))
rare_result <- rare_helpers$cv_rare(
  y = y_train, X = X_train, A = A, foldid = foldid
)
rare_loss <- mean((y_test - as.numeric(X_test %*% rare_result$beta))^2)

# Retain the grouped-LASSO comparison: sum columns within treeFA coefficient groups.
# Refit at this grouped model's own selected lambda.
group_id <- as.integer(factor(base_result$our.beta))
group_count <- length(unique(group_id))
new_lasso_beta <- new_lasso_loss <- new_lasso_param <- NA_real_
new_lasso_cv <- NULL
if (group_count > 1L) {
  H <- matrix(0, nrow = ncol(X), ncol = group_count)
  H[cbind(seq_len(ncol(X)), group_id)] <- 1
  X_group_train <- X_train %*% H
  X_group_test <- X_test %*% H
  set.seed(2 * uu)
  new_lasso_cv <- glmnet::cv.glmnet(
    X_group_train, y_train, family = "gaussian", alpha = 1,
    intercept = FALSE, standardize = TRUE, nfolds = nfolds
  )
  new_lasso_param <- new_lasso_cv$lambda.min
  new_lasso_fit <- glmnet::glmnet(
    X_group_train, y_train, family = "gaussian", alpha = 1,
    intercept = FALSE, standardize = TRUE, lambda = new_lasso_param
  )
  if (new_lasso_fit$jerr != 0L || length(new_lasso_fit$lambda) != 1L) {
    stop("The selected grouped-LASSO refit did not finish.", call. = FALSE)
  }
  new_lasso_beta <- as.numeric(new_lasso_fit$beta)
  new_lasso_loss <- mean((y_test - as.numeric(X_group_test %*% new_lasso_beta))^2)
}

# Keep established result keys for existing summaries; 'our' denotes treeFA.
result <- list(
  lasso.beta = base_result$lasso.beta, lasso.loss = base_result$lasso.loss,
  lasso.param = base_result$lasso.param,
  ridge.beta = base_result$ridge.beta, ridge.loss = base_result$ridge.loss,
  ridge.param = base_result$ridge.param,
  rare.beta = rare_result$beta, rare.loss = rare_loss, rare.param = rare_result$selected.param,
  our.beta = base_result$our.beta, our.loss = base_result$our.loss,
  our.param = base_result$our.param,
  new.lasso.beta = new_lasso_beta, new.lasso.loss = new_lasso_loss,
  new.lasso.param = new_lasso_param
)
config <- list(
  analysis = "wang_2020_continuous", outcome = "BMI", family = "gaussian",
  split.seed = uu, replication = uu, nreps = nreps, weight.order = weight.order,
  nfolds = nfolds, nval = nval, data_file = data_file, response_source = "data$y",
  treeFA_Mmratio = treeFA_Mmratio, rare_Mmratio = rare_Mmratio,
  code_dir = code_dir, helper_file = helper_file, treeFA_max_iter = 10000L, rare_maxite = 1000000L,
  treeFA_stop_rule = "coef", treeFA_coarest_rule = "revised", rare_alpha = 1,
  packages = setNames(vapply(required_packages, function(pkg)
    as.character(utils::packageVersion(pkg)), character(1)), required_packages)
)
cv_details <- list(RARE = rare_result,
                   treeFA = list(lambda = base_result$our.lambda,
                                 valid.error = base_result$our.valid.error),
                   new.lasso = if (is.null(new_lasso_cv)) NULL else
                     list(lambda = new_lasso_cv$lambda, cvm = new_lasso_cv$cvm))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create result directory: ", out_dir)
savename <- file.path(out_dir, paste0(
  "wang_2020_continuous_weight", weight.order, "_uu", uu, ".RData"
))
session_info <- sessionInfo()
save(result, config, cv_details, train.index, test.index, foldid, tree_df, A,
     session_info, file = savename)
log_msg("Saved ", savename)
