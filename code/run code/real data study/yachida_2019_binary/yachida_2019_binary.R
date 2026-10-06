#!/usr/bin/env Rscript
# Yachida 2019 binary response (CRC): LASSO, Ridge, RARE and treeFA.
# Use the revised treeFA package and code/core_code/rare_logistic_helpers.R.

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

# This script lives in code/run code/real data study/yachida_2019_binary/.
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
  unset = file.path(project_root, "output", "real data study", "yachida_2019_binary")
)
data_file <- file.path(data_dir, "yachida_2019_binary_data.RData")
helper_file <- file.path(code_dir, "rare_logistic_helpers.R")
r_lib <- Sys.getenv("TREEFA_R_LIB", unset = file.path(code_dir, "rlib"))
.libPaths(c(r_lib, .libPaths()))

uu <- suppressWarnings(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
nreps <- 200L
weight.order <- -0.5
nfolds <- 5L
nval <- 50L
treeFA_Mmratio <- 100
rare_Mmratio <- 10000  # The original RARE grid is independent of treeFA Mmratio.
if (length(uu) != 1L || !is.finite(uu) || uu != floor(uu) || uu < 1 || uu > nreps) {
  stop("SLURM_ARRAY_TASK_ID must be an integer in 1:200.", call. = FALSE)
}
required_packages <- c("treeFA", "Matrix", "glmnet", "pROC")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
required_tree_api <- c("cv.logistic", "negtv_lglkh", "tree_leaf_ids",
                       "gather_leaf_nodes_per_non_leaf")
if (!all(required_tree_api %in% getNamespaceExports("treeFA")) ||
    !all(c("stop_rule", "coarest_rule") %in% names(formals(treeFA::cv.logistic)))) {
  stop("Install the revised treeFA package with the current logistic API.", call. = FALSE)
}
required_files <- c(data_file, helper_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing input file(s): ", paste(missing_files, collapse = "; "), call. = FALSE)
}
rare_helpers <- new.env(parent = globalenv())
source(helper_file, local = rare_helpers)
required_helpers <- c("cv.rare.logistic", "rarefit.logistic", "df_to_A")
if (!all(vapply(required_helpers, exists, logical(1), envir = rare_helpers,
                mode = "function", inherits = FALSE))) {
  stop("The revised rare_logistic_helpers.R is missing required functions.", call. = FALSE)
}

data_env <- new.env(parent = emptyenv())
load(data_file, envir = data_env)
if (!exists("data", envir = data_env, inherits = FALSE) ||
    !all(c("X", "y", "tree_df") %in% names(data_env$data))) {
  stop("yachida_2019_binary_data.RData must contain data$X, data$y and data$tree_df.", call. = FALSE)
}
X <- as.matrix(data_env$data$X)
y <- data_env$data$y
if (!is.numeric(X) || !is.numeric(y) || length(y) != nrow(X) ||
    any(!is.finite(X)) || any(!is.finite(y)) || !all(y %in% c(0, 1)) ||
    nrow(X) - nval < 2L * nfolds) {
  stop("X must be finite numeric data and y must be matching numeric 0/1 labels.",
       call. = FALSE)
}
y <- as.numeric(y)

# Preserve the normalized full-tree weights using the revised treeFA helpers.
tree_df <- as.data.frame(data_env$data$tree_df)
if (!identical(as.integer(treeFA::tree_leaf_ids(tree_df)), seq_len(ncol(X)))) {
  stop("The tree leaves must be numbered 1:ncol(X).", call. = FALSE)
}
tree_df$weight <- 0
groups <- treeFA::gather_leaf_nodes_per_non_leaf(tree_df)$df_result
if (!nrow(groups) || sum(groups$depth == 0) != 1L) {
  stop("This analysis requires a full tree with one internal root.", call. = FALSE)
}
raw_weights <- lengths(groups$leaves)^weight.order
tree_df$weight[match(groups$node, tree_df$node)] <- raw_weights / mean(raw_weights)

# Preserve the held-out split and each method's original CV seed and fold rule.
set.seed(uu)
test.index <- sample(seq_len(nrow(X)), nval)
train.index <- setdiff(seq_len(nrow(X)), test.index)
X_train <- X[train.index, , drop = FALSE]
X_test <- X[test.index, , drop = FALSE]
y_train <- y[train.index]
y_test <- y[test.index]
if (length(unique(y_train)) != 2L || length(unique(y_test)) != 2L) {
  stop("The fixed train/test split must contain both classes in each set.", call. = FALSE)
}
set.seed(2 * uu)
tree_rare_foldid <- as.integer(cut(sample(seq_along(y_train)), breaks = nfolds,
                                 labels = FALSE))
if (any(vapply(seq_len(nfolds), function(k)
  length(unique(y_train[tree_rare_foldid != k])) != 2L, logical(1)))) {
  stop("A fixed RARE/treeFA CV training fold has only one class.", call. = FALSE)
}
cat("Task=", uu, "; outcome=CRC; weight=", weight.order,
    "; train/test=", length(y_train), "/", length(y_test), "\n", sep = "")
cat("treeFA: ", as.character(utils::packageVersion("treeFA")), " at ",
    find.package("treeFA"), "\nRARE helpers: ", helper_file, "\n", sep = "")

baseline <- list()
baseline_alpha <- c(lasso = 1, ridge = 0)
for (key in names(baseline_alpha)) {
  cat("Fitting ", key, "\n", sep = "")
  set.seed(2 * uu)
  cv <- glmnet::cv.glmnet(
    X_train, y_train, family = "binomial", alpha = baseline_alpha[[key]],
    intercept = FALSE, standardize = TRUE, nfolds = nfolds, type.measure = "deviance"
  )
  fit <- glmnet::glmnet(
    X_train, y_train, family = "binomial", alpha = baseline_alpha[[key]],
    intercept = FALSE, standardize = TRUE, lambda = cv$lambda.min
  )
  if (fit$jerr != 0L || length(fit$lambda) != 1L) {
    stop("The selected ", key, " refit did not finish.", call. = FALSE)
  }
  baseline[[key]] <- list(beta = as.numeric(fit$beta), param = cv$lambda.min,
                         cv = list(lambda = cv$lambda, cvm = cv$cvm, cvsd = cv$cvsd))
}

cat("Fitting RARE with the shared logistic CV helper\n")
set.seed(2 * uu)
rare_cv <- rare_helpers$cv.rare.logistic(
  Y = y_train, X = X_train, tree_df = tree_df, folds = nfolds,
  thresh = 1e-5, intercept = FALSE, Mmratio = rare_Mmratio, verbose = TRUE
)
cat("Fitting treeFA with the revised package\n")
set.seed(2 * uu)
tree_cv <- treeFA::cv.logistic(
  Y = y_train, X = X_train, tree_df = tree_df, folds = nfolds,
  thresh = 1e-5, ridge.param = 0, intercept = FALSE, Mmratio = treeFA_Mmratio,
  max_iter = 10000L, stop_rule = "coef", coarest_rule = "revised", verbose = TRUE
)

# Preserve the original 16 result keys and their order; 'our' denotes treeFA.
betas <- list(lasso = baseline$lasso$beta, ridge = baseline$ridge$beta,
              rare = as.numeric(rare_cv$beta$beta), our = as.numeric(tree_cv$beta))
params <- c(lasso = baseline$lasso$param, ridge = baseline$ridge$param,
            rare = rare_cv$selected.param, our = tree_cv$selected.param)
result <- list()
auc_direction <- setNames(character(length(betas)), names(betas))
for (key in names(betas)) {
  beta <- betas[[key]]
  if (length(beta) != ncol(X) || any(!is.finite(beta))) {
    stop("Invalid final coefficients for ", key, ".", call. = FALSE)
  }
  probability <- as.numeric(plogis(X_test %*% beta))
  # Keep the original pROC automatic direction; pass a vector, not a matrix.
  roc_fit <- pROC::roc(response = y_test, predictor = probability,
                       levels = c(0, 1), direction = "auto", quiet = TRUE)
  result[[paste0(key, ".beta")]] <- beta
  result[[paste0(key, ".loss")]] <- as.numeric(treeFA::negtv_lglkh(y_test, X_test, beta))
  result[[paste0(key, ".auc")]] <- as.numeric(pROC::auc(roc_fit))
  result[[paste0(key, ".param")]] <- unname(params[[key]])
  auc_direction[[key]] <- roc_fit$direction
}

config <- list(
  analysis = "yachida_2019_binary", outcome = "CRC", family = "binomial",
  split.seed = uu, replication = uu, nreps = nreps, weight.order = weight.order,
  nval = nval, nfolds = nfolds, intercept = FALSE, ridge_param = 0,
  treeFA_Mmratio = treeFA_Mmratio, rare_Mmratio = rare_Mmratio,
  data_file = data_file, code_dir = code_dir, helper_file = helper_file,
  treeFA_max_iter = 10000L, rare_maxite = 1000000L,
  treeFA_stop_rule = "coef", treeFA_coarest_rule = "revised",
  rare_solver = "rarefit.logistic: glmnet binomial on X %*% A",
  loss = "mean negative binomial log-likelihood",
  auc_direction = "auto (original convention)", actual_auc_direction = auc_direction,
  packages = setNames(vapply(required_packages, function(pkg)
    as.character(utils::packageVersion(pkg)), character(1)), required_packages)
)
cv_details <- list(RARE = rare_cv, treeFA = tree_cv,
                   lasso = baseline$lasso$cv, ridge = baseline$ridge$cv)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create result directory: ", out_dir)
savename <- file.path(out_dir, paste0(
  "yachida_2019_binary_weight", weight.order, "_uu", uu, ".RData"
))
session_info <- sessionInfo()
save(result, config, cv_details, train.index, test.index, tree_rare_foldid,
     tree_df, session_info, file = savename)
cat("Saved: ", savename, "\n", sep = "")
