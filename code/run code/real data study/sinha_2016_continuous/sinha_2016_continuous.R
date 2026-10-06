#!/usr/bin/env Rscript
# Sinha 2016, continuous response BMI. All model fits use installed package APIs.
# No external R scripts are sourced.

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

# This script lives in code/run code/real data study/sinha_2016_continuous/.
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
  unset = file.path(project_root, "output", "real data study", "sinha_2016_continuous")
)
data_file <- file.path(data_dir, "sinha_2016_data.RData")
metadata_file <- file.path(data_dir, "sinha_2016_metadata.tsv")
r_lib <- Sys.getenv("TREEFA_R_LIB", unset = file.path(code_dir, "rlib"))
.libPaths(c(r_lib, .libPaths()))

uu <- suppressWarnings(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
nreps <- 200L
weight.order <- -0.5
nfolds <- 5L
nval <- 31L
if (length(uu) != 1L || !is.finite(uu) || uu != floor(uu) || uu < 1 || uu > nreps) {
  stop("SLURM_ARRAY_TASK_ID must be an integer in 1:200.", call. = FALSE)
}

required_packages <- c("treeFA", "rare", "Matrix", "glmnet", "readr")
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
missing_inputs <- c(data_file, metadata_file)[!file.exists(c(data_file, metadata_file))]
if (length(missing_inputs)) {
  stop("Missing input file(s): ", paste(missing_inputs, collapse = "; "), call. = FALSE)
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

# The supplied RData contains data$X and data$tree_df. Preserve metadata row order.
data_env <- new.env(parent = emptyenv())
load(data_file, envir = data_env)
if (!exists("data", envir = data_env, inherits = FALSE) ||
    !all(c("X", "tree_df") %in% names(data_env$data))) {
  stop("sinha_2016_data.RData must contain data$X and data$tree_df.", call. = FALSE)
}
metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE)
if (!"BMI" %in% names(metadata)) stop("The metadata have no BMI column.", call. = FALSE)
X <- as.matrix(data_env$data$X)
y <- metadata$BMI
if (!is.numeric(X) || !is.numeric(y) || length(y) != nrow(X) ||
    any(!is.finite(X)) || any(!is.finite(y)) || nrow(X) - nval < 2L * nfolds) {
  stop("X and BMI must be finite numeric data with matching row counts and enough training rows.",
       call. = FALSE)
}
tree_objects <- make_tree_objects(data_env$data$tree_df, ncol(X), weight.order)
tree_df <- tree_objects$tree_df
A <- tree_objects$A

log_msg("Task=", uu, "; weight=", weight.order, "; outcome=BMI")
log_msg("treeFA ", as.character(utils::packageVersion("treeFA")), " at ", find.package("treeFA"))
log_msg("rare ", as.character(utils::packageVersion("rare")), " at ", find.package("rare"))
log_msg("Fitting treeFA, LASSO and Ridge through the new treeFA package")
base_result <- treeFA::real_data_one_round.linear(
  y = y, X = X, tree_df = tree_df, split.seed = uu, nfolds = nfolds, nval = nval,
  ridge_param = 0, Mmratio = 10000, max_iter = 10000L,
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
rare_result <- cv_rare(y_train, X_train, A, foldid)
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
  analysis = "sinha_2016_continuous", outcome = "BMI", family = "gaussian",
  split.seed = uu, replication = uu, nreps = nreps, weight.order = weight.order,
  nfolds = nfolds, nval = nval, data_file = data_file, metadata_file = metadata_file,
  code_dir = code_dir, treeFA_max_iter = 10000L, rare_maxite = 1000000L,
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
  "sinha_2016_continuous_weight", weight.order, "_uu", uu, ".RData"
))
session_info <- sessionInfo()
save(result, config, cv_details, train.index, test.index, foldid, tree_df, A,
     session_info, file = savename)
log_msg("Saved ", savename)
