# Joint aggregation/selection paths. These interfaces deliberately do not change
# the existing treeFA paths or their return conventions.

.selection_parameter_grid <- function(param_grid, lambda_seq, tau_seq) {
  if (!is.null(param_grid)) {
    if (!is.null(lambda_seq) || !is.null(tau_seq)) {
      stop("Supply param_grid or lambda_seq and tau_seq, not both.", call. = FALSE)
    }
    if (!is.data.frame(param_grid) && !is.matrix(param_grid)) {
      stop("param_grid must be a data frame or matrix with lambda and tau columns.", call. = FALSE)
    }
    param_grid <- as.data.frame(param_grid)
    if (anyDuplicated(names(param_grid)) || !all(c("lambda", "tau") %in% names(param_grid))) {
      stop("param_grid must have uniquely named lambda and tau columns.", call. = FALSE)
    }
    if (!nrow(param_grid)) stop("param_grid must contain at least one row.", call. = FALSE)
  } else {
    if (is.null(lambda_seq) || is.null(tau_seq)) {
      stop("Supply an explicit param_grid or both lambda_seq and tau_seq.", call. = FALSE)
    }
    for (name in c("lambda_seq", "tau_seq")) {
      value <- get(name)
      if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) || !length(value) ||
          any(!is.finite(value)) || any(value < 0)) {
        stop(name, " must be a nonempty numeric vector of finite non-negative values.", call. = FALSE)
      }
    }
    param_grid <- expand.grid(lambda = as.numeric(lambda_seq), tau = as.numeric(tau_seq),
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }
  for (name in c("lambda", "tau")) {
    value <- param_grid[[name]]
    if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) || any(!is.finite(value)) || any(value < 0)) {
      stop("param_grid$", name, " must contain finite non-negative numeric values.", call. = FALSE)
    }
    param_grid[[name]] <- as.numeric(value)
  }
  rownames(param_grid) <- NULL
  param_grid
}

.selection_path_integer <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || is.complex(value) || length(value) != 1L || !is.finite(value) ||
      value < minimum || value != floor(value) || value > .Machine$integer.max) {
    stop(name, " must be a finite integer at least ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

.selection_path_controls <- function(abs_tol, rel_tol, max_iter, warm_start) {
  for (name in c("abs_tol", "rel_tol")) {
    value <- get(name)
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L || !is.finite(value) || value < 0) {
      stop(name, " must be a finite non-negative scalar.", call. = FALSE)
    }
  }
  if (abs_tol == 0 && rel_tol == 0) {
    stop("At least one of abs_tol and rel_tol must be positive.", call. = FALSE)
  }
  if (!is.logical(warm_start) || length(warm_start) != 1L || is.na(warm_start)) {
    stop("warm_start must be TRUE or FALSE.", call. = FALSE)
  }
  .selection_path_integer(max_iter, "max_iter")
}

.selection_fit_with_retry <- function(prepared, lambda, tau, abs_tol, rel_tol,
                                      max_iter, init_z = NULL, retry_max_iter = NULL) {
  fit <- .selection_fit_prepared(
    prepared, lambda = lambda, tau = tau, abs_tol = abs_tol, rel_tol = rel_tol,
    max_iter = max_iter, init_z = init_z, keep_history = FALSE, warn = FALSE
  )
  fit$attempts <- 1L
  if (!isTRUE(fit$converged) && !is.null(retry_max_iter)) {
    first_iter <- fit$iter
    fit <- .selection_fit_prepared(
      prepared, lambda = lambda, tau = tau, abs_tol = abs_tol, rel_tol = rel_tol,
      max_iter = retry_max_iter, init_z = fit$state$z, keep_history = FALSE, warn = FALSE
    )
    fit$iter <- as.double(first_iter) + fit$iter
    fit$attempts <- 2L
  }
  fit
}

.selection_run_path <- function(prepared, param_grid, abs_tol, rel_tol, max_iter,
                                warm_start, retry_max_iter = NULL) {
  n_candidates <- nrow(param_grid)
  beta <- NULL
  beta0 <- numeric(n_candidates)
  diagnostics <- data.frame(
    converged = rep(FALSE, n_candidates), iter = numeric(n_candidates),
    status = rep(NA_character_, n_candidates), residual = rep(NA_real_, n_candidates),
    residual_scaled = rep(NA_real_, n_candidates), objective = rep(NA_real_, n_candidates),
    attempts = integer(n_candidates), stringsAsFactors = FALSE
  )
  z <- NULL
  for (i in seq_len(n_candidates)) {
    fit <- .selection_fit_with_retry(
      prepared, lambda = param_grid$lambda[i], tau = param_grid$tau[i],
      abs_tol = abs_tol, rel_tol = rel_tol, max_iter = max_iter,
      init_z = if (warm_start) z else NULL, retry_max_iter = retry_max_iter
    )
    if (is.null(beta)) {
      beta <- matrix(0, nrow = length(fit$beta), ncol = n_candidates,
                     dimnames = list(names(fit$beta), NULL))
    }
    beta[, i] <- fit$beta
    beta0[i] <- fit$beta0
    for (name in names(diagnostics)) diagnostics[[name]][i] <- fit[[name]]
    z <- fit$state$z
  }
  list(beta = beta, beta0 = beta0, param_grid = param_grid, diagnostics = diagnostics,
       n_nonzero = colSums(beta != 0))
}

.selection_grid <- function(Y, X, tree_df, param_grid, lambda_seq, tau_seq, family,
                            intercept, ridge_param, step_size, abs_tol, rel_tol,
                            max_iter, warm_start) {
  param_grid <- .selection_parameter_grid(param_grid, lambda_seq, tau_seq)
  max_iter <- .selection_path_controls(abs_tol, rel_tol, max_iter, warm_start)
  prepared <- .selection_prepare_model(Y, X, tree_df, family = family,
                                       intercept = intercept, ridge_param = ridge_param,
                                       step_size = step_size)
  path <- .selection_run_path(prepared, param_grid, abs_tol, rel_tol, max_iter, warm_start)
  if (any(!path$diagnostics$converged)) {
    warning(sum(!path$diagnostics$converged),
            " parameter candidate(s) did not converge; inspect diagnostics.", call. = FALSE)
  }
  path
}

#' Joint tree aggregation and feature selection paths
#'
#' Fit an explicitly specified collection of (lambda, tau) pairs. Rows and
#' duplicates in param_grid are preserved. With sequences, lambda varies fastest
#' in the Cartesian product. Tree preprocessing and the step-size bound are
#' computed once for the whole path. Warm starts carry the splitting state z.
#' Neither function standardizes columns of X.
#'
#' @param Y Response vector, binary 0/1 for the logistic model.
#' @param X Numeric design matrix.
#' @param tree_df Tree data frame with node, parent, and weight columns.
#' @param param_grid Data frame or matrix with finite non-negative lambda and tau columns.
#' @param lambda_seq,tau_seq Both numeric sequences, used instead of param_grid.
#' @param intercept Whether to fit an unpenalized intercept.
#' @param ridge_param Non-negative ridge parameter, with penalty ridge_param * sum(beta^2)/(2*n).
#' @param step_size Optional fixed step size; otherwise a safe value is computed.
#' @param abs_tol,rel_tol Absolute and relative fixed-point stopping tolerances.
#' @param max_iter Maximum number of iterations per candidate.
#' @param warm_start Whether to initialize each candidate using the previous z.
#' @return A list with beta (one column per candidate), beta0, param_grid,
#'   diagnostics, and n_nonzero. Diagnostics include convergence, iteration count,
#'   status, fixed-point residuals, full objective, and number of attempts.
#' @export
grid.select_linear <- function(Y, X, tree_df, param_grid = NULL, lambda_seq = NULL,
                               tau_seq = NULL, intercept = FALSE, ridge_param = 0,
                               step_size = NULL, abs_tol = 1e-7, rel_tol = 1e-6,
                               max_iter = 10000L, warm_start = TRUE) {
  .selection_grid(Y, X, tree_df, param_grid, lambda_seq, tau_seq, "gaussian",
                  intercept, ridge_param, step_size, abs_tol, rel_tol, max_iter, warm_start)
}

#' @rdname grid.select_linear
#' @export
grid.select_logistic <- function(Y, X, tree_df, param_grid = NULL, lambda_seq = NULL,
                                 tau_seq = NULL, intercept = FALSE, ridge_param = 0,
                                 step_size = NULL, abs_tol = 1e-7, rel_tol = 1e-6,
                                 max_iter = 10000L, warm_start = TRUE) {
  .selection_grid(Y, X, tree_df, param_grid, lambda_seq, tau_seq, "binomial",
                  intercept, ridge_param, step_size, abs_tol, rel_tol, max_iter, warm_start)
}
