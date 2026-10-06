# Cross-validation for the separate joint aggregation/selection interfaces.

.selection_cv_foldid <- function(Y, folds, foldid, family) {
  n <- length(Y)
  if (is.null(foldid)) {
    folds <- .selection_path_integer(folds, "folds", minimum = 2L)
    if (folds > n) stop("folds cannot exceed the number of observations.", call. = FALSE)
    if (family == "binomial") {
      counts <- tabulate(as.integer(Y) + 1L, nbins = 2L)
      if (any(counts < 2L)) {
        stop("Binomial cross-validation requires at least two observations in each class.", call. = FALSE)
      }
      foldid <- integer(n)
      offset <- 0L
      for (value in c(0, 1)) {
        indices <- which(Y == value)
        shuffled <- indices[sample.int(length(indices))]
        foldid[shuffled] <- (seq_along(shuffled) - 1L + offset) %% folds + 1L
        offset <- (offset + length(shuffled)) %% folds
      }
    } else {
      foldid <- sample(rep(seq_len(folds), length.out = n))
    }
  } else {
    if (!is.numeric(foldid) || is.complex(foldid) || !is.null(dim(foldid)) || length(foldid) != n ||
        any(!is.finite(foldid)) || any(foldid < 1) || any(foldid != floor(foldid)) ||
        any(foldid > .Machine$integer.max)) {
      stop("foldid must contain one finite positive integer per observation.", call. = FALSE)
    }
    unique_ids <- sort(unique(foldid))
    if (length(unique_ids) < 2L || !identical(as.numeric(unique_ids), as.numeric(seq_along(unique_ids)))) {
      stop("foldid must contain at least two nonempty folds numbered consecutively from 1.", call. = FALSE)
    }
    folds <- length(unique_ids)
  }
  foldid <- as.integer(foldid)
  sizes <- tabulate(foldid, nbins = folds)
  if (any(sizes == 0L) || any(sizes >= n)) {
    stop("Every fold must have nonempty training and validation sets.", call. = FALSE)
  }
  if (family == "binomial") {
    for (i in seq_len(folds)) {
      if (length(unique(Y[foldid != i])) != 2L) {
        stop("Every binomial training fold must contain both response classes.", call. = FALSE)
      }
    }
  }
  foldid
}

.selection_cv <- function(Y, X, tree_df, param_grid, lambda_seq, tau_seq, family,
                          intercept, ridge_param, step_size, abs_tol, rel_tol,
                          max_iter, warm_start, folds, foldid, retry_max_iter, verbose) {
  param_grid <- .selection_parameter_grid(param_grid, lambda_seq, tau_seq)
  max_iter <- .selection_path_controls(abs_tol, rel_tol, max_iter, warm_start)
  if (is.null(retry_max_iter)) retry_max_iter <- min(2 * as.double(max_iter), .Machine$integer.max)
  retry_max_iter <- .selection_path_integer(retry_max_iter, "retry_max_iter")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || is.complex(X) || !is.numeric(Y) || is.complex(Y) ||
      nrow(X) < 2L || ncol(X) < 1L || length(Y) != nrow(X) ||
      any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("Cross-validation requires finite numeric X and Y with matching lengths, at least two observations, and one feature.",
         call. = FALSE)
  }
  Y <- as.numeric(Y)
  if (family == "binomial") .check_binary_y(Y)
  foldid <- .selection_cv_foldid(Y, folds, foldid, family)
  folds <- max(foldid)
  fold_sizes <- tabulate(foldid, nbins = folds)
  n_candidates <- nrow(param_grid)
  fold_error <- matrix(NA_real_, nrow = folds, ncol = n_candidates,
                       dimnames = list(as.character(seq_len(folds)), as.character(seq_len(n_candidates))))
  fold_diagnostics <- vector("list", folds)
  names(fold_diagnostics) <- as.character(seq_len(folds))

  for (i in seq_len(folds)) {
    if (verbose) message("fold ", i, "/", folds)
    train <- foldid != i
    test <- !train
    prepared <- .selection_prepare_model(
      Y[train], X[train, , drop = FALSE], tree_df, family = family,
      intercept = intercept, ridge_param = ridge_param, step_size = step_size
    )
    path <- .selection_run_path(prepared, param_grid, abs_tol, rel_tol, max_iter,
                                warm_start, retry_max_iter = retry_max_iter)
    eta <- sweep(X[test, , drop = FALSE] %*% path$beta, 2L, path$beta0, "+")
    if (family == "gaussian") {
      heldout_loss <- colMeans((Y[test] - eta)^2)
    } else {
      # log(1 + exp(eta)) - y*eta, evaluated without subtracting two
      # potentially large numbers when y is one.
      signed_eta <- sweep(eta, 1L, 1 - 2 * Y[test], "*")
      heldout_loss <- colMeans(pmax(signed_eta, 0) + log1p(exp(-abs(signed_eta))))
    }
    score_finite <- colSums(!is.finite(eta)) == 0L & is.finite(heldout_loss)
    available <- path$diagnostics$converged & score_finite
    fold_error[i, available] <- heldout_loss[available]
    path$diagnostics$score_finite <- score_finite
    fold_diagnostics[[i]] <- path$diagnostics
  }

  eligible <- colSums(is.finite(fold_error)) == folds
  cv_error <- rep(Inf, n_candidates)
  if (!any(eligible)) {
    stop("No parameter candidate converged with a finite validation loss in every fold; increase max_iter or retry_max_iter.",
         call. = FALSE)
  }
  cv_error[eligible] <- colSums(sweep(fold_error[, eligible, drop = FALSE], 1L, fold_sizes, "*")) / length(Y)
  selected <- which.min(cv_error)
  if (any(!eligible)) {
    warning(sum(!eligible), " parameter candidate(s) were excluded because at least one fold did not converge or had a non-finite validation loss.",
            call. = FALSE)
  }
  if (verbose) message("refitting candidate ", selected, " on all observations")
  prepared <- .selection_prepare_model(Y, X, tree_df, family = family,
                                       intercept = intercept, ridge_param = ridge_param,
                                       step_size = step_size)
  fit <- .selection_fit_with_retry(
    prepared, lambda = param_grid$lambda[selected], tau = param_grid$tau[selected],
    abs_tol = abs_tol, rel_tol = rel_tol, max_iter = max_iter,
    retry_max_iter = retry_max_iter
  )
  if (!isTRUE(fit$converged)) {
    stop("The selected parameter candidate did not converge when refitted on all observations; increase max_iter or retry_max_iter.",
         call. = FALSE)
  }
  list(selected.lambda = param_grid$lambda[selected], selected.tau = param_grid$tau[selected],
       selected.index = selected, beta = fit$beta, beta0 = fit$beta0, fit = fit,
       param_grid = param_grid, foldid = foldid, fold.error = fold_error, cv.error = cv_error,
       eligible = eligible, fold.diagnostics = fold_diagnostics)
}

#' Cross-validation for joint tree aggregation and feature selection
#'
#' Every fold uses the same explicit parameter grid. Centering and step-size
#' calculations use only training observations. Validation scores are squared
#' error for the linear model and stable negative binomial log-likelihood for
#' the logistic model, averaged per observation, including for unequal folds.
#'
#' Nonconverged fits are retried once from their last splitting state with
#' retry_max_iter additional iterations. Candidates unavailable in any fold get
#' cv.error = Inf and eligible = FALSE; they are never scored on only successful
#' folds. The first row wins exact ties. The selected full-data refit must converge.
#'
#' @inheritParams grid.select_linear
#' @param folds Number of folds, at least two, when foldid is not supplied.
#' @param foldid Optional positive consecutive integer fold labels starting at one.
#'   When supplied, its distinct labels determine the fold count and folds is ignored.
#' @param retry_max_iter Additional iterations for one continuation of a failed
#'   fit; NULL defaults to twice max_iter. Also applied to the full-data refit.
#' @param verbose Whether to print fold and refit progress.
#' @return A list containing selected.lambda, selected.tau, selected.index, beta,
#'   beta0, the converged full-data fit, param_grid, foldid, fold.error (fold by
#'   candidate), cv.error, eligible, and fold.diagnostics. Diagnostic iter counts
#'   include both attempts when a retry was needed. Default logistic folds are
#'   stratified, and every training fold must contain both response classes.
#' @export
cv.select_linear <- function(Y, X, tree_df, param_grid = NULL, lambda_seq = NULL,
                             tau_seq = NULL, intercept = FALSE, ridge_param = 0,
                             step_size = NULL, abs_tol = 1e-7, rel_tol = 1e-6,
                             max_iter = 10000L, warm_start = TRUE, folds = 5L,
                             foldid = NULL, retry_max_iter = NULL, verbose = FALSE) {
  .selection_cv(Y, X, tree_df, param_grid, lambda_seq, tau_seq, "gaussian",
                intercept, ridge_param, step_size, abs_tol, rel_tol, max_iter, warm_start,
                folds, foldid, retry_max_iter, verbose)
}

#' @rdname cv.select_linear
#' @export
cv.select_logistic <- function(Y, X, tree_df, param_grid = NULL, lambda_seq = NULL,
                               tau_seq = NULL, intercept = FALSE, ridge_param = 0,
                               step_size = NULL, abs_tol = 1e-7, rel_tol = 1e-6,
                               max_iter = 10000L, warm_start = TRUE, folds = 5L,
                               foldid = NULL, retry_max_iter = NULL, verbose = FALSE) {
  .selection_cv(Y, X, tree_df, param_grid, lambda_seq, tau_seq, "binomial",
                intercept, ridge_param, step_size, abs_tol, rel_tol, max_iter, warm_start,
                folds, foldid, retry_max_iter, verbose)
}
