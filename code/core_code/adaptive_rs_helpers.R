# Adaptive gamma-grid wrapper for RS-DL2 and RS-CL2.

make_log_gamma_grid <- function(min_gamma, max_gamma, length.out, decreasing = TRUE) {
  if (!is.finite(min_gamma) || !is.finite(max_gamma) || min_gamma <= 0 || max_gamma <= 0) {
    stop("Gamma bounds must be positive finite values.", call. = FALSE)
  }
  lo <- min(min_gamma, max_gamma)
  hi <- max(min_gamma, max_gamma)
  grid <- exp(seq(log(hi), log(lo), length.out = length.out))
  if (!decreasing) grid <- rev(grid)
  grid
}

rs_selected_boundary <- function(rs_fit, method, edge = 2L) {
  method_idx <- match(method, rs_fit$summary$method)
  if (is.na(method_idx)) stop("Method not found in RS summary: ", method, call. = FALSE)

  selected <- rs_fit$summary$gamma[method_idx]
  cv_file <- file.path(rs_fit$output_dir, paste0(method, "_cv.csv"))
  if (!file.exists(cv_file) || !is.finite(selected) || selected <= 0) {
    return(list(side = "unknown", selected = selected, index = NA_integer_, n_grid = NA_integer_))
  }

  cv <- utils::read.csv(cv_file, header = FALSE)
  gamma_values <- as.numeric(cv[[1]])
  gamma_values <- gamma_values[is.finite(gamma_values) & gamma_values > 0]
  if (!length(gamma_values)) {
    return(list(side = "unknown", selected = selected, index = NA_integer_, n_grid = NA_integer_))
  }

  ordered <- sort(gamma_values)
  selected_rank <- which.min(abs(log(ordered) - log(selected)))
  edge <- min(as.integer(edge), length(ordered))

  side <- "interior"
  if (selected_rank <= edge) side <- "lower"
  if (selected_rank >= length(ordered) - edge + 1L) side <- "upper"

  list(side = side, selected = selected, index = selected_rank, n_grid = length(ordered))
}

expand_gamma_range <- function(gamma_range, side, factor = 10, length.out = length(gamma_range)) {
  gamma_range <- as.numeric(gamma_range)
  gamma_range <- gamma_range[is.finite(gamma_range) & gamma_range > 0]
  if (!length(gamma_range)) stop("gamma_range must contain positive finite values.", call. = FALSE)

  min_gamma <- min(gamma_range)
  max_gamma <- max(gamma_range)
  if (identical(side, "lower")) min_gamma <- min_gamma / factor
  if (identical(side, "upper")) max_gamma <- max_gamma * factor

  make_log_gamma_grid(min_gamma, max_gamma, length.out = length.out, decreasing = TRUE)
}

safe_method_name <- function(method) {
  gsub("[^A-Za-z0-9]+", "_", method)
}

extract_rs_method_fit <- function(rs_fit, method) {
  method_idx <- match(method, rs_fit$summary$method)
  if (is.na(method_idx)) stop("Method not found in RS summary: ", method, call. = FALSE)
  list(
    summary = rs_fit$summary[method_idx, , drop = FALSE],
    beta = rs_fit$beta[, method],
    fit = rs_fit
  )
}

run_rs_method_adaptive <- function(
  method,
  X_train,
  y_train,
  X_valid,
  y_valid,
  X_test,
  y_test,
  tree_df,
  base_gamma_range,
  work_dir,
  seed = 1,
  selection = "validation",
  model = "article_aligned",
  normalize_rows = FALSE,
  mu = NULL,
  adaptive = TRUE,
  edge = 2L,
  expand_factor = 10,
  max_expansions = 3L,
  grid_length = length(base_gamma_range),
  verbose = TRUE
) {
  gamma_range <- as.numeric(base_gamma_range)
  history <- list()
  best <- NULL

  max_round <- if (adaptive) as.integer(max_expansions) else 0L
  for (round_id in 0:max_round) {
    round_dir <- file.path(work_dir, paste0(safe_method_name(method), "_round", round_id))
    rs_fit <- run_rs_matlab_once(
      X_train = X_train,
      y_train = y_train,
      X_valid = X_valid,
      y_valid = y_valid,
      X_test = X_test,
      y_test = y_test,
      tree_df = tree_df,
      tree_hclust = NULL,
      work_dir = round_dir,
      seed = seed,
      gamma_range = gamma_range,
      selection = selection,
      model = model,
      normalize_rows = normalize_rows,
      mu = mu,
      methods = method,
      keep_files = TRUE
    )

    method_fit <- extract_rs_method_fit(rs_fit, method)
    boundary <- rs_selected_boundary(rs_fit, method, edge = edge)
    row <- method_fit$summary
    history[[length(history) + 1L]] <- data.frame(
      method = method,
      round = round_id,
      min_gamma = min(gamma_range),
      max_gamma = max(gamma_range),
      selected_gamma = row$gamma,
      valid_mse = row$valid_mse,
      boundary = boundary$side,
      selected_index = boundary$index,
      n_grid = boundary$n_grid,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat(
        "  ", method, " adaptive round ", round_id,
        ": selected gamma = ", signif(row$gamma, 5),
        ", valid_mse = ", signif(row$valid_mse, 5),
        ", boundary = ", boundary$side, "\n",
        sep = ""
      )
    }

    if (is.null(best) || row$valid_mse < best$summary$valid_mse[1]) {
      best <- method_fit
    }
    if (!adaptive || identical(boundary$side, "interior") || identical(boundary$side, "unknown")) break
    gamma_range <- expand_gamma_range(
      gamma_range,
      side = boundary$side,
      factor = expand_factor,
      length.out = grid_length
    )
  }

  best$history <- do.call(rbind, history)
  best
}

run_rs_matlab_adaptive_once <- function(
  X_train,
  y_train,
  X_valid,
  y_valid,
  X_test,
  y_test,
  tree_df,
  work_dir,
  seed = 1,
  gamma_range,
  selection = "validation",
  model = "article_aligned",
  normalize_rows = FALSE,
  mu = NULL,
  adaptive_methods = c("RS-DL2", "RS-CL2"),
  edge = 2L,
  expand_factor = 10,
  max_expansions = 3L,
  grid_length = length(gamma_range),
  verbose = TRUE
) {
  method_order <- c("RS-DL2", "RS-CL2", "RS-L1")
  method_fits <- vector("list", length(method_order))
  names(method_fits) <- method_order

  for (method in method_order) {
    method_fits[[method]] <- run_rs_method_adaptive(
      method = method,
      X_train = X_train,
      y_train = y_train,
      X_valid = X_valid,
      y_valid = y_valid,
      X_test = X_test,
      y_test = y_test,
      tree_df = tree_df,
      base_gamma_range = gamma_range,
      work_dir = work_dir,
      seed = seed,
      selection = selection,
      model = model,
      normalize_rows = normalize_rows,
      mu = mu,
      adaptive = method %in% adaptive_methods,
      edge = edge,
      expand_factor = expand_factor,
      max_expansions = max_expansions,
      grid_length = grid_length,
      verbose = verbose
    )
  }

  summary <- do.call(rbind, lapply(method_fits, `[[`, "summary"))
  rownames(summary) <- NULL
  beta <- do.call(cbind, lapply(method_fits, `[[`, "beta"))
  colnames(beta) <- method_order
  adaptive_history <- do.call(rbind, lapply(method_fits, `[[`, "history"))
  rownames(adaptive_history) <- NULL

  list(
    summary = summary,
    beta = beta,
    adaptive_history = adaptive_history,
    method_fits = method_fits,
    work_dir = work_dir
  )
}
