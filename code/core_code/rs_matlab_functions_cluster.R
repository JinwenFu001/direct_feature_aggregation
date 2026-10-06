# R wrappers for calling the original MATLAB relative-shift code.

find_project_root <- function() {
  env_root <- Sys.getenv("TREEFA_CLUSTER_CODE_DIR", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  wd <- normalizePath(getwd(), mustWork = TRUE)
  current <- wd
  repeat {
    has_rs <- dir.exists(file.path(current, "RS-code")) || dir.exists(file.path(current, "RS-Code"))
    has_treefa <- dir.exists(file.path(current, "treeFA")) ||
      file.exists(file.path(current, "treeFA_0.0.1.tar.gz"))
    if (has_rs && has_treefa) return(current)

    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }

  stop(
    "Cannot find cluster code root. Set TREEFA_CLUSTER_CODE_DIR to the directory containing RS-Code/ and treeFA/.",
    call. = FALSE
  )
}

rs_find_matlab <- function() {
  env_bin <- Sys.getenv("MATLAB_BIN", unset = "")
  candidates <- c(
    env_bin,
    Sys.which("matlab"),
    "/apps/matlab/bin/matlab",
    "/opt/matlab/bin/matlab",
    "/Applications/MATLAB_R2026a.app/bin/matlab",
    "/Applications/MATLAB_R2025b.app/bin/matlab",
    "/Applications/MATLAB_R2025a.app/bin/matlab",
    "/Applications/MATLAB_R2024b.app/bin/matlab",
    "/Applications/MATLAB_R2024a.app/bin/matlab"
  )
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop(
      "Cannot find MATLAB. Load the MATLAB module or set MATLAB_BIN to the full matlab executable path.",
      call. = FALSE
    )
  }
  normalizePath(candidates[1], mustWork = TRUE)
}

rs_matlab_quote <- function(x) {
  paste0("'", gsub("'", "''", normalizePath(x, mustWork = FALSE), fixed = TRUE), "'")
}

rs_write_matrix <- function(x, file) {
  x <- as.matrix(x)
  utils::write.table(x, file = file, sep = ",", row.names = FALSE, col.names = FALSE)
  invisible(file)
}

rs_normalize_rows <- function(X) {
  X <- as.matrix(X)
  row_sums <- rowSums(X)
  if (any(row_sums <= 0)) {
    warning("Some rows have non-positive sums; those rows are left unchanged.", call. = FALSE)
  }
  ok <- row_sums > 0
  X[ok, ] <- X[ok, , drop = FALSE] / row_sums[ok]
  X
}

rs_default_gamma_range <- function() {
  exp(seq(log(1e-4), log(5e-2), length.out = 50))
}

rs_rare_scaled_gamma_range <- function(
  X_train,
  y_train,
  intercept = FALSE,
  nlam = 50,
  lam.min.ratio = 1e-4
) {
  X_use <- as.matrix(X_train)
  y_use <- as.numeric(y_train)
  n <- nrow(X_use)

  if (intercept) {
    X_use <- sweep(X_use, 2, colMeans(X_use), "-")
    y_use <- y_use - mean(y_use)
  }

  lambda <- max(abs(as.numeric(crossprod(X_use, y_use)))) / n *
    exp(seq(0, log(lam.min.ratio), length.out = nlam))
  n * lambda
}

rs_tree_df_to_expansion_matrix <- function(tree_df, p = NULL) {
  tree_df <- as.data.frame(tree_df)
  required <- c("node", "parent")
  missing <- setdiff(required, names(tree_df))
  if (length(missing)) stop("tree_df is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)

  tree_df$node <- as.integer(tree_df$node)
  tree_df$parent <- as.integer(tree_df$parent)
  node_order <- sort(tree_df$node)
  if (is.null(p)) {
    parent_nodes <- unique(tree_df$parent[!is.na(tree_df$parent)])
    leaves <- sort(setdiff(tree_df$node, parent_nodes))
    p <- length(leaves)
  }
  leaves <- seq_len(p)
  root <- tree_df$node[is.na(tree_df$parent)]
  if (length(root) != 1L) stop("tree_df must contain exactly one root node.", call. = FALSE)
  if (!all(leaves %in% tree_df$node)) stop("Expected leaf nodes 1:p in tree_df.", call. = FALSE)

  parent_by_node <- stats::setNames(tree_df$parent, tree_df$node)
  A <- matrix(0, nrow = p, ncol = length(node_order))
  colnames(A) <- as.character(node_order)

  for (leaf in leaves) {
    current <- leaf
    visited <- integer()
    repeat {
      if (current %in% visited) stop("Cycle detected in tree_df.", call. = FALSE)
      visited <- c(visited, current)
      A[leaf, match(current, node_order)] <- 1
      parent <- parent_by_node[[as.character(current)]]
      if (is.na(parent)) break
      current <- as.integer(parent)
    }
  }

  A
}

rs_groups_to_penalty <- function(groups, q) {
  groups <- lapply(groups, function(g) sort(unique(as.integer(g))))
  groups <- groups[vapply(groups, length, integer(1)) > 0L]
  sizes <- vapply(groups, length, integer(1))
  total <- sum(sizes)

  C <- matrix(0, nrow = total, ncol = q)
  g_idx <- matrix(0, nrow = length(groups), ncol = 3)
  row_start <- 1L
  for (i in seq_along(groups)) {
    rows <- row_start:(row_start + sizes[i] - 1L)
    C[cbind(rows, groups[[i]])] <- 1
    g_idx[i, ] <- c(row_start, row_start + sizes[i] - 1L, sizes[i])
    row_start <- row_start + sizes[i]
  }

  c_norm <- if (total > 0L) max(colSums(C^2)) else 0
  list(C = C, g_idx = g_idx, CNorm = c_norm)
}

rs_method_to_penalty_keys <- function(methods = NULL) {
  method_map <- c("RS-DL2" = "RS_DL2", "RS-CL2" = "RS_CL2", "RS-L1" = "RS_L1")
  if (is.null(methods)) return(unname(method_map))

  methods <- unique(as.character(methods))
  methods <- methods[nzchar(methods)]
  keys <- ifelse(methods %in% names(method_map), unname(method_map[methods]), methods)
  bad <- setdiff(keys, unname(method_map))
  if (length(bad)) stop("Unknown RS method(s): ", paste(methods[keys %in% bad], collapse = ", "), call. = FALSE)
  keys
}

rs_article_penalty_matrices <- function(tree_df, A = NULL, p = NULL, methods = NULL) {
  tree_df <- as.data.frame(tree_df)
  required <- c("node", "parent")
  missing <- setdiff(required, names(tree_df))
  if (length(missing)) stop("tree_df is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)

  tree_df$node <- as.integer(tree_df$node)
  tree_df$parent <- as.integer(tree_df$parent)
  node_order <- sort(tree_df$node)
  root <- tree_df$node[is.na(tree_df$parent)]
  if (length(root) != 1L) stop("tree_df must contain exactly one root node.", call. = FALSE)

  if (is.null(A)) A <- rs_tree_df_to_expansion_matrix(tree_df, p = p)
  A <- as.matrix(A)
  if (ncol(A) != length(node_order)) {
    stop("A must have one column per tree node.", call. = FALSE)
  }

  q <- ncol(A)
  parent_nodes <- unique(tree_df$parent[!is.na(tree_df$parent)])
  internal_nodes <- sort(parent_nodes)
  children <- split(tree_df$node[!is.na(tree_df$parent)], tree_df$parent[!is.na(tree_df$parent)])

  get_descendants <- function(node) {
    frontier <- as.integer(children[[as.character(node)]])
    out <- integer()
    while (length(frontier)) {
      out <- c(out, frontier)
      next_frontier <- unlist(children[as.character(frontier)], use.names = FALSE)
      frontier <- as.integer(next_frontier)
    }
    out
  }

  to_columns <- function(nodes) {
    nodes <- sort(unique(nodes[nodes != root]))
    match(nodes, node_order)
  }

  penalty_keys <- rs_method_to_penalty_keys(methods)
  penalties <- list()

  if ("RS_DL2" %in% penalty_keys) {
    descendant_groups <- lapply(internal_nodes, function(node) to_columns(get_descendants(node)))
    penalties$RS_DL2 <- rs_groups_to_penalty(descendant_groups, q)
  }
  if ("RS_CL2" %in% penalty_keys) {
    child_groups <- lapply(internal_nodes, function(node) to_columns(as.integer(children[[as.character(node)]])))
    penalties$RS_CL2 <- rs_groups_to_penalty(child_groups, q)
  }
  if ("RS_L1" %in% penalty_keys) {
    l1_groups <- as.list(match(setdiff(node_order, root), node_order))
    penalties$RS_L1 <- rs_groups_to_penalty(l1_groups, q)
  }

  list(
    A = A,
    penalties = penalties,
    root_node = root,
    node_order = node_order
  )
}

rs_write_article_inputs <- function(input_dir, tree_df, tree_hclust = NULL, p = NULL, methods = NULL) {
  if (!is.null(tree_hclust)) {
    if (!requireNamespace("rare", quietly = TRUE)) {
      stop("Package 'rare' is required to build the article-aligned RS expansion matrix.", call. = FALSE)
    }
    A <- as.matrix(rare::tree.matrix(tree_hclust))
  } else {
    A <- NULL
  }

  article <- rs_article_penalty_matrices(tree_df = tree_df, A = A, p = p, methods = methods)
  rs_write_matrix(article$A, file.path(input_dir, "A_full.csv"))
  rs_write_matrix(article$root_node, file.path(input_dir, "root_node.csv"))
  rs_write_matrix(article$node_order, file.path(input_dir, "node_order.csv"))

  for (key in names(article$penalties)) {
    penalty <- article$penalties[[key]]
    rs_write_matrix(penalty$C, file.path(input_dir, paste0("C_", key, ".csv")))
    rs_write_matrix(penalty$g_idx, file.path(input_dir, paste0("gidx_", key, ".csv")))
    rs_write_matrix(penalty$CNorm, file.path(input_dir, paste0("cnorm_", key, ".csv")))
  }

  invisible(article)
}

rs_tree_df_to_taxonomy <- function(tree_df, p = NULL) {
  tree_df <- as.data.frame(tree_df)
  required <- c("node", "parent")
  missing <- setdiff(required, names(tree_df))
  if (length(missing)) stop("tree_df is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)

  tree_df$node <- as.integer(tree_df$node)
  if (is.null(p)) {
    parent_nodes <- unique(tree_df$parent[!is.na(tree_df$parent)])
    leaf_nodes <- sort(setdiff(tree_df$node, parent_nodes))
    p <- length(leaf_nodes)
  }
  leaves <- seq_len(p)
  if (!all(leaves %in% tree_df$node)) {
    stop("Expected leaf nodes 1:p to appear in tree_df$node.", call. = FALSE)
  }

  paths <- vector("list", p)
  for (i in seq_along(leaves)) {
    current <- leaves[i]
    path <- current
    visited <- integer()
    repeat {
      if (current %in% visited) stop("Cycle detected in tree_df.", call. = FALSE)
      visited <- c(visited, current)
      parent <- tree_df$parent[match(current, tree_df$node)]
      if (is.na(parent)) break
      current <- as.integer(parent)
      path <- c(path, current)
    }
    paths[[i]] <- path
  }

  max_len <- max(vapply(paths, length, integer(1)))
  taxonomy <- matrix(NA_integer_, nrow = p, ncol = max_len)
  for (i in seq_along(paths)) {
    path <- paths[[i]]
    taxonomy[i, seq_along(path)] <- path
    if (length(path) < max_len) {
      taxonomy[i, (length(path) + 1):max_len] <- path[length(path)]
    }
  }
  taxonomy
}

run_rs_matlab_once <- function(
  X_train,
  y_train,
  X_valid,
  y_valid,
  X_test,
  y_test,
  tree_df,
  work_dir = NULL,
  rs_code_dir = NULL,
  matlab_bin = NULL,
  gamma_range = rs_default_gamma_range(),
  seed = 1,
  selection = c("validation", "cv"),
  model = c("article_aligned", "paper_rootless"),
  tree_hclust = NULL,
  normalize_rows = FALSE,
  mu = NULL,
  methods = NULL,
  keep_files = TRUE
) {
  selection <- match.arg(selection)
  model <- match.arg(model)
  allowed_methods <- c("RS-DL2", "RS-CL2", "RS-L1")
  if (!is.null(methods)) {
    methods <- unique(as.character(methods))
    methods <- methods[nzchar(methods)]
    bad_methods <- setdiff(methods, allowed_methods)
    if (length(bad_methods)) {
      stop("Unknown RS method(s): ", paste(bad_methods, collapse = ", "), call. = FALSE)
    }
    if (!length(methods)) stop("methods must contain at least one RS method.", call. = FALSE)
  }
  root <- find_project_root()
  if (is.null(rs_code_dir)) {
    rs_code_dir <- file.path(root, "RS-Code")
    if (!dir.exists(rs_code_dir)) rs_code_dir <- file.path(root, "RS-code")
  }
  if (is.null(matlab_bin)) matlab_bin <- rs_find_matlab()
  if (is.null(work_dir)) {
    run_id <- paste0("rs_", format(Sys.time(), "%Y%m%d%H%M%OS6"), "_", Sys.getpid())
    run_id <- gsub("[^A-Za-z0-9_]", "", run_id)
    work_dir <- file.path(root, "test", "rs_runs", run_id)
  }

  input_dir <- file.path(work_dir, "input")
  output_dir <- file.path(work_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  X_train <- as.matrix(X_train)
  X_valid <- as.matrix(X_valid)
  X_test <- as.matrix(X_test)
  if (normalize_rows) {
    X_train <- rs_normalize_rows(X_train)
    X_valid <- rs_normalize_rows(X_valid)
    X_test <- rs_normalize_rows(X_test)
  }

  p <- ncol(X_train)
  taxonomy <- rs_tree_df_to_taxonomy(tree_df, p = p)
  if (identical(model, "article_aligned")) {
    rs_write_article_inputs(input_dir, tree_df = tree_df, tree_hclust = tree_hclust, p = p, methods = methods)
  }

  rs_write_matrix(X_train, file.path(input_dir, "X_train.csv"))
  rs_write_matrix(as.numeric(y_train), file.path(input_dir, "y_train.csv"))
  rs_write_matrix(X_valid, file.path(input_dir, "X_valid.csv"))
  rs_write_matrix(as.numeric(y_valid), file.path(input_dir, "y_valid.csv"))
  rs_write_matrix(X_test, file.path(input_dir, "X_test.csv"))
  rs_write_matrix(as.numeric(y_test), file.path(input_dir, "y_test.csv"))
  rs_write_matrix(taxonomy, file.path(input_dir, "taxonomy.csv"))
  rs_write_matrix(as.numeric(gamma_range), file.path(input_dir, "gamma_range.csv"))
  rs_write_matrix(seed, file.path(input_dir, "seed.csv"))
  writeLines(selection, file.path(input_dir, "selection.txt"))
  writeLines(model, file.path(input_dir, "model_mode.txt"))
  if (!is.null(methods)) writeLines(methods, file.path(input_dir, "method_list.txt"))
  if (!is.null(mu)) rs_write_matrix(mu, file.path(input_dir, "mu.csv"))

  matlab_cmd <- paste0(
    "addpath(genpath(", rs_matlab_quote(rs_code_dir), ")); ",
    "rs_fit_once(", rs_matlab_quote(input_dir), ", ", rs_matlab_quote(output_dir), ");"
  )
  matlab_args <- c("-batch", shQuote(matlab_cmd))
  matlab_command <- matlab_bin
  if (identical(Sys.info()[["sysname"]], "Darwin") &&
      file.exists(file.path(dirname(matlab_bin), "maca64"))) {
    matlab_args <- c("-maca64", matlab_args)
  }

  log_file <- file.path(work_dir, "matlab.log")
  cat("  MATLAB command:", matlab_command, paste(matlab_args, collapse = " "), "\n")
  cat("  MATLAB log file:", log_file, "\n")
  flush.console()

  log <- system2(matlab_command, matlab_args, stdout = log_file, stderr = log_file)
  status <- if (is.numeric(log)) log else attr(log, "status")
  if (is.null(status)) status <- 0L
  matlab_log <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character()
  if (!is.null(status) && status != 0) {
    stop("MATLAB RS run failed:\n", paste(matlab_log, collapse = "\n"), call. = FALSE)
  }

  summary_file <- file.path(output_dir, "rs_summary.csv")
  beta_file <- file.path(output_dir, "rs_betas.csv")
  if (!file.exists(summary_file) || !file.exists(beta_file)) {
    stop("MATLAB did not write expected RS output files.\n", paste(matlab_log, collapse = "\n"), call. = FALSE)
  }

  summary <- utils::read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  beta <- as.matrix(utils::read.csv(beta_file, check.names = FALSE))
  colnames(beta) <- summary$method

  if (identical(model, "article_aligned") && "RS-L1" %in% summary$method) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' is required for article-aligned exact RS-L1.", call. = FALSE)
    }

    A_full <- as.matrix(utils::read.csv(file.path(input_dir, "A_full.csv"), header = FALSE))
    lambda_l1 <- as.numeric(gamma_range) / nrow(X_train)
    fit_l1 <- glmnet::glmnet(
      X_train %*% A_full,
      as.numeric(y_train),
      lambda = lambda_l1,
      standardize = FALSE,
      intercept = FALSE,
      penalty.factor = c(rep(1, ncol(A_full) - 1L), 0),
      thresh = 1e-6,
      maxit = 1e6
    )
    beta_l1_path <- A_full %*% as.matrix(fit_l1$beta)
    valid_l1 <- colMeans((as.numeric(y_valid) - X_valid %*% beta_l1_path)^2)
    best_l1 <- which.min(valid_l1)
    beta_l1 <- beta_l1_path[, best_l1]

    l1_idx <- match("RS-L1", summary$method)
    summary$gamma[l1_idx] <- as.numeric(gamma_range)[best_l1]
    summary$valid_mse[l1_idx] <- valid_l1[best_l1]
    summary$test_mse[l1_idx] <- mean((as.numeric(y_test) - X_test %*% beta_l1)^2)
    summary$intercept[l1_idx] <- 0
    summary$status[l1_idx] <- "ok (article_aligned exact L1 from glmnet/RARE)"
    beta[, "RS-L1"] <- beta_l1
  }

  bad <- !grepl("^ok", summary$status)
  if (any(bad)) {
    warning(
      "Some RS methods failed: ",
      paste(summary$method[bad], summary$status[bad], sep = "=", collapse = "; "),
      call. = FALSE
    )
  }

  if (!keep_files) unlink(work_dir, recursive = TRUE, force = TRUE)

  list(
    summary = summary,
    beta = beta,
    matlab_log = matlab_log,
    matlab_log_file = log_file,
    work_dir = work_dir,
    input_dir = input_dir,
    output_dir = output_dir
  )
}
