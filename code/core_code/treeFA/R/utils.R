.require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but is not installed.", pkg), call. = FALSE)
  }
  invisible(TRUE)
}

.check_required_columns <- function(df, cols, arg_name = "df") {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(arg_name, " is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_tree_df <- function(tree_df, require_weight = FALSE) {
  required_cols <- c("node", "parent")
  if (isTRUE(require_weight)) required_cols <- c(required_cols, "weight")
  .check_required_columns(tree_df, required_cols, "tree_df")

  if (anyDuplicated(tree_df$node)) {
    stop("tree_df$node must contain unique node IDs.", call. = FALSE)
  }
  parent_values <- tree_df$parent[!is.na(tree_df$parent)]
  if (length(parent_values) > 0L && !all(parent_values %in% tree_df$node)) {
    stop("Every non-NA parent in tree_df$parent must also appear in tree_df$node.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_yx <- function(Y, X) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  if (length(Y) != nrow(X)) {
    stop("length(Y) must equal nrow(X).", call. = FALSE)
  }
  if (!all(is.finite(Y)) || !all(is.finite(X))) {
    stop("Y and X must contain only finite values.", call. = FALSE)
  }
  list(Y = Y, X = X)
}
