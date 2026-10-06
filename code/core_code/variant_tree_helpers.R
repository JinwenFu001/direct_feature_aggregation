# Helpers for running RARE and RS penalties on the zero-weight tree variants.

find_project_root_variant <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "treeFA")) &&
        (dir.exists(file.path(current, "RS-code")) || dir.exists(file.path(current, "RS-Code")))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Cannot find project root with treeFA/ and RS-code or RS-Code/.", call. = FALSE)
}

collapse_zero_weight_tree <- function(tree_df, zero_tol = 0, keep_root = TRUE) {
  tree_df <- as.data.frame(tree_df)
  required <- c("node", "parent", "weight")
  missing <- setdiff(required, names(tree_df))
  if (length(missing)) {
    stop("tree_df is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  tree_df$node <- as.integer(tree_df$node)
  tree_df$parent <- as.integer(tree_df$parent)
  parent_nodes <- unique(tree_df$parent[!is.na(tree_df$parent)])
  leaves <- setdiff(tree_df$node, parent_nodes)
  root <- tree_df$node[is.na(tree_df$parent)]
  if (length(root) != 1L) stop("tree_df must contain exactly one root node.", call. = FALSE)

  is_zero_internal <- tree_df$node %in% parent_nodes &
    !is.na(tree_df$weight) &
    abs(tree_df$weight) <= zero_tol
  keep <- tree_df$node %in% leaves | !is_zero_internal
  if (keep_root) keep[tree_df$node == root] <- TRUE
  keep_nodes <- tree_df$node[keep]

  parent_by_node <- stats::setNames(tree_df$parent, tree_df$node)
  nearest_kept_parent <- function(node) {
    parent <- parent_by_node[[as.character(node)]]
    while (!is.na(parent) && !(parent %in% keep_nodes)) {
      parent <- parent_by_node[[as.character(parent)]]
    }
    as.integer(parent)
  }

  out <- tree_df[tree_df$node %in% keep_nodes, , drop = FALSE]
  out$parent <- vapply(out$node, function(node) {
    if (node == root) return(NA_integer_)
    nearest_kept_parent(node)
  }, integer(1))

  out <- out[order(out$node), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "removed_zero_nodes") <- sort(setdiff(tree_df$node, keep_nodes))
  out
}

tree_variant_summary <- function(original_tree_df, collapsed_tree_df) {
  original_parent_nodes <- unique(original_tree_df$parent[!is.na(original_tree_df$parent)])
  collapsed_parent_nodes <- unique(collapsed_tree_df$parent[!is.na(collapsed_tree_df$parent)])
  data.frame(
    original_nodes = nrow(original_tree_df),
    collapsed_nodes = nrow(collapsed_tree_df),
    original_internal = length(original_parent_nodes),
    collapsed_internal = length(collapsed_parent_nodes),
    removed_zero_nodes = length(attr(collapsed_tree_df, "removed_zero_nodes")),
    stringsAsFactors = FALSE
  )
}

make_sparse_expansion_from_tree_df <- function(tree_df, p) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }
  A <- rs_tree_df_to_expansion_matrix(tree_df, p = p)
  methods::as(Matrix::Matrix(A, sparse = TRUE), "dgCMatrix")
}
