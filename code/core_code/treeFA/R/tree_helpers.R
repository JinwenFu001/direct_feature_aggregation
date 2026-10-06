tree_leaf_ids <- function(tree_df) {
  .validate_tree_df(tree_df, require_weight = FALSE)
  sort(tree_df$node[!(tree_df$node %in% tree_df$parent[!is.na(tree_df$parent)])])
}

calculate_depth <- function(node, df) {
  .validate_tree_df(df, require_weight = FALSE)
  if (!node %in% df$node) stop("node is not present in df$node.", call. = FALSE)

  current <- node
  depth <- 0L
  visited <- character(0)

  repeat {
    current_key <- as.character(current)
    if (current_key %in% visited) stop("Cycle detected in tree_df.", call. = FALSE)
    visited <- c(visited, current_key)

    parent <- df$parent[match(current, df$node)]
    if (is.na(parent)) break

    current <- parent
    depth <- depth + 1L
  }

  depth
}

detect_non_leaf_nodes <- function(df) {
  .validate_tree_df(df, require_weight = FALSE)
  non_leaf_nodes <- unique(df$parent[!is.na(df$parent)])
  non_leaf_depths <- vapply(non_leaf_nodes, calculate_depth, integer(1), df = df)
  data.frame(node = non_leaf_nodes, depth = as.integer(non_leaf_depths))
}

gather_leaves <- function(node, df) {
  .validate_tree_df(df, require_weight = FALSE)
  if (!node %in% df$node) stop("node is not present in df$node.", call. = FALSE)

  leaves <- c()
  queue <- node
  visited <- character(0)

  while (length(queue) > 0L) {
    current <- queue[1L]
    queue <- queue[-1L]

    current_key <- as.character(current)
    if (current_key %in% visited) stop("Cycle detected in tree_df.", call. = FALSE)
    visited <- c(visited, current_key)

    children <- df$node[!is.na(df$parent) & df$parent == current]
    if (length(children) == 0L) {
      leaves <- c(leaves, current)
    } else {
      queue <- c(queue, children)
    }
  }

  sort(leaves)
}

gather_direct_latent_nodes <- function(node, df) {
  .validate_tree_df(df, require_weight = FALSE)
  if (!node %in% df$node) stop("node is not present in df$node.", call. = FALSE)

  children <- df$node[!is.na(df$parent) & df$parent == node]
  latent_children <- c()

  for (child in children) {
    grand_children <- df$node[!is.na(df$parent) & df$parent == child]
    if (length(grand_children) > 0L) latent_children <- c(latent_children, child)
  }

  sort(latent_children)
}

make_list <- function(group_df) {
  weight_list <- list()
  group_list <- list()
  children_list <- list()
  status_list <- list()

  if (!nrow(group_df)) {
    return(list(groups = group_list, weights = weight_list, children = children_list, status = status_list))
  }

  depth_values <- sort(unique(group_df$depth), decreasing = TRUE)

  for (i in seq_along(depth_values)) {
    temp_df <- group_df[group_df$depth == depth_values[i], , drop = FALSE]
    node_names <- as.character(temp_df$node)

    status_list[[i]] <- rep(1, length(temp_df$node))
    names(status_list[[i]]) <- node_names

    group_list[[i]] <- temp_df$leaves
    names(group_list[[i]]) <- node_names

    weight_list[[i]] <- temp_df$weight
    names(weight_list[[i]]) <- node_names

    children_list[[i]] <- temp_df$latent_children
    names(children_list[[i]]) <- node_names
  }

  list(groups = group_list, weights = weight_list, children = children_list, status = status_list)
}

gather_leaf_nodes_per_non_leaf <- function(df) {
  .validate_tree_df(df, require_weight = TRUE)
  non_leaf_nodes_df <- detect_non_leaf_nodes(df)

  if (!nrow(non_leaf_nodes_df)) {
    df_result <- data.frame(
      node = integer(),
      depth = integer(),
      leaves = I(list()),
      latent_children = I(list()),
      weight = numeric()
    )
    return(list(list_result = make_list(df_result), df_result = df_result))
  }

  non_leaf_nodes_df$leaves <- I(lapply(non_leaf_nodes_df$node, gather_leaves, df = df))
  non_leaf_nodes_df$latent_children <- I(lapply(non_leaf_nodes_df$node, gather_direct_latent_nodes, df = df))
  non_leaf_nodes_df$weight <- vapply(
    non_leaf_nodes_df$node,
    function(node) df$weight[match(node, df$node)],
    numeric(1)
  )

  result <- make_list(non_leaf_nodes_df)
  list(list_result = result, df_result = non_leaf_nodes_df)
}

find_p <- function(tree_df) {
  tree_df <- as.data.frame(tree_df)
  .validate_tree_df(tree_df, require_weight = FALSE)
  leaves <- tree_leaf_ids(tree_df)
  if (!length(leaves)) return(0L)

  expected <- seq_len(length(leaves))
  if (!all(leaves == expected)) {
    stop("Leaf nodes must be indexed as 1:p to match the columns of X.", call. = FALSE)
  }

  length(leaves)
}

.add_singleton_rows <- function(coarest_set, template, p) {
  covered <- if (nrow(coarest_set)) {
    sort(unique(unlist(coarest_set$leaves, use.names = FALSE)))
  } else {
    integer(0)
  }
  covered <- covered[!is.na(covered)]
  indiv_nodes <- setdiff(seq_len(p), covered)

  if (length(indiv_nodes) == 0L) {
    rownames(coarest_set) <- NULL
    return(coarest_set)
  }

  singleton_rows <- lapply(indiv_nodes, function(z) {
    if (nrow(template) > 0L) {
      row <- template[1, , drop = FALSE]
      row$node <- z
      row$depth <- 1L
      row$leaves <- I(list(as.integer(z)))
      row$latent_children <- I(list(NA_integer_))
      row$weight <- 0
      row
    } else {
      data.frame(
        node = z,
        depth = 1L,
        leaves = I(list(as.integer(z))),
        latent_children = I(list(NA_integer_)),
        weight = 0
      )
    }
  })

  singleton_rows <- do.call(rbind, singleton_rows)
  singleton_rows <- singleton_rows[, names(template), drop = FALSE]
  coarest_set <- rbind(coarest_set, singleton_rows)
  rownames(coarest_set) <- NULL
  coarest_set
}

find_coarest <- function(df, result) {
  df <- as.data.frame(df)
  .validate_tree_df(df, require_weight = TRUE)

  required <- c("node", "depth", "leaves", "latent_children", "weight")
  .check_required_columns(result, required, "result")

  p <- find_p(df)
  if (p == 0L) stop("tree_df contains no leaves.", call. = FALSE)
  if (!nrow(result)) return(.add_singleton_rows(result, result, p))

  result <- result[order(result$depth, result$node), , drop = FALSE]
  root_set <- result[result$depth == 0, , drop = FALSE]
  if (!nrow(root_set)) stop("result must contain a root row with depth = 0.", call. = FALSE)

  root_nonzero <- !is.na(root_set$weight) & root_set$weight != 0
  if (any(root_nonzero)) {
    return(.add_singleton_rows(root_set[root_nonzero, , drop = FALSE], result, p))
  }

  coarest_set <- root_set[FALSE, , drop = FALSE]
  max_depth <- max(result$depth, na.rm = TRUE)
  all_root_leaves <- sort(unique(unlist(root_set$leaves, use.names = FALSE)))

  for (current_depth in seq_len(max_depth)) {
    sub_result <- result[result$depth == current_depth, , drop = FALSE]

    if (nrow(coarest_set)) {
      parent_of_sub <- df$parent[match(sub_result$node, df$node)]
      sub_result <- sub_result[!(parent_of_sub %in% coarest_set$node), , drop = FALSE]
    }

    if (!nrow(sub_result)) {
      return(.add_singleton_rows(coarest_set, result, p))
    }

    sub_nonzero <- !is.na(sub_result$weight) & sub_result$weight != 0
    if (any(sub_nonzero)) {
      coarest_set <- rbind(coarest_set, sub_result[sub_nonzero, , drop = FALSE])
    }

    covered <- if (nrow(coarest_set)) {
      sort(unique(unlist(coarest_set$leaves, use.names = FALSE)))
    } else {
      integer(0)
    }
    covered <- covered[!is.na(covered)]

    if (length(covered) == length(all_root_leaves) && all(all_root_leaves %in% covered)) {
      return(.add_singleton_rows(coarest_set, result, p))
    }
  }

  .add_singleton_rows(coarest_set, result, p)
}
