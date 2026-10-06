#!/usr/bin/env Rscript

# Reproduce manuscript Figure 3 from the 200 saved experiment1 results.
# Plotting and result-aggregation code is copied from the Figure 3 section of
# "final figures and real data.ipynb" (cells 10-13, zero-based).
# The preliminary simulated-tree example in cell 10 is unnecessary: cell 12
# restores the actual tree and true groups from the saved results.
#
# From the project root:
#   Rscript --vanilla "code/run code/simulations/experiment1/plot_figure3.R"
# Input:  output/experiment1/*.RData
# Output: manuscript/figures/fig2joint.pdf
# The PDF basename and 18 x 6 inch dimensions are kept from the notebook.

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) {
    stop("Run this script with Rscript as shown above.", call. = FALSE)
  }
  file_path <- sub("^--file=", "", file_arg[1])
  # Rscript may encode spaces as ~+~ in commandArgs(), even for quoted paths.
  if (!file.exists(file_path)) file_path <- gsub("~+~", " ", file_path, fixed = TRUE)
  dirname(normalizePath(file_path, mustWork = TRUE))
}

main <- function() {
project_root <- Sys.getenv(
  "TREEFA_PROJECT_ROOT",
  unset = file.path(script_dir(), "..", "..", "..", "..")
)
project_root <- normalizePath(project_root, mustWork = TRUE)
code_dir <- file.path(project_root, "code", "core_code")
figure_dir <- file.path(project_root, "manuscript", "figures")
r_lib <- Sys.getenv("TREEFA_R_LIB", unset = file.path(code_dir, "rlib"))
.libPaths(c(r_lib, .libPaths()))

required <- c(
  "treeFA", "dplyr", "tidyr", "purrr", "tibble", "RColorBrewer",
  "colorspace", "igraph", "ggplot2", "ggplotify", "patchwork"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing or unloadable R package(s): ", paste(missing, collapse = ", "),
    ". Install these packages for the current R/platform before plotting.",
    call. = FALSE
  )
}
if (!"outliers" %in% names(formals(ggplot2::geom_boxplot))) {
  stop("The notebook requires a ggplot2 version with geom_boxplot(outliers = FALSE).",
       call. = FALSE)
}

# This existing helper identifies the same latent nodes used by the notebook.
# Its treeFA dependency only traverses the saved tree; no estimator is fitted.
source(file.path(code_dir, "other_functions_treeFA.R"), local = environment())

# ---- Notebook cell 10: tree plotting function only ----
visualize_tree_figure2 <- function(tree_df, true_group_indices, circle_size = 8, font_size = 1,main_title  = "full tree T[0]") {
  library(igraph)
  
  p = length(true_group_indices)
  
  # Identify which nodes should have the "deeper" color (previously green)
  latent_nodes = unlist(find_latent_nodes_with_exact_leaf_groups(tree_df, true_group_indices))
  
  # 1) Identify root (for layout)
  root_idx <- which(is.na(tree_df$parent))
  
  # 2) Build edge list (parent-child pairs), ignoring NA parent
  edges <- tree_df %>%
    dplyr::filter(!is.na(parent)) %>%
    dplyr::select(parent, node)
  
  # 3) Create an undirected graph from these edges
  g <- graph_from_data_frame(
    d        = edges,
    directed = FALSE,
    vertices = tree_df
  )
  
  # 4) Convert igraph vertex names (often strings) to numeric indices
  #    so we can manage colors, etc.
  vertex_ids <- (V(g)$name)
  
  # 5) Assign black-and-white fill colors
  #    - "white" for normal nodes
  #    - "gray40" (dark gray) for 'latent_nodes'
  #      that were previously "green".
  pastel <- brewer.pal(8, "Pastel1")
  vertex_colors <- rep("white", length(vertex_ids))
  # idx_latent <- which(vertex_ids %in% latent_nodes)
  vertex_colors[latent_nodes] <- pastel[1]
  
  # 6) Determine dashed edges if parent's weight == 0
  weight_by_node <- setNames(tree_df$weight, tree_df$node)
  edge_lty <- sapply(seq_len(nrow(edges)), function(i) {
    parent_id <- edges$parent[i]
    if (!is.na(parent_id) && weight_by_node[parent_id] == 0) {
      return("dashed")
    } else {
      return("solid")
    }
  })
  
  # 7) Only label leaf nodes. A leaf is any vertex with degree == 1
  leaf_idx <- which(degree(g) == 1)
  # Create a label vector: blank for non-leaf, numeric ID for leaf
  node_labels <- rep("", length(vertex_ids))
  node_labels[leaf_idx] <- vertex_ids[leaf_idx]
    
  if (font_size == 0) {
      node_labels <- rep("", length(vertex_ids))
    }
  
  # 8) Plot
  plot(
    g,
    layout             = layout_as_tree(g, root = root_idx),
    
    # (A) Node appearance
    vertex.size        = circle_size,
    vertex.color       = vertex_colors,
    vertex.frame.color = "black",  # black boundary
    
    # (B) Labels: only on leaves, placed below the circle
    vertex.label       = node_labels,
    vertex.label.color = "black",
    vertex.label.dist  = 1.2,      # move label below
    vertex.label.degree= pi/2,     # place underneath (adjust if needed)
    vertex.label.cex   = font_size,
    
    # (C) Edges
    edge.lty           = edge_lty,
    edge.color         = "black",
    
    # (D) Other
    asp                = 0.0,
    #main               = expression("Full tree"~T[0]),  # LaTeX-style T₀
    #cex.main           = 14 / 12,   # ggplot’s 14 pt ≈ 1.17× default
    #font.main          = 2,         # 2 = bold
  )
}

# ---- Notebook cell 11: boxplot functions ----
saturate <- function(col, p = 0.4) {
  hcl <- coords(as(hex2RGB(col), "polarLUV"))
  hcl[, "C"] <- pmin(100, hcl[, "C"] * (1 + p))   # raise chroma
  hex(polarLUV(hcl))
}


boxplotFromList <- function(df_list, 
                            measure_cols = NULL,      
                            dimension_indices = NULL, 
                            y_limits = NULL,
                            plot_title = NULL,
                            x_lab=NULL,
                            y_lab=NULL,
                            names_lab=parse(text = paste0("T[", 0:5, "]")),
                            legend_labels      = NULL,   # <- NEW
                            show_legend        = TRUE,
                            col_vec=(1:length(measure_cols))+1,
                            box_width=0.75,dodge_width=0.8,legend_title="Methods") {

    
  is_digits_only <- function(x) grepl("^[0-9]+$", x)


  
  long_list <- map(df_list, 
    ~ {
      df_temp <- .x
      if (!is.null(dimension_indices)) {
        df_temp <- df_temp[dimension_indices, , drop = FALSE]
      }
      
      if (!is.null(measure_cols)) {
        df_temp <- df_temp[, measure_cols, drop = FALSE]
      }
      
      df_temp %>%
        rownames_to_column("dimension") %>%
        # Convert dimension to numeric, then make it a factor with sorted levels
        mutate(
          dimension_numeric=ifelse(!is_digits_only(dimension),1:length(dimension),as.numeric(dimension)),
          #if(!is_digits_only(dimension[1])) dimension_numeric=1:length(dimension),
          #if(is_digits_only(dimension[1])) dimension_numeric = as.numeric((dimension)),
          
          dimension = factor(dimension_numeric,
                             levels = sort(unique(dimension_numeric)))
        ) %>%
        # We keep 'dimension' as the factor. 'dimension_numeric' is optional.
        pivot_longer(
          cols = c(-dimension, -dimension_numeric),  # pivot everything except dimension columns
          names_to = "measure",
          values_to = "value"
        )
        
        #print((df_temp))
    }
  )
    #print(long_list[[1]])
  
  long_all <- bind_rows(long_list, .id = "dfID")
    
    #if(is.null(names_lab)) legend_names=
  pastel <- brewer.pal(9, "Pastel1")
  mellow <- darken(saturate(pastel, 0.10), 0.10)
  if(length(col_vec)==2) col_vec=mellow[1:length(col_vec)]
  else col_vec=mellow[c(3,4,1,2)]
  
  p <- ggplot(long_all, aes(
    x    = dimension,  # dimension is a factor with numeric-sorted levels
    y    = value,
    fill = measure
  )) +
    geom_boxplot(position = position_dodge(width = dodge_width),outliers = F,width = box_width) +
    scale_fill_manual(
    #values = c("grey90","grey60","grey30","white")
    #values = c("grey30","white")
    values=col_vec,
        labels = legend_labels,
        drop=FALSE
  )+
    scale_x_discrete(labels = names_lab) +
    labs(
      x     = x_lab,
      y     = y_lab,
      fill  = legend_title,
      title = plot_title
    ) +
    theme_bw(base_size = 14)+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 14),axis.title.x = element_text(size = 18))
    
    if (!show_legend)                      # legend switch
    p <- p + guides(fill = "none")
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}
# ---- Notebook cell 12: saved-result reader ----
## 补齐绘图函数需要的 packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(RColorBrewer)
  library(colorspace)
})

## ── 1. 文件位置：weight = -0.5，200 次重复 ────────────────────────
result_dir <- file.path(project_root, "output", "experiment1")

uu <- 1:200

files <- file.path(
  result_dir,
  sprintf(
    "Valid_BienTree_n50_p60_k10_treeChange_weight-0.5_uu%d_treeFA_RARE_RS_modifiedTree.RData",
    uu
  )
)

## 文件未收齐时停止，避免不知情地少用部分重复
missing_uu <- uu[!file.exists(files)]
if (length(missing_uu)) {
  stop("以下重复的结果文件尚不存在：", paste(missing_uu, collapse = ", "))
}

## ── 2. 固定树和方法的顺序 ────────────────────────────────────────
variant_order <- c(
  "tree_df0", "tree1_l0", "tree2_l0",
  "tree3_l0", "tree1_l1", "tree2_l1"
)

method_order <- c(
  "treeFA", "RARE", "RS-DL2", "RS-CL2", "oLS", "oRidge"
)

required_keys <- unlist(
  lapply(method_order, function(m) paste(variant_order, m, sep = "::")),
  use.names = FALSE
)

required_cols <- c(
  "variant", "method", "test_mse", "test_mse_ideal",
  "rand", "rand_ideal"
)

## result.all[[1]][[j]]：第 j 次重复，6 行对应 T0,...,T5
result.all <- list(
  "weight-0.5" = setNames(vector("list", length(uu)), paste0("uu", uu))
)
long_list <- vector("list", length(uu))

## ── 3. 逐文件整理，保持旧绘图代码所用的列编号 ────────────────────
for (i in seq_along(files)) {
  e <- new.env(parent = baseenv())
  load(files[i], envir = e)

  tab <- if (exists("combined_result", envir = e, inherits = FALSE)) {
    e$combined_result
  } else {
    e$result
  }

  if (!all(required_cols %in% names(tab))) {
    stop("结果字段不完整：", basename(files[i]))
  }

  keys <- paste(tab$variant, tab$method, sep = "::")

  if (anyDuplicated(keys)) {
    stop("存在重复的 variant–method 组合：", basename(files[i]))
  }

  missing_keys <- setdiff(required_keys, keys)
  if (length(missing_keys)) {
    stop(
      basename(files[i]), " 缺少：",
      paste(missing_keys, collapse = ", ")
    )
  }

  if (any(!is.finite(tab$test_mse[match(required_keys, keys)]))) {
    stop("test_mse 含 NA/Inf，请检查：", basename(files[i]))
  }

  ## 根据名称匹配，避免依赖结果文件中的行顺序
  get_metric <- function(method, metric) {
    idx <- match(paste(variant_order, method, sep = "::"), keys)
    as.numeric(tab[[metric]][idx])
  }

  result.all[[1]][[i]] <- data.frame(
    our.minloss       = get_metric("treeFA", "test_mse"),        #  1
    rare.minloss      = get_metric("RARE",   "test_mse"),        #  2
    our.rand          = get_metric("treeFA", "rand"),            #  3
    rare.rand         = get_metric("RARE",   "rand"),            #  4
    our.minloss.ideal = get_metric("treeFA", "test_mse_ideal"),  #  5
    rare.minloss.ideal= get_metric("RARE",   "test_mse_ideal"),  #  6
    our.rand.ideal    = get_metric("treeFA", "rand_ideal"),      #  7
    rare.rand.ideal   = get_metric("RARE",   "rand_ideal"),      #  8
    oracle.ls.loss    = get_metric("oLS",    "test_mse"),        #  9
    ls.loss          = rep(NA_real_, 6),                       # 10：未运行
    ridge.loss       = rep(NA_real_, 6),                       # 11：未运行
    oracle.ridge.loss = get_metric("oRidge", "test_mse"),        # 12
    rs.dl2.minloss    = get_metric("RS-DL2", "test_mse"),        # 13
    rs.cl2.minloss    = get_metric("RS-CL2", "test_mse"),        # 14
    rs.dl2.rand       = get_metric("RS-DL2", "rand"),            # 15
    rs.cl2.rand       = get_metric("RS-CL2", "rand"),            # 16
    row.names = as.character(seq_along(variant_order))
  )

  ## 另保留完整长表，包括调参值、验证误差等
  tab$uu <- uu[i]
  long_list[[i]] <- tab

  ## 从保存的树和真实分组矩阵恢复树图所需对象
  A <- as.matrix(e$data$A)

  if (i == 1L) {
    stopifnot(
      all(A %in% c(0, 1)),
      all(rowSums(A) == 1)
    )

    tree_df0 <- e$df_list[[1]]
    A_reference <- A
    trees <- list(group = max.col(A, ties.method = "first"))
  } else {
    if (!isTRUE(all.equal(A, A_reference)) ||
        !isTRUE(all.equal(e$df_list[[1]], tree_df0))) {
      stop("真实分组或原始树与其他重复不一致：", basename(files[i]))
    }
  }
}

result.long <- do.call(rbind, long_list)
rownames(result.long) <- NULL

## ── 4. 检查并准备绘图目录 ────────────────────────────────────────
stopifnot(
  length(result.all[[1]]) == length(uu),
  all(vapply(result.all[[1]], nrow, integer(1)) == 6L)
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cat("已整理", length(uu), "次重复；每次重复包含 6 棵树。\n")
print(result.all[[1]][[1]][, c(9, 12, 1, 2, 3, 4)])

## 现在可直接运行你给出的绘图代码：
## Prediction Error：measure_cols = c(9, 12, 1, 2)
## Group Accuracy：  measure_cols = c(3, 4)
##
## 上述预测误差列均来自 test_mse：普通调参 + 无噪声测试误差。
## RS 结果保留在第 13:16 列；你当前的绘图调用不会画出 RS。
# Retain the notebook's preview print without creating a stray Rplots.pdf.
# The final export below still uses the notebook's original PDF device settings.
grDevices::pdf(file = NULL)
preview_device <- grDevices::dev.cur()
on.exit({
  if (preview_device %in% grDevices::dev.list()) {
    grDevices::dev.off(preview_device)
  }
}, add = TRUE)

# ---- Notebook cell 13: panels and PDF export ----
## ── packages ────────────────────────────────────────────────────────
library(igraph)
library(ggplot2)
library(ggplotify)
library(patchwork)

## ── 1. tree：内部标注字号保持不变 ─────────────────────────────────────
p_tree <- as.ggplot(function() {
  visualize_tree_figure2(
    tree_df0, trees$group,
    circle_size = 5,
    font_size   = 0,
    main_title  = NULL
  )
}) +
  ggtitle(expression("(a) Full Tree" ~ T[0]))

## ── 2. two boxplots ─────────────────────────────────────────────────
p_pred <- boxplotFromList(
  result.all[[1]],
  legend_labels    = c("O-LS", "O-Ridge", "treeFA", "RARE"),
  box_width        = 0.7,
  measure_cols     = c(9, 12, 1, 2),
  dimension_indices = 1:6,
  x_lab            = "Trees",
  y_lab            = "",
  plot_title       = "(b) Prediction Error",
  col_vec          = c(2, 4, 5, 3)
) +
  theme(legend.position = "bottom")

p_group <- boxplotFromList(
  result.all[[1]],
  box_width        = 0.35,
  dodge_width      = 0.4,
  measure_cols     = c(3, 4),
  dimension_indices = 1:6,
  x_lab            = "Trees",
  y_lab            = "",
  plot_title       = "(c) Group Accuracy",
  col_vec          = c(2, 4),
  show_legend      = FALSE
)

title_theme <- theme(
  plot.title = element_text(
    hjust = 0.5,
    face  = "plain",
    size  = 18
  ),
  plot.margin = margin(t = 8, r = 5, b = 5, l = 5)
)

p_tree  <- p_tree  + title_theme
p_pred  <- p_pred  + title_theme
p_group <- p_group + title_theme

## ── 3. 所有 theme 字体增大为原来的 1.5 倍 ─────────────────────────────
## 包括标题、坐标轴和 legend；不影响 tree 内部标注
scale_fonts <- function(p, factor = 1.2) {
  th <- theme_get() + p$theme

  for (nm in names(th)) {
    el <- th[[nm]]

    if (inherits(el, c("element_text", "ggplot2::element_text")) &&
        !is.null(el$size) &&
        !inherits(el$size, "rel")) {
      el$size <- el$size * factor
      th[[nm]] <- el
    }
  }

  p + th
}

p_tree  <- scale_fonts(p_tree)
p_pred  <- scale_fonts(p_pred)
p_group <- scale_fonts(p_group)

## ── 4. assemble panels and collect legend ───────────────────────────
panel <- (p_tree | p_pred | p_group) +
  plot_layout(
    guides = "collect",
    widths = c(1.4, 1, 1)
  ) &
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 18 * 1.2),
    legend.key.size = grid::unit(1, "cm")
  )

print(panel)

## ── 5. save ─────────────────────────────────────────────────────────
pdf(file.path(figure_dir, "fig2joint.pdf"), width = 18, height = 6)
print(panel)
dev.off()
cat("Figure 3 saved to: ", file.path(figure_dir, "fig2joint.pdf"), "\n", sep = "")
}

main()
