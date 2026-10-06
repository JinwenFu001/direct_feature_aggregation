#!/usr/bin/env Rscript

# Reproduce manuscript Figure 4 from the 200 saved experiment2 results.
# Result aggregation, regression, drawing and layout are copied from
# "final figures and real data.ipynb": shared cell 11 and cells 15-18 (zero-based).
# Paths and the requested output filename are adapted to this repository.
#
# From the project root (with dependencies available):
#   Rscript --vanilla "code/run code/simulations/experiment2/Figure 4.R"
# Local configured environment on this machine:
#   .tools/bin/plot-figure4
# Input:  output/experiment2/*.RData
# Output: manuscript/figures/Figure 4.pdf (18 x 6 inches)

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

# Reuse the core tree conversion helper without executing a simulation runner.
source(file.path(code_dir, "other_functions_treeFA.R"), local = environment())

# Tree constructor copied verbatim from run_cluster_experiment2.R.
# For p0 = 1, all within-block distances are zero, so the displayed 20-leaf
# meta-structure is deterministic. No simulation repetitions or models are fitted.
simulate_combined_tree <- function(k, p0) {
  total_leaves <- k * p0
  group_labels <- integer()
  combined_dist <- matrix(0, nrow = total_leaves, ncol = total_leaves)

  start_index <- 1L
  for (i in seq_len(k)) {
    latent_values <- stats::rnorm(p0, mean = i, sd = 0.05)
    subtree_dist <- as.matrix(stats::dist(latent_values))

    end_index <- start_index + p0 - 1L
    combined_dist[start_index:end_index, start_index:end_index] <- subtree_dist
    group_labels <- c(group_labels, rep(i, p0))
    start_index <- end_index + 1L
  }

  combined_dist[combined_dist == 0] <- Inf
  block_indices <- seq(1L, total_leaves, by = p0)

  for (i in seq_len(k - 1L)) {
    start1 <- block_indices[i]
    end1 <- start1 + p0 - 1L
    start2 <- block_indices[i + 1L]
    end2 <- start2 + p0 - 1L

    combined_dist[start1:end1, start2:end2] <- 10
    combined_dist[start2:end2, start1:end1] <- 10
  }

  combined_dist[is.infinite(combined_dist)] <- max(combined_dist[is.finite(combined_dist)]) * 2
  list(tree = stats::hclust(as.dist(combined_dist), method = "average"), group = group_labels)
}

# ---- Notebook cell 11: shared boxplot functions ----
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
# ---- Notebook cell 15: meta-structure ----
visualize_tree_figure3 <- function(tree_df,
                                   highlight_indices = 1:20,   # nodes to colour
                                   circle_size = 8,
                                   font_size   = 1,
                                   main_title  = NULL) {       # optional title
  library(igraph)
  library(RColorBrewer)

  ## ── 1.  graph backbone ───────────────────────────────────────────────
  root_idx <- which(is.na(tree_df$parent))
  edges    <- dplyr::filter(tree_df, !is.na(parent)) %>%
              dplyr::select(parent, node)

  g <- graph_from_data_frame(edges, directed = FALSE, vertices = tree_df)
  vertex_ids <- V(g)$name

  ## ── 2.  vertex fill colours (white vs pastel) ────────────────────────
  pastel         <- brewer.pal(8, "Pastel1")
  vertex_colors  <- rep("white", length(vertex_ids))
  vertex_colors[highlight_indices] <- pastel[1]

  ## ── 3.  dashed/solid edges by parent weight ──────────────────────────
  weight_map <- setNames(tree_df$weight, tree_df$node)
  edge_lty   <- sapply(seq_len(nrow(edges)), function(i) {
                   ifelse(weight_map[edges$parent[i]] == 0, "dashed", "solid")
                 })

  ## ── 4.  leaf labels only ─────────────────────────────────────────────
  leaf_idx    <- which(degree(g) == 1)
  node_labels <- rep("", length(vertex_ids))
  node_labels[leaf_idx] <- vertex_ids[leaf_idx]
  if (font_size == 0) node_labels <- rep("", length(vertex_ids))

  ## ── 5.  draw ─────────────────────────────────────────────────────────
  plot(g,
       layout             = layout_as_tree(g, root = root_idx),
       vertex.size        = circle_size,
       vertex.color       = vertex_colors,
       vertex.frame.color = "black",
       vertex.label       = node_labels,
       vertex.label.color = "black",
       vertex.label.dist  = 1.2,
       vertex.label.degree= pi/2,
       vertex.label.cex   = font_size,
       edge.lty           = edge_lty,
       edge.color         = "black",
       asp                = 0,
       main               = if (!is.null(main_title)) parse(text = main_title),
       cex.main           = if (!is.null(main_title)) 14/12 else 1,
       font.main          = if (!is.null(main_title)) 2 else 1)
}
       
trees = simulate_combined_tree(20, 1)
tree_df = hclust_to_df(trees$tree)
# The preliminary notebook preview is omitted; the tree is drawn in the final panel.
# ---- Notebook cell 16: saved-result reader ----
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(colorspace)
library(RColorBrewer)

## ── 1. 路径与实验设置 ──────────────────────────────────────────────
result_dir <- file.path(project_root, "output", "experiment2")

p0s <- c(3, 5, 10, 20, 30, 40, 50)
uu_use <- 1:200
method_names <- c("treeFA", "RARE", "RS-DL2", "RS-CL2", "oLS", "oRidge")

files <- file.path(result_dir, paste0(
  "Valid_Tree1_n50_p0Incre_k20_sigma1_idealError_weight-1_uu",
  uu_use, "_treeFA_RARE_RS_runtime.RData"
))

missing_uu <- uu_use[!file.exists(files)]
if (length(missing_uu)) {
  stop("缺少结果文件，uu = ", paste(missing_uu, collapse = ", "))
}

## ── 2. 读取并检查：每次重复应包含 7 × 6 条结果 ───────────────────────
long_list <- lapply(seq_along(files), function(i) {
  e <- new.env(parent = baseenv())
  load(files[i], envir = e)

  if (!isTRUE(e$weight.order == -1) || !isTRUE(e$use_weight_order) ||
      !isTRUE(e$sigma == 1) || !isTRUE(e$k == 20) ||
      !isTRUE(e$uu_rep == uu_use[i]) ||
      !identical(as.numeric(e$p0s), p0s)) {
    stop("实验设置不一致：", basename(files[i]))
  }

  d <- e$combined_result
  required <- c(
    "p0", "method", "test_mse", "test_mse_ideal", "rand", "rand_ideal"
  )

  if (!is.data.frame(d) || !all(required %in% names(d))) {
    stop("结果格式不符：", basename(files[i]))
  }

  expected <- expand.grid(p0 = p0s, method = method_names)
  key <- paste(d$p0, d$method, sep = "|")
  expected_key <- paste(expected$p0, expected$method, sep = "|")

  if (anyDuplicated(key) || !setequal(key, expected_key)) {
    stop("存在重复或缺失的 p0/method：", basename(files[i]))
  }

  if (any(!is.finite(d$test_mse)) || any(d$test_mse < 0) ||
      any(!is.finite(d$rand[d$method %in% c("treeFA", "RARE")]))) {
    stop("误差或分组指标包含无效值：", basename(files[i]))
  }

  d$replication <- uu_use[i]
  d$weight.order <- e$weight.order
  d
})

names(long_list) <- paste0("uu", uu_use)

# 完整长表：保留所有方法和误差指标
result.long <- do.call(rbind, long_list)
rownames(result.long) <- NULL

## ── 3. 转为旧绘图代码需要的 result.all 格式 ─────────────────────────
# result.all[[1]][[rep]]：7 行对应 p0s；前 12 列沿用旧列顺序。
# 普通调参后的无噪声测试误差取 test_mse，不取 test_mse_ideal。

result.all <- list(`weight-1` = lapply(long_list, function(d) {

  get_value <- function(method, metric = "test_mse") {
    z <- d[d$method == method, , drop = FALSE]
    as.numeric(z[[metric]][match(p0s, z$p0)])
  }

  out <- data.frame(
    our.minloss       = get_value("treeFA"),                   # 1
    rare.minloss      = get_value("RARE"),                     # 2
    our.rand          = get_value("treeFA", "rand"),           # 3
    rare.rand         = get_value("RARE", "rand"),             # 4
    our.minloss.ideal = get_value("treeFA", "test_mse_ideal"), # 5
    rare.minloss.ideal = get_value("RARE", "test_mse_ideal"),  # 6
    our.rand.ideal    = get_value("treeFA", "rand_ideal"),     # 7
    rare.rand.ideal   = get_value("RARE", "rand_ideal"),       # 8
    oracle.ls.loss    = get_value("oLS"),                      # 9
    ls.loss           = NA_real_,                             # 10：未运行
    ridge.loss        = NA_real_,                             # 11：未运行
    oracle.ridge.loss = get_value("oRidge"),                   # 12
    row.names = as.character(p0s)
  )

  if (any(out$rare.minloss <= 0)) {
    stop("RARE 的测试误差为零，无法计算比值；uu = ", d$replication[1])
  }

  # 同一次重复、同一个 p0 下的误差比
  out$ratio <- out$our.minloss / out$rare.minloss               # 13

  if (any(!is.finite(out$ratio))) {
    stop("误差比包含非有限值；uu = ", d$replication[1])
  }

  # RS 结果追加在后面，不占用 ratio 的第 13 列
  out$rs.dl2.minloss <- get_value("RS-DL2")                    # 14
  out$rs.cl2.minloss <- get_value("RS-CL2")                    # 15
  out$rs.dl2.rand <- get_value("RS-DL2", "rand")               # 16
  out$rs.cl2.rand <- get_value("RS-CL2", "rand")               # 17

  out
}))

## ── 4. 白点：各 p0 下，200 次重复的 ratio 均值 ──────────────────────
# 先计算每次重复的比值，再取均值；不是两种方法平均误差的比值。
ratio_matrix <- do.call(
  cbind,
  lapply(result.all[[1]], function(d) d$ratio)
)
ratio_mean <- rowMeans(ratio_matrix)

obs_df_overlay <- data.frame(
  dimFactor = factor(p0s, levels = p0s),
  pred = ratio_mean
)

## ── 5. 黑线：沿用原来的均值拟合 ────────────────────────────────────
# x 是总维数 p = 20 * p0。
# 这是经验拟合曲线，不是理论误差界。
df_fit <- data.frame(
  x = 20 * p0s,
  y = ratio_mean
)

fit <- lm(y ~ 1 + I(1 / sqrt(log(x))), data = df_fit)

pred_df_overlay <- data.frame(
  dimFactor = factor(p0s, levels = p0s),
  pred = as.numeric(predict(fit, newdata = df_fit))
)

## ── 6. 检查读取情况，并准备图形输出目录 ────────────────────────────
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cat(
  "已读取", length(result.all[[1]]),
  "次重复；weight.order = -1。\n"
)

print(data.frame(
  p0 = p0s,
  p = 20 * p0s,
  mean_ratio = ratio_mean
))

# 接下来直接运行你贴的绘图代码。
# ---- Notebook cell 17: final regression and overlays ----
order=1
for(i in 1:length(result.all[[order]])){
    result.all[[order]][[i]]$ratio=result.all[[order]][[i]]$our.minloss/result.all[[order]][[i]]$rare.minloss
}
ps=as.numeric(rownames(result.all[[order]][[1]]))*20
means=colSums(matrix(unlist(lapply(result.all[[order]],function(x) x$ratio/length(result.all[[order]]))),ncol=7,byrow=T))
ps
means
df_fit <- data.frame(
  x = ps,
  y = means
)

# Fit the model y ~ a/log(x) by using: y ~ 0 + I(1/log(x))
fit <- lm(y ~ 1+ I(1/sqrt(log(x))), data = df_fit)
summary(fit)

# Extract fitted values at each original x
df_fit$pred <- predict(fit, newdata = df_fit)

pred_df_overlay <- data.frame(
  dimFactor = factor(ps/20),
  pred      = df_fit$pred
)

obs_df_overlay <- data.frame(
  dimFactor = factor(ps/20),
  pred      = means
)
# Retain the notebook preview without creating an unrelated Rplots.pdf.
grDevices::pdf(file = NULL)
preview_device <- grDevices::dev.cur()
on.exit({
  if (preview_device %in% grDevices::dev.list()) {
    grDevices::dev.off(preview_device)
  }
}, add = TRUE)

# ---- Notebook cell 18: final panels and PDF export ----
#pdf("figures/fig3b_ratio.pdf")

## ── packages ──────────────────────────────────────────────────────────
library(igraph)
library(ggplot2)
library(ggplotify)      # as.ggplot()
library(patchwork)      # patchwork syntax
library(grid)           # for unit()

## ---------------------------------------------------------------------
##  1.  TREE  (visualize_tree_figure3)  → as ggplot grob
## ---------------------------------------------------------------------
p_tree <- as.ggplot(function() {
  visualize_tree_figure3(
    tree_df,                     # <- combine_tree result
    highlight_indices = trees$group,
    circle_size  = 5,
    font_size    = 0,
    main_title   = NULL           # no base main; ggtitle will be added
  )
}) +
  ggtitle("(a) The Meta-structure")   # ggplot title → themeable

## ---------------------------------------------------------------------
##  2.  OVERLAY PLOT  (“a” + white line / points)
## ---------------------------------------------------------------------
a=boxplotFromList(result.all[[1]],box_width = 0.5,names_lab = c(3,5,10,20,30,40,50),show_legend = TRUE,legend_labels=c("Error Ratio"),legend_title = "(b) Method",
                  measure_cols = c(13),dimension_indices = 1:7,plot_title = "Prediction Error Ratio",y_limits = c(0.2,1.5),x_lab = expression(p[s]),y_lab = "")
p_overlay <- a +                             # your original ggplot object
  geom_line(data = pred_df_overlay,
            aes(x = dimFactor, y = pred, group = 1),
            colour = "black", inherit.aes = FALSE) +
  geom_point(data = obs_df_overlay,
             aes(x = dimFactor, y = pred),
             colour = "white", size = 2, inherit.aes = FALSE) +
  ggtitle("(b) Prediction Error Ratio")    # suppress any legend here

## ---------------------------------------------------------------------
##  3.  BOXPLOT  (boxplotFromList)  – keep its legend
## ---------------------------------------------------------------------
p_box <- boxplotFromList(result.all[[1]],names_lab = c(3,5,10,20,30),legend_labels = c("treeFA","RARE"),legend_title = "(c) Methods",
                         measure_cols      = c(3, 4),
                         dimension_indices = 1:5,
                         plot_title = "(c) Group Accuracy",
                         y_lab      = "",
                         x_lab      = expression(p[s]),
                         col_vec    = c(2, 4),
                         box_width  = 0.5,
                         dodge_width= 0.40,
                         show_legend= TRUE)        # keep legend

font_scale <- 1.2

## 4. 标题样式，并统一放大字体
title_theme <- theme(
  plot.title = element_text(
    hjust = 0.5, face = "plain", size = 18
  ),
  plot.margin = margin(t = 8, r = 5, b = 5, l = 5)
)

scale_fonts <- function(p, factor = 1.2) {
  th <- theme_get() + p$theme

  for (nm in names(th)) {
    el <- th[[nm]]

    # 只放大绝对字号；相对字号通过继承自动放大
    if (inherits(el, c("element_text", "ggplot2::element_text")) &&
        !is.null(el$size) &&
        !inherits(el$size, "rel")) {
      el$size <- el$size * factor
      th[[nm]] <- el
    }
  }

  p + th
}

# 使用新对象，避免重复运行时字号继续累乘
p_tree_large <- scale_fonts(p_tree + title_theme, font_scale)
p_overlay_large <- scale_fonts(p_overlay + title_theme, font_scale)
p_box_large <- scale_fonts(p_box + title_theme, font_scale)

## 5. 拼图
panel <- (p_tree_large | p_overlay_large | p_box_large) +
  plot_layout(
    guides = "collect",
    widths = c(1.4, 1, 1)
  ) &
  theme(
    legend.position = "bottom",
    # 此处会覆盖各面板的 legend 设置，因此也乘以 1.2
    legend.title = element_text(size = 18 * font_scale),
    legend.text = element_text(size = 14 * font_scale),
    legend.key.size = unit(0.9, "cm")
  )

## 6. 显示和保存
print(panel)

ggsave(
  file.path(figure_dir, "Figure 4.pdf"),
  panel,
  width = 18,
  height = 6
)
cat("Figure 4 saved to: ", file.path(figure_dir, "Figure 4.pdf"), "\n", sep = "")
}

main()
