#!/usr/bin/env Rscript

# Reproduce the final taxonomic tree figure (manuscript Figure 6).
# Analysis and plotting: cell 68 (zero-based) of
# "final figures and real data (1).ipynb".
# The notebook's calculations, random seeds, grouping, and plot settings
# are retained; setup, paths, progress messages and output names are adapted.
# Legend order is explicit to reproduce the notebook PDF across R versions.
#
# Run from any working directory:
#   Rscript --vanilla "<project>/code/run code/real data study/sinha_2016_binary/Figure 6.R"
# Local configured environment on this machine: .tools/bin/plot-figure6
#
# Inputs (directly under data/):
#   sinha_2016_data.RData  (identical to the notebook's data1.RData)
#   sinha_2016_genera.tsv (original taxonomy columns, in their original order)
# Outputs:
#   manuscript/figures/Figure 6.pdf  (10 x 10 inches)
#   output/real data study/sinha_2016_binary/CRC_inference.csv
#   output/real data study/sinha_2016_binary/CRC_Fusobacterium_contrasts_BH.csv
# This analysis refits the models using data fission; it does not use the
# saved repeated train/test results used for Table 1.

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) {
    stop("Run this script with Rscript as shown above.", call. = FALSE)
  }
  file_path <- sub("^--file=", "", file_arg[1])
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
input_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "output", "real data study", "sinha_2016_binary")
figure_dir <- file.path(project_root, "manuscript", "figures")
r_lib <- Sys.getenv("TREEFA_R_LIB", unset = file.path(code_dir, "rlib"))
.libPaths(c(r_lib, .libPaths()))

# Check inputs before starting the computationally expensive cross-validation.
input_files <- file.path(input_dir, c("sinha_2016_data.RData", "sinha_2016_genera.tsv"))
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"),
       "\nKeep the original taxonomy column order in sinha_2016_genera.tsv.", call. = FALSE)
}
required <- c(
  "glmnet", "car", "readr", "dplyr", "tidyr", "igraph", "tidygraph",
  "ggraph", "patchwork", "ggplot2", "Rcpp", "RcppEigen"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing or unloadable R package(s): ", paste(missing, collapse = ", "),
       ". Install these packages for the current R/platform before running.", call. = FALSE)
}

## 1. Data fission and model selection ---------------------------------
suppressPackageStartupMessages({
  library(glmnet)
  library(car)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(patchwork)
})

# Use the original notebook's function definitions and proximal operator.
# A working C++ toolchain is needed for Rcpp::sourceCpp (once per R session).
Rcpp::sourceCpp(
  file.path(code_dir, "figure6_notebook_algos.cpp"),
  env = environment(), rebuild = FALSE
)
source(file.path(code_dir, "figure6_notebook_helpers.R"), local = environment())

data_env <- new.env()
load(file.path(input_dir, "sinha_2016_data.RData"), envir = data_env)
crc <- data_env$data

X <- as.matrix(crc$X)
y <- crc$y
tree_df <- assign_weight(
  crc$tree_df, p = ncol(X), weight.order = -0.5
)
stopifnot(length(y) == nrow(X), all(y %in% c(0, 1)))

eps <- 0.9
set.seed(123)
z <- rbinom(n = length(y), size = 1, prob = eps)
y1 <- (1 - z) * y + z * (1 - y)
y2 <- y
offset.vec <- ifelse(
  y1,
  log(1 - eps) - log(eps),
  -log(1 - eps) + log(eps)
)

set.seed(5)
# 保留原 notebook 调用：CV 默认 Gaussian，随后拟合 binomial LASSO。
lasso.result <- cv.glmnet(
  x = X, y = y1, nfolds = 5, type.measure = "deviance"
)
lasso.beta <- as.numeric(coef(glmnet(
  x = X, y = y1,
  family = "binomial", intercept = FALSE,
  alpha = 1, lambda = lasso.result$lambda.min
))[-1, 1])
lasso.selected <- which(lasso.beta != 0)

message("Fitting treeFA with the notebook settings (5-fold CV)...")
set.seed(5)
cv.our.result <- cv.logistic(
  y1, X, tree_df = tree_df, folds = 5, intercept = FALSE
)

message("Fitting RARE with the notebook settings (5-fold CV)...")
set.seed(5)
cv.rare.result <- cv.rare.logistic(
  y1, X, tree_df = tree_df, folds = 5, intercept = FALSE
)

# 与原 notebook 一致：按完全相同的估计系数分组。
our.groups <- as.integer(factor(as.numeric(cv.our.result$beta)))
rare.groups <- as.integer(factor(as.numeric(cv.rare.result$beta$beta)))

our.members <- split(seq_len(ncol(X)), our.groups)
rare.members <- split(seq_len(ncol(X)), rare.groups)

# 保证组、系数和标签都按 1, ..., K 排列。
our.members <- our.members[as.character(seq_len(max(our.groups)))]
rare.members <- rare.members[as.character(seq_len(max(rare.groups)))]

if (length(our.members) != 12L || length(rare.members) != 2L) {
  stop("当前拟合未复现正文中的 treeFA 12 组、RARE 2 组，请核对函数版本。")
}

aggregate_X <- function(members) {
  vapply(
    members,
    function(j) rowSums(X[, j, drop = FALSE]),
    numeric(nrow(X))
  )
}

fit_inference <- function(Z) {
  glm(y2 ~ Z + 0, family = binomial(), offset = offset.vec)
}

newfit.our <- fit_inference(aggregate_X(our.members))
newfit.rare <- fit_inference(aggregate_X(rare.members))

stopifnot(length(lasso.selected) > 0L)
lasso.mod <- fit_inference(X[, lasso.selected, drop = FALSE])


## 2. Taxonomy and mapping fitted groups to the plotted tree ------------
genera <- read_tsv(
  file.path(input_dir, "sinha_2016_genera.tsv"),
  show_col_types = FALSE,
  name_repair = "minimal"
)

# 保留原 notebook 的列顺序：去掉第一列 Sample 和最后一列。
OTU_names <- names(genera)[2:(ncol(genera) - 1L)]
stopifnot(length(OTU_names) == ncol(X))

rank_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

tax_df <- as.data.frame(
  t(vapply(OTU_names, function(s) {
    out <- setNames(
      rep(NA_character_, 6),
      c("d", "p", "c", "o", "f", "g")
    )
    for (part in strsplit(s, ";", fixed = TRUE)[[1]]) {
      item <- strsplit(part, "__", fixed = TRUE)[[1]]
      if (length(item) >= 2L &&
          item[1] %in% names(out) &&
          nzchar(item[2])) {
        out[item[1]] <- item[2]
      }
    }
    unname(out)
  }, character(6))),
  stringsAsFactors = FALSE
)
names(tax_df) <- rank_cols

# 原 notebook 中用于绘图的 taxonomy 修正。
tax_df[63, "Order"] <- "Christensenellales_A"
tax_df[c(38, 39), "Order"] <- "Erysipelotrichales_A"
tax_df[39, "Family"] <- "Erysipelotrichaceae_A"
tax_df[10, "Phylum"] <- tax_df[11, "Phylum"]

nodes_all <- tax_df %>%
  pivot_longer(
    everything(), names_to = "layer", values_to = "taxon"
  ) %>%
  filter(!is.na(taxon)) %>%
  distinct(layer, taxon) %>%
  arrange(layer, taxon) %>%
  mutate(node_id = row_number())

edges_df <- tax_df %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(
    all_of(rank_cols), names_to = "layer", values_to = "taxon"
  ) %>%
  filter(!is.na(taxon)) %>%
  arrange(row_id, match(layer, rank_cols)) %>%
  group_by(row_id) %>%
  mutate(next_taxon = lead(taxon), next_layer = lead(layer)) %>%
  filter(!is.na(next_taxon)) %>%
  ungroup() %>%
  left_join(nodes_all, by = c("layer", "taxon")) %>%
  rename(from = node_id) %>%
  left_join(
    nodes_all,
    by = c("next_layer" = "layer", "next_taxon" = "taxon")
  ) %>%
  rename(to = node_id) %>%
  dplyr::select(from, to) %>%
  distinct()

graph_tax <- tbl_graph(
  nodes = nodes_all, edges = edges_df,
  directed = TRUE, node_key = "node_id"
)

leaf_nodes <- vapply(seq_len(nrow(tax_df)), function(i) {
  j <- tail(which(!is.na(tax_df[i, ])), 1)
  hit <- which(
    nodes_all$layer == rank_cols[j] &
      nodes_all$taxon == tax_df[i, j]
  )
  stopifnot(length(hit) == 1L)
  hit
}, integer(1))

leaf_names <- nodes_all$taxon[leaf_nodes]
tree_sets <- gather_leaf_nodes_per_non_leaf(tree_df)$df_result

group_node <- function(members) {
  if (length(members) == 1L) return(leaf_nodes[members])

  hit <- which(vapply(
    tree_sets$leaves, setequal, logical(1), y = members
  ))
  if (!length(hit)) return(NA_integer_)  # 该组不对应单棵子树

  stopifnot(length(hit) == 1L)
  taxon <- tree_df$name[match(tree_sets$node[hit], tree_df$node)]
  id <- which(nodes_all$taxon == taxon)
  stopifnot(length(id) == 1L)
  id
}

our.nodes <- vapply(our.members, group_node, integer(1))
rare.nodes <- vapply(rare.members, group_node, integer(1))

stopifnot(!anyNA(our.nodes), sum(is.na(rare.nodes)) == 1L)

our.names <- nodes_all$taxon[our.nodes]
rare.names <- ifelse(
  is.na(rare.nodes), "Others", nodes_all$taxon[rare.nodes]
)


## 3. Coefficients, Wald intervals and BH contrast tests ----------------
coef_table <- function(fit, labels, method, sizes) {
  s <- coef(summary(fit))

  if (!isTRUE(fit$converged) ||
      nrow(s) != length(labels) ||
      any(!is.finite(s))) {
    stop("请检查 inference fit：", method)
  }

  data.frame(
    Method = method,
    Group = labels,
    N_features = sizes,
    Estimate = s[, 1],
    SE = s[, 2],
    CI_low = s[, 1] - qnorm(0.975) * s[, 2],
    CI_high = s[, 1] + qnorm(0.975) * s[, 2],
    Pvalue = s[, 4],
    row.names = NULL
  )
}

our.table <- coef_table(
  newfit.our, our.names, "treeFA", lengths(our.members)
)
rare.table <- coef_table(
  newfit.rare, rare.names, "RARE", lengths(rare.members)
)
lasso.table <- coef_table(
  lasso.mod, leaf_names[lasso.selected], "LASSO", 1L
)

inference_table <- bind_rows(our.table, rare.table, lasso.table)
print(inference_table, row.names = FALSE)

fuso <- which(tax_df$Genus == "Fusobacterium")
stopifnot(length(fuso) == 1L)

fuso_group <- our.groups[fuso]
stopifnot(identical(
  as.integer(our.members[[fuso_group]]),
  as.integer(fuso)
))

other_groups <- setdiff(seq_along(our.members), fuso_group)

# 与原 notebook 相同的 11 个 Wald Chi-square 检验。
contrast_table <- bind_rows(lapply(other_groups, function(j) {
  L <- matrix(0, nrow = 1, ncol = length(coef(newfit.our)))
  colnames(L) <- names(coef(newfit.our))
  L[1, fuso_group] <- 1
  L[1, j] <- -1

  test <- car::linearHypothesis(
    newfit.our, L, rhs = 0, test = "Chisq"
  )

  data.frame(
    Contrast = paste("Fusobacterium -", our.names[j]),
    Pvalue = test[2, "Pr(>Chisq)"]
  )
}))

contrast_table$Pvalue_BH <- p.adjust(
  contrast_table$Pvalue, method = "BH"
)
print(contrast_table, row.names = FALSE)

cat(
  "\nNonzero Fusobacterium observations:",
  sum(X[, fuso] != 0), "/", nrow(X), "\n"
)
cat(
  "All 11 BH-adjusted p-values < 0.05:",
  all(contrast_table$Pvalue_BH < 0.05), "\n"
)

# 图中的显著性使用单个系数未经多重校正的 p 值，与原图一致。
cat("\nSignificant groups (unadjusted p < 0.05):\n")
print(
  inference_table[inference_table$Pvalue < 0.05, ],
  row.names = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  inference_table,
  file.path(output_dir, "CRC_inference.csv"),
  row.names = FALSE
)
write.csv(
  contrast_table,
  file.path(output_dir, "CRC_Fusobacterium_contrasts_BH.csv"),
  row.names = FALSE
)


## 4. Figure used in the manuscript ------------------------------------
our.colors <- c(
  "#A65628", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#E41A1C", "#A6CEE3", "#7570B3",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"
)
rare.colors <- ifelse(is.na(rare.nodes), "#8DA0CB", "#E41A1C")

make_panel <- function(anchors, tab, colors, title, legend_title) {
  node_group <- rep("empty", nrow(nodes_all))
  label <- rep("", nrow(nodes_all))
  nonclade <- which(is.na(anchors))

  if (length(nonclade)) node_group[] <- tab$Group[nonclade]

  for (j in which(!is.na(anchors))) {
    descendants <- as.integer(igraph::subcomponent(
      graph_tax, v = anchors[j], mode = "out"
    ))
    node_group[descendants] <- tab$Group[j]
    label[anchors[j]] <- as.character(round(tab$Estimate[j], 2))
  }

  if (length(nonclade)) {
    # 保留原图 Others 系数的标注位置；132 不是该组的定义。
    stopifnot(node_group[132] == tab$Group[nonclade])
    label[132] <- as.character(round(tab$Estimate[nonclade], 2))
  }

  pvalue <- tab$Pvalue[match(node_group, tab$Group)]
  significance <- ifelse(
    !is.na(pvalue) & pvalue < 0.05, "sig", "nsig"
  )

  edge_type <- if (length(nonclade)) {
    ifelse(
      node_group[edges_df$from] == node_group[edges_df$to],
      "dashed", "solid"
    )
  } else {
    ifelse(node_group[edges_df$from] != "empty", "dashed", "solid")
  }

  g <- graph_tax %>%
    activate(nodes) %>%
    mutate(
      subtree_group = node_group,
      signif_group = significance,
      effect_label = label
    ) %>%
    activate(edges) %>%
    mutate(edge_type = edge_type)

  palette <- setNames(colors, tab$Group)
  legend_labels <- setNames(tab$Group, tab$Group)

  if (!length(nonclade)) {
    palette <- c(palette, empty = "grey50")
    legend_labels <- c(legend_labels, empty = "Not aggregated")
  }

  ggraph(g, layout = "tree") +
    geom_edge_bend(aes(linetype = edge_type), colour = "grey40") +
    scale_edge_linetype_manual(
      name = "Edge Type",
      values = c(solid = "solid", dashed = "dashed"),
      labels = c(solid = "not aggregated", dashed = "aggregated"),
      # Fix the notebook PDF's legend order across R/package versions.
      guide = guide_legend(order = 2)
    ) +
    geom_node_point(
      aes(colour = subtree_group, shape = signif_group), size = 2
    ) +
    geom_node_text(
      aes(label = effect_label),
      hjust = -0.3, vjust = 0.04, size = 3
    ) +
    scale_colour_manual(
      name = legend_title, values = palette, labels = legend_labels,
      guide = guide_legend(order = if (length(nonclade)) 4 else 3)
    ) +
    scale_shape_manual(
      name = "Significance (<0.05)",
      values = c(nsig = 16, sig = 2),
      labels = c(nsig = "not significant", sig = "significant"),
      drop = FALSE,
      guide = guide_legend(order = 1)
    ) +
    scale_y_reverse(
      limits = c(5, 0), breaks = 5:0, labels = rank_cols
    ) +
    coord_flip(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = grid::unit(0.6, "cm"),
      legend.spacing = grid::unit(0.5, "cm"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "left",
      plot.margin = margin(t = 10, r = 40, b = 10, l = 40)
    )
}

p_our <- make_panel(
  our.nodes, our.table, our.colors,
  "(a) treeFA Aggregation", "treeFA Groups"
)
p_rare <- make_panel(
  rare.nodes, rare.table, rare.colors,
  "(b) RARE Aggregation", "RARE Groups"
)

combined_plot <- (p_our | p_rare) +
  plot_layout(guides = "collect", widths = c(3, 3), heights = 0.5)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(figure_dir, "Figure 6.pdf"),
  combined_plot, width = 10, height = 10, dpi = 300
)
message("Saved figure: ", file.path(figure_dir, "Figure 6.pdf"))
message("Saved inference and contrast tables: ", output_dir)
}

main()
