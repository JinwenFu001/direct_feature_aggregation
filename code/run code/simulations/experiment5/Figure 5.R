#!/usr/bin/env Rscript

# Reproduce the complete manuscript Figure 5 from saved experiment3-6 results.
# All aggregation, validation and plotting logic is copied from cell 20
# (zero-based) of "final figures and real data.ipynb"; only paths, filenames
# and the script wrapper are adapted for this repository.
# Identical copies of this complete script are supplied in experiment3-6.
# Run any one copy; it reads all four experiments and produces the same figure.
#
# From the project root (with dependencies available), for example:
#   Rscript --vanilla "code/run code/simulations/experiment3/Figure 5.R"
# Local configured environment on this machine:
#   .tools/bin/plot-figure5
# Input:  output/experiment3/ through output/experiment6/
# Output: manuscript/figures/Figure 5.pdf and Figure 5.png (26 x 12 inches)

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
r_lib <- Sys.getenv("TREEFA_R_LIB", unset = file.path(project_root, "code", "core_code", "rlib"))
.libPaths(c(r_lib, .libPaths()))
required <- c("ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing or unloadable R package(s): ", paste(missing, collapse = ", "),
       ". Install these packages for the current R/platform before plotting.",
       call. = FALSE)
}

# ---- Notebook cell 20: Figure 5 ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ================================================================
## 1. Paths and settings
## ================================================================
results_dir <- file.path(project_root, "output")
font_scale <- 1.2

scenario1_pattern <- paste0(
  "^Valid_BienTree_n50_p100_kIncre_weight-0\\.5_uu[0-9]+_",
  "treeFA_RARE_RS\\.RData$"
)
scenario2_pattern <- paste0(
  "^Valid_BienTree_n500_pIncre_kpRatio025_sigma1_idealError_",
  "weight-0\\.5_uu[0-9]+_treeFA_RARE_RS\\.RData$"
)
scenario3_pattern <- paste0(
  "^Valid_BienTree_Logistic_n50_pIncre_k20_weightChange_",
  "weight-0\\.5_uu[0-9]+\\.RData$"
)
scenario4_pattern <- paste0(
  "^MisspecifiedTree_n50_p200_q10_outer10_inner40_weight-0\\.5_",
  "uu[0-9]+_treeFA_RARE_RS\\.RData$"
)

choose_result_dir <- function(paths, pattern, label) {
  found <- paths[vapply(paths, function(path) {
    dir.exists(path) &&
      length(list.files(path, pattern = pattern)) > 0L
  }, logical(1))]
  if (!length(found)) {
    stop("Cannot find new weight=-0.5 results for ", label,
         ". Checked:\n", paste(paths, collapse = "\n"))
  }
  cat(label, " directory:\n  ", found[1], "\n", sep = "")
  found[1]
}

scenario1_dir <- choose_result_dir(
  file.path(results_dir, "experiment3"),
  scenario1_pattern, "Scenario 1"
)
scenario2_dir <- choose_result_dir(
  file.path(results_dir, "experiment4"),
  scenario2_pattern, "Scenario 2"
)
scenario3_dir <- choose_result_dir(
  file.path(results_dir, "experiment5"),
  scenario3_pattern, "Scenario 3"
)
scenario4_dir <- choose_result_dir(
  file.path(results_dir, "experiment6"),
  scenario4_pattern, "Scenario 4"
)

figure_dir <- file.path(project_root, "manuscript", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
pdf_file <- file.path(figure_dir, "Figure 5.pdf")
png_file <- file.path(figure_dir, "Figure 5.png")

## Retain the four methods and colors already present in this figure.
method_levels <- c("O-LS", "O-Ridge", "treeFA", "RARE")
method_colors <- c(
  "O-LS"    = "#B8D4B1",
  "O-Ridge" = "#C8B4CC",
  "treeFA"  = "#DDA19C",
  "RARE"    = "#A1B6C8"
)
method_map <- c(
  "oLS" = "O-LS", "O-LS" = "O-LS",
  "oRidge" = "O-Ridge", "O-Ridge" = "O-Ridge",
  "treeFA" = "treeFA", "RARE" = "RARE"
)
c_levels <- c(0, 0.10, 0.20, 0.35, 0.55, 0.80, 1.10)

extract_uu <- function(path) {
  as.integer(sub(".*_uu([0-9]+).*", "\\1", basename(path)))
}

result_files <- function(result_dir, pattern, label, expected_reps) {
  files <- list.files(result_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No new result files found for ", label)
  ids <- vapply(files, extract_uu, integer(1))
  if (anyNA(ids) || anyDuplicated(ids) ||
      any(!ids %in% seq_len(expected_reps))) {
    stop(label, ": invalid or duplicate replication IDs.")
  }
  files <- files[order(ids)]
  cat(label, ": loaded ", length(files), " replications\n", sep = "")
  if (length(files) != expected_reps) {
    warning(label, ": expected ", expected_reps, " replications, found ",
            length(files), ". Using the available replications.", call. = FALSE)
  }
  files
}

load_result <- function(file, object_name) {
  env <- new.env(parent = emptyenv())
  load(file, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Missing ", object_name, " in: ", file)
  }
  if (exists("weight.order", envir = env, inherits = FALSE)) {
    weight <- as.numeric(env$weight.order)
    if (length(weight) != 1L || !is.finite(weight) ||
        abs(weight + 0.5) > 1e-10) {
      stop("Unexpected weight.order in: ", file)
    }
  }
  if (exists("uu_rep", envir = env, inherits = FALSE) &&
      !isTRUE(as.numeric(env$uu_rep) == extract_uu(file))) {
    stop("Saved uu_rep disagrees with filename: ", file)
  }
  as.data.frame(env[[object_name]])
}

check_columns <- function(data, required_cols, file) {
  missing <- setdiff(required_cols, names(data))
  if (length(missing)) {
    stop("Missing columns in ", file, ": ", paste(missing, collapse = ", "))
  }
}

as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))

check_panel_rows <- function(data, x_values, methods, file) {
  x_id <- match(round(data$x, 8), round(x_values, 8))
  if (anyNA(x_id) || any(!data$method %in% methods)) {
    stop("Unexpected x values or methods in: ", file)
  }
  data$x <- x_values[x_id]
  counts <- table(factor(data$x, levels = x_values),
                  factor(data$method, levels = methods))
  if (any(counts != 1L)) {
    stop("Missing or duplicate x/method rows in: ", file)
  }
  if (any(!is.finite(data$value))) {
    stop("Non-finite plotted values in: ", file)
  }
  data
}

## ================================================================
## 2. Read new Scenarios 1–3
## ================================================================
## Linear results use combined_result. test_mse is the saved noiseless
## error after tuning on observed Y; no ideal-tuned columns are used.
read_linear_scenario <- function(result_dir, pattern, x_column, x_values,
                                 label, expected_reps = 200L,
                                 group_column = "rand", with_oracle = FALSE) {
  files <- result_files(result_dir, pattern, label, expected_reps)
  loaded <- lapply(files, function(file) {
    result <- load_result(file, "combined_result")
    required <- c(x_column, "method", "test_mse", group_column)
    if (with_oracle) required <- c(required, "oracle_rand_true_group")
    check_columns(result, required, file)
    mapped_method <- unname(method_map[as.character(result$method)])
    keep <- !is.na(mapped_method)
    result <- result[keep, , drop = FALSE]
    result$method <- mapped_method[keep]
    replication <- extract_uu(file)

    prediction <- data.frame(
      replication = rep(replication, nrow(result)),
      x = as_number(result[[x_column]]),
      method = result$method,
      value = as_number(result$test_mse)
    )
    group_rows <- result[result$method %in% c("treeFA", "RARE"), , drop = FALSE]
    grouping <- data.frame(
      replication = rep(replication, nrow(group_rows)),
      x = as_number(group_rows[[x_column]]),
      method = group_rows$method,
      value = as_number(group_rows[[group_column]])
    )
    out <- list(
      prediction = check_panel_rows(prediction, x_values, method_levels, file),
      grouping = check_panel_rows(grouping, x_values, c("treeFA", "RARE"), file)
    )
    if (with_oracle) {
      oracle <- result[result$method == "treeFA", , drop = FALSE]
      out$oracle <- data.frame(
        replication = rep(replication, nrow(oracle)),
        x = as_number(oracle[[x_column]]),
        value = as_number(oracle$oracle_rand_true_group)
      )
      if (any(!is.finite(out$oracle$value))) {
        stop("Non-finite oracle Rand indices in: ", file)
      }
    }
    out
  })
  out <- list(
    prediction = do.call(rbind, lapply(loaded, `[[`, "prediction")),
    grouping = do.call(rbind, lapply(loaded, `[[`, "grouping"))
  )
  if (with_oracle) {
    oracle <- do.call(rbind, lapply(loaded, `[[`, "oracle"))
    out$oracle_points <- stats::aggregate(
      value ~ x, data = oracle,
      FUN = function(z) mean(z, na.rm = TRUE)
    )
  }
  out
}

## Logistic results still save the original wide-format result object.
## Retain its four loss columns and the two ordinary-tuned Rand columns.
read_logistic_scenario <- function(result_dir, pattern, x_values, label) {
  files <- result_files(result_dir, pattern, label, expected_reps = 200L)
  loaded <- lapply(files, function(file) {
    result <- load_result(file, "result")
    check_columns(result, c(
      "oracle.ls.loss", "oracle.ridge.loss", "our.minloss", "rare.minloss",
      "our.rand", "rare.rand"
    ), file)
    x <- as_number(rownames(result))
    replication <- extract_uu(file)
    prediction <- data.frame(
      replication = rep(replication, nrow(result) * 4L),
      x = rep(x, times = 4L),
      method = rep(method_levels, each = nrow(result)),
      value = c(as_number(result$oracle.ls.loss),
                as_number(result$oracle.ridge.loss),
                as_number(result$our.minloss),
                as_number(result$rare.minloss))
    )
    grouping <- data.frame(
      replication = rep(replication, nrow(result) * 2L),
      x = rep(x, times = 2L),
      method = rep(c("treeFA", "RARE"), each = nrow(result)),
      value = c(as_number(result$our.rand), as_number(result$rare.rand))
    )
    list(
      prediction = check_panel_rows(prediction, x_values, method_levels, file),
      grouping = check_panel_rows(grouping, x_values, c("treeFA", "RARE"), file)
    )
  })
  list(
    prediction = do.call(rbind, lapply(loaded, `[[`, "prediction")),
    grouping = do.call(rbind, lapply(loaded, `[[`, "grouping"))
  )
}

s1 <- read_linear_scenario(
  scenario1_dir, scenario1_pattern, x_column = "k",
  x_values = c(10, 20, 30, 40, 50), label = "Scenario 1"
)
s2 <- read_linear_scenario(
  scenario2_dir, scenario2_pattern, x_column = "p",
  x_values = c(50, 100, 200, 400, 600, 800, 1000), label = "Scenario 2"
)
## Display only p = 400, 600, 800, 1000, as in the original figure.
s2$prediction <- s2$prediction[
  s2$prediction$x %in% c(400, 600, 800, 1000), , drop = FALSE
]
s2$grouping <- s2$grouping[
  s2$grouping$x %in% c(400, 600, 800, 1000), , drop = FALSE
]
s3 <- read_logistic_scenario(
  scenario3_dir, scenario3_pattern,
  x_values = c(50, 100, 200, 400, 600, 800, 1000), label = "Scenario 3"
)

## ================================================================
## 3. Read Scenario 4 and oracle Rand index
## ================================================================
## Retain rand_true_group and the original white-point reference.
s4 <- read_linear_scenario(
  scenario4_dir, scenario4_pattern, x_column = "target_c",
  x_values = c_levels, label = "Scenario 4", expected_reps = 400L,
  group_column = "rand_true_group", with_oracle = TRUE
)
cat("\nScenario 4 maximum possible Rand indices:\n")
print(s4$oracle_points)

## ================================================================
## 4. Plotting function
## ================================================================
make_boxplot <- function(
  data,
  x_levels,
  plot_title,
  x_label,
  y_label = NULL
) {
  data$x_factor <- factor(
    data$x,
    levels = x_levels
  )

  data$method <- factor(
    data$method,
    levels = method_levels
  )

  ggplot(
    data,
    aes(
      x = x_factor,
      y = value,
      fill = method
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.78),
      width = 0.68,
      outlier.shape = NA,
      color = "#3B3B3B",
      linewidth = 0.4
    ) +
    scale_fill_manual(
      values = method_colors,
      limits = method_levels,
      breaks = method_levels,
      drop = FALSE
    ) +
    labs(
      x = x_label,
      y = y_label,
      fill = "Methods",
      title = plot_title
    ) +
    theme_bw(base_size = 18) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 7, 6, 7)
    )
}

## ================================================================
## 5. Construct eight panels
## ================================================================
p1_pred <- make_boxplot(
  s1$prediction,
  c(10, 20, 30, 40, 50),
  "(a) Scenario 1: Increasing K",
  "K",
  "Prediction Error"
)

p2_pred <- make_boxplot(
  s2$prediction,
  c(400, 600, 800, 1000),
  "(b) Scenario 2: Fix K/p = 0.25",
  "p"
)

p3_pred <- make_boxplot(
  s3$prediction,
  c(50, 100, 200, 400, 600, 800, 1000),
  "(c) Scenario 3: Binary outcome",
  "p"
)

p4_pred <- make_boxplot(
  s4$prediction,
  c_levels,
  "(d) Scenario 4: Misspecified tree",
  expression(gamma)
)

p1_group <- make_boxplot(
  s1$grouping,
  c(10, 20, 30, 40, 50),
  "(e) Scenario 1: Increasing K",
  "K",
  "Group Accuracy"
)

p2_group <- make_boxplot(
  s2$grouping,
  c(400, 600, 800, 1000),
  "(f) Scenario 2: Fix K/p = 0.25",
  "p"
)

p3_group <- make_boxplot(
  s3$grouping,
  c(50, 100, 200, 400, 600, 800, 1000),
  "(g) Scenario 3: Binary outcome",
  "p"
)

p4_group <- make_boxplot(
  s4$grouping,
  c_levels,
  "(h) Scenario 4: Misspecified tree",
  expression(gamma)
)

## ================================================================
## 6. Correct prediction-error display ranges
## ================================================================
p1_pred <- p1_pred +
  scale_y_continuous(
    breaks = c(0, 4, 8, 12, 16),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  coord_cartesian(ylim = c(0, 16))

p3_pred <- p3_pred +
  scale_y_continuous(
    breaks = c(0.2, 0.4, 0.6, 0.8),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(ylim = c(0.18, 0.82))

p4_pred <- p4_pred +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 6))

## ================================================================
## 7. Add Scenario 4 oracle Rand-index points
## ================================================================
oracle_rand_points <- s4$oracle_points

oracle_rand_points$x_factor <- factor(
  oracle_rand_points$x,
  levels = c_levels
)

p4_group <- p4_group +
  geom_point(
    data = oracle_rand_points,
    aes(
      x = x_factor,
      y = value,
      shape = "Maximum possible Rand index"
    ),
    inherit.aes = FALSE,
    fill = "white",
    color = "black",
    size = 3.2,
    stroke = 0.9
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "Maximum possible Rand index" = 21
    )
  )

## ================================================================
## 8. Keep one method legend and one white-point legend
## ================================================================
p1_pred <- p1_pred +
  guides(
    fill = guide_legend(
      title = "Methods",
      order = 1
    )
  )

p2_pred  <- p2_pred  + guides(fill = "none")
p3_pred  <- p3_pred  + guides(fill = "none")
p4_pred  <- p4_pred  + guides(fill = "none")
p1_group <- p1_group + guides(fill = "none")
p2_group <- p2_group + guides(fill = "none")
p3_group <- p3_group + guides(fill = "none")

p4_group <- p4_group +
  guides(
    fill = "none",
    shape = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        fill = "white",
        color = "black",
        size = 3.2
      )
    )
  )

## ================================================================
## 9. Assemble, display, and save
## ================================================================
font18_theme <- theme(
  text = element_text(size = 22 * font_scale),
  plot.title = element_text(
    size = 22 * font_scale,
    face = "plain",
    hjust = 0.5,
    margin = margin(b = 8)
  ),
  axis.title = element_text(size = 22 * font_scale),
  axis.text = element_text(size = 22 * font_scale),
  legend.title = element_text(
    size = 22 * font_scale,
    face = "bold"
  ),
  legend.text = element_text(size = 22 * font_scale),
  legend.key.size = grid::unit(0.9, "cm"),
  legend.key.width = grid::unit(1.1, "cm"),
  legend.spacing.x = grid::unit(0.35, "cm")
)

joint_figure <- wrap_plots(
  p1_pred,
  p2_pred,
  p3_pred,
  p4_pred,
  p1_group,
  p2_group,
  p3_group,
  p4_group,
  ncol = 4,
  guides = "collect"
) &
  font18_theme &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

options(
  repr.plot.width = 26,
  repr.plot.height = 12
)

# Keep the notebook preview without creating an unrelated Rplots.pdf.
grDevices::pdf(file = NULL)
preview_device <- grDevices::dev.cur()
on.exit({
  if (preview_device %in% grDevices::dev.list()) {
    grDevices::dev.off(preview_device)
  }
}, add = TRUE)
print(joint_figure)

ggsave(
  filename = pdf_file,
  plot = joint_figure,
  width = 26,
  height = 12,
  units = "in"
)

ggsave(
  filename = png_file,
  plot = joint_figure,
  width = 26,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat(
  "\nSaved figures:\n",
  pdf_file, "\n",
  png_file, "\n",
  sep = ""
)

}

main()
