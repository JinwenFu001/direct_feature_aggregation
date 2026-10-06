#!/usr/bin/env Rscript

# Reproduce Table 1 from all 2,200 saved real-data results, without refitting.
# Reference: "final figures and real data (1).ipynb", zero-based cells
# 39, 43, 46, 48, 51, 53, 56, 58, 60, 62 and 64.
# Run from any working directory:
#   Rscript --vanilla "/path/to/project/code/run code/real data study/Table 1.R"
# Requires only R's standard packages (including grid); no package installation.
# Inputs: output/real data study/<analysis>/<analysis>_weight-0.5_uu1..200.RData
# Outputs: manuscript/tables/Table 1.pdf
#          output/real data study/Table 1 summaries.rds
#
# IMPORTANT: The notebook retains vector-valued new.lasso.beta in continuous
# results. as.data.frame() therefore repeats scalar values across multiple rows
# in a replication; apply(..., 2, mean) averages those expanded rows. We preserve
# that behavior exactly, including NA propagation. This is NOT an equal-weight
# mean across the 200 replications. Binary results have one row per replication.
# S-N ratios are not computed; the PDF prints "-" in that column.

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) stop("Run this file with Rscript.", call. = FALSE)
  file_path <- sub("^--file=", "", file_arg[1])
  if (!file.exists(file_path)) file_path <- gsub("~+~", " ", file_path, fixed = TRUE)
  dirname(normalizePath(file_path, mustWork = TRUE))
}

read_analysis <- function(results_root, stem, continuous, notebook_cell) {
  reps <- seq_len(200L)
  filenames <- paste0(stem, "_weight-0.5_uu", reps, ".RData")
  paths <- file.path(results_root, stem, filenames)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing result file(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  beta_positions <- if (continuous) c(1L, 4L, 7L, 10L) else c(1L, 5L, 9L, 13L)
  beta_names <- c("lasso.beta", "ridge.beta", "rare.beta", "our.beta")
  required_metrics <- c("our.loss", "rare.loss", "lasso.loss", "ridge.loss")
  if (!continuous) required_metrics <- c(required_metrics, "our.auc", "rare.auc")
  sizes <- matrix(NA_integer_, nrow = 200L, ncol = 3L,
                  dimnames = list(NULL, c("n", "p", "n_test")))
  expanded_rows <- integer(200L)
  result.all <- vector("list", 200L)

  for (j in reps) {
    env <- new.env(parent = baseenv())
    load(paths[j], envir = env)
    needed <- c("result", "config", "train.index", "test.index")
    if (!all(vapply(needed, exists, logical(1), envir = env, inherits = FALSE))) {
      stop("Missing saved object(s) in ", paths[j], call. = FALSE)
    }
    result <- env$result
    if (!identical(names(result)[beta_positions], beta_names) ||
        !all(required_metrics %in% names(result))) {
      stop("Unexpected result schema in ", paths[j], call. = FALSE)
    }
    if (!isTRUE(as.numeric(env$config$replication) == j) ||
        !isTRUE(as.numeric(env$config$weight.order) == -0.5)) {
      stop("Unexpected replication or weight in ", paths[j], call. = FALSE)
    }
    if (!all(vapply(result[required_metrics], function(x) {
      is.numeric(x) && length(x) == 1L && is.finite(x)
    }, logical(1)))) {
      stop("Invalid Table 1 metric in ", paths[j], call. = FALSE)
    }
    sizes[j, ] <- c(length(env$train.index) + length(env$test.index),
                    length(result$our.beta), length(env$test.index))

    # Exact notebook transformations; do not compress new.lasso.beta or
    # replace the following data.frame/rbind/apply sequence by scalar means.
    if (continuous) {
      result[[1]] = length((result[[1]]))
      result[[4]] = length(table(result[[4]]))
      result[[7]] = length(table(result[[7]]))
      result[[10]] = length(table(result[[10]]))
    } else {
      result[[1]] = length(table(result[[1]]))
      result[[5]] = length(table(result[[5]]))
      result[[9]] = length(table(result[[9]]))
      result[[13]] = length(table(result[[13]]))
    }
    result.all[[j]] = as.data.frame(result)
    expanded_rows[j] <- nrow(result.all[[j]])
  }

  if (any(apply(sizes, 2L, function(x) length(unique(x))) != 1L)) {
    stop("Inconsistent data dimensions across replications for ", stem, call. = FALSE)
  }
  df_merged <- do.call(rbind, result.all)
  notebook_means <- apply(df_merged, 2, mean)
  cat("\n", stem, ": 200 files, ", nrow(df_merged), " notebook rows\n", sep = "")
  print(notebook_means)
  list(
    analysis = stem, notebook_cell = notebook_cell, continuous = continuous,
    n = sizes[1L, "n"], p = sizes[1L, "p"], n_test = sizes[1L, "n_test"],
    replications = reps, weight = -0.5, input_files = file.path(stem, filenames),
    expanded_rows = expanded_rows, notebook_means = notebook_means
  )
}

draw_table_pdf <- function(rows, pdf_file) {
  # One landscape A4 page. Use only standard PDF fonts and grid primitives.
  page_width <- 841.89
  page_height <- 595.28
  grDevices::pdf(pdf_file, width = page_width / 72, height = page_height / 72,
                 family = "Times", title = "Table 1 - Microbiome prediction results",
                 onefile = TRUE, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()

  left <- 32
  widths <- c(170, 62, 90, 138, 79, 79, 79, 79)
  edges <- left + c(0, cumsum(widths))
  centers <- (head(edges, -1L) + tail(edges, -1L)) / 2
  txt <- function(label, x, top, size = 12, face = "plain", just = "centre") {
    grid::grid.text(label, x = grid::unit(x / 72, "inches"),
                    y = grid::unit((page_height - top) / 72, "inches"),
                    just = just,
                    gp = grid::gpar(fontfamily = "Times", fontsize = size,
                                    fontface = face, lineheight = 1.12))
  }
  rule <- function(top, lwd = 0.7) {
    grid::grid.lines(x = grid::unit(c(edges[1L], tail(edges, 1L)) / 72, "inches"),
                     y = grid::unit(rep((page_height - top) / 72, 2), "inches"),
                     gp = grid::gpar(lwd = lwd))
  }
  cell <- function(label, column, top, header = FALSE) {
    size <- if (header) 11.4 else 12
    face <- if (header) "bold" else "plain"
    for (line in strsplit(label, "\n", fixed = TRUE)[[1L]]) {
      grob <- grid::textGrob(line, gp = grid::gpar(fontfamily = "Times",
                             fontsize = size, fontface = face))
      width <- grid::convertWidth(grid::grobWidth(grob), "inches", valueOnly = TRUE) * 72
      if (width > widths[column] - 8) stop("Table cell is too wide: ", label)
    }
    if (column == 1L && !header) txt(label, edges[1L] + 5, top, size, face, "left")
    else txt(label, centers[column], top, size, face)
  }

  txt("Table 1. Prediction results on microbiome datasets", left, 35,
      size = 16, face = "bold", just = "left")
  txt("Continuous and binary outcomes", left, 56, size = 12, just = "left")

  section <- function(title, indices, top, fourth_header) {
    rule(top, 1.0)
    txt(title, mean(range(edges)), top + 14, 13, "bold")
    rule(top + 28)
    labels <- c("Data Source", "Outcome", "Data Size", fourth_header,
                "treeFA\nLoss", "Rare\nLoss", "Lasso\nLoss", "Ridge\nLoss")
    for (j in seq_along(labels)) cell(labels[j], j, top + 47, TRUE)
    rule(top + 66)
    for (i in seq_along(indices)) {
      vals <- as.character(unlist(rows[indices[i], c(
        "source", "outcome", "data_size", "sn_or_auc", "treeFA_loss",
        "RARE_loss", "LASSO_loss", "Ridge_loss")], use.names = FALSE))
      for (j in seq_along(vals)) cell(vals[j], j, top + 80 + (i - 1L) * 25)
    }
    bottom <- top + 66 + length(indices) * 25
    rule(bottom, 1.0)
    bottom
  }
  bottom <- section("Continuous Outcomes", 1:4, 80, "Est. S-N Ratio")
  bottom <- section("Binary Outcomes", 5:11, bottom + 18, "treeFA AUC /\nRare AUC")
}

main <- function() {
  project_root <- Sys.getenv("TREEFA_PROJECT_ROOT",
                             unset = file.path(script_dir(), "..", "..", ".."))
  project_root <- normalizePath(project_root, mustWork = TRUE)
  results_root <- file.path(project_root, "output", "real data study")
  table_dir <- file.path(project_root, "manuscript", "tables")

  spec <- data.frame(
    analysis = c("sinha_2016_continuous", "Erawijantari_2020_continuous",
                 "wang_2020_continuous", "yachida_2019_continuous",
                 "sinha_2016_binary", "kim_2020_binary", "jacobs_2016_binary",
                 "Erawijantari_2020_binary", "franzosa_2019_binary",
                 "wang_2020_binary", "yachida_2019_binary"),
    source = c("Sinha et al. (2016)", "Erawijantari et al. (2020)",
               "Wang et al. (2020)", "Yachida et al. (2019)",
               "Sinha et al. (2016)", "Kim et al. (2020)", "Jacobs et al. (2016)",
               "Erawijantari et al. (2020)", "Franzosa et al. (2019)",
               "Wang et al. (2020)", "Yachida et al. (2019)"),
    outcome = c(rep("BMI", 4), "CRC", "CRC", "IBD", "GC*", "IBD*", "ESRD*", "CRC"),
    continuous = c(rep(TRUE, 4), rep(FALSE, 7)),
    notebook_cell = c(46L, 51L, 56L, 39L, 48L, 60L, 64L, 53L, 62L, 58L, 43L),
    stringsAsFactors = FALSE
  )
  summaries <- lapply(seq_len(nrow(spec)), function(i) {
    read_analysis(results_root, spec$analysis[i], spec$continuous[i], spec$notebook_cell[i])
  })
  names(summaries) <- spec$analysis
  rows <- do.call(rbind, lapply(seq_along(summaries), function(i) {
    s <- summaries[[i]]
    means <- s$notebook_means
    data.frame(
      source = spec$source[i], outcome = spec$outcome[i],
      data_size = paste(s$n, "\u00d7", s$p),
      sn_or_auc = if (s$continuous) "-" else sprintf("%.3f/%.3f", means["our.auc"], means["rare.auc"]),
      treeFA_loss = sprintf("%.3f", means["our.loss"]),
      RARE_loss = sprintf("%.3f", means["rare.loss"]),
      LASSO_loss = sprintf("%.3f", means["lasso.loss"]),
      Ridge_loss = sprintf("%.3f", means["ridge.loss"]),
      stringsAsFactors = FALSE
    )
  }))
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(table_dir, "Table 1.pdf")
  summary_file <- file.path(results_root, "Table 1 summaries.rds")
  draw_table_pdf(rows, pdf_file)
  saveRDS(list(
    notebook = "final figures and real data (1).ipynb", spec = spec,
    summaries = summaries, displayed_table = rows,
    aggregation = "Exact notebook beta transformations, as.data.frame, rbind, apply(..., 2, mean); no replication exclusions or na.rm.",
    sn_ratio = "Not calculated; displayed as -", session_info = sessionInfo()
  ), summary_file)
  cat("\nTable 1 (all 2,200 result files):\n")
  print(rows, row.names = FALSE)
  cat("\nSaved PDF: ", pdf_file, "\nSaved full summaries: ", summary_file, "\n", sep = "")
}

main()
