\name{treeFA-core}
\alias{treeFA-core}
\alias{acc_prox_simple_linear}
\alias{calculate_depth}
\alias{cv.simple_linear}
\alias{detect_non_leaf_nodes}
\alias{find_coarest}
\alias{find_max_param.linear}
\alias{find_p}
\alias{gather_direct_latent_nodes}
\alias{gather_leaf_nodes_per_non_leaf}
\alias{gather_leaves}
\alias{grid.simple_linear}
\alias{make_list}
\alias{real_data_one_round.linear}
\alias{tree_leaf_ids}
\title{Tree-guided feature aggregation core functions}
\description{
Core utilities and model-fitting routines for tree-guided feature aggregation.
The main direct-penalty linear solver is implemented in C++ and called by
\code{acc_prox_simple_linear()}.
}
\usage{
acc_prox_simple_linear(Y, X, tree_result, lambda, warm_start = FALSE,
  init_beta = NULL, intercept = FALSE, ridge_param = 0, thresh = 1e-6,
  max_iter = 10000, algos_cpp = NULL,
  stop_rule = c("objective", "coef"))

grid.simple_linear(Y, X, tree_df, true_beta = NULL, ridge.param = 0,
  seqc = NULL, thresh = 1e-5, intercept = FALSE, max_iter = 10000,
  algos_cpp = NULL, stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy"))

cv.simple_linear(Y, X, tree_df, folds = 5, seqc = NULL, thresh = 1e-5,
  ridge.param = 0, intercept = FALSE, Mmratio = 10000,
  max_iter = 10000, verbose = FALSE, algos_cpp = NULL,
  stop_rule = c("objective", "coef"),
  coarest_rule = c("revised", "legacy"))
}
\arguments{
\item{Y}{Numeric outcome vector.}
\item{X}{Numeric design matrix.}
\item{tree_result}{List produced by \code{gather_leaf_nodes_per_non_leaf(tree_df)$list_result}.}
\item{lambda}{Non-negative direct-penalty tuning parameter.}
\item{warm_start}{Whether to initialize from \code{init_beta}.}
\item{init_beta}{Optional warm-start coefficient vector.}
\item{intercept}{Whether to center \code{Y} and columns of \code{X} and return an intercept.}
\item{ridge_param, ridge.param}{Non-negative ridge parameter.}
\item{thresh}{Stopping threshold.}
\item{max_iter}{Maximum number of proximal-gradient iterations.}
\item{algos_cpp}{Compatibility argument from the script version. The package uses its compiled C++ code and does not need an external \file{algos.cpp}.}
\item{stop_rule}{Stopping rule. \code{"objective"} matches the revised function; \code{"coef"} matches the older \file{scaled_functions.R} coefficient-change rule.}
\item{coarest_rule}{Coarsest-set rule. \code{"revised"} uses the robust package logic; \code{"legacy"} matches the older \file{scaled_functions.R} simulation logic.}
\item{tree_df}{Tree data frame with \code{node}, \code{parent}, and \code{weight}.}
\item{true_beta}{Optional coefficient vector used only to compute path loss.}
\item{seqc}{Optional lambda sequence.}
\item{folds}{Number of cross-validation folds.}
\item{Mmratio}{Ratio used to construct the default cross-validation lambda grid.}
\item{verbose}{Whether to print fold progress.}
}
\value{
Model-fitting functions return coefficient vectors or lists containing fitted
coefficient paths, selected tuning parameters, and validation errors.
}
