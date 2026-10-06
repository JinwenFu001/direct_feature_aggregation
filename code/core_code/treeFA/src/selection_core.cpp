#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

namespace {

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct SelectionTree {
  std::vector<std::vector<int> > groups;
  std::vector<double> weights;
};

void check_nonnegative_scalar(double value, const char* name) {
  if (!std::isfinite(value) || value < 0.0) {
    Rcpp::stop("%s must be a finite non-negative scalar.", name);
  }
}

SelectionTree read_selection_tree(const Rcpp::List& groups,
                                  const Rcpp::NumericVector& weights,
                                  int p) {
  if (groups.size() != weights.size()) {
    Rcpp::stop("groups and weights must have the same length.");
  }
  SelectionTree tree;
  tree.groups.reserve(groups.size());
  tree.weights.reserve(weights.size());
  std::vector<int> seen(p, -1);
  for (int k = 0; k < groups.size(); ++k) {
    check_nonnegative_scalar(weights[k], "weights");
    SEXP group_sexp = groups[k];
    if (TYPEOF(group_sexp) != INTSXP && TYPEOF(group_sexp) != REALSXP) {
      Rcpp::stop("Each group must be a numeric vector of feature indices.");
    }
    const Rcpp::NumericVector indices = Rcpp::as<Rcpp::NumericVector>(group_sexp);
    std::vector<int> group;
    group.reserve(indices.size());
    for (int j = 0; j < indices.size(); ++j) {
      const double index = indices[j];
      if (!std::isfinite(index) || index < 1.0 || index > p ||
          std::floor(index) != index) {
        Rcpp::stop("Group indices must be integers between 1 and ncol(X) (or length(x)).");
      }
      const int feature = static_cast<int>(index) - 1;
      if (seen[feature] == k) Rcpp::stop("A group cannot contain duplicate feature indices.");
      seen[feature] = k;
      group.push_back(feature);
    }
    tree.groups.push_back(group);
    tree.weights.push_back(weights[k]);
  }
  return tree;
}

// The R interface validates the laminar groups and supplies descendants first.
// Offset leaves an optional, unpenalized intercept outside all groups.
double group_mean(const VectorXd& x, const std::vector<int>& group, int offset) {
  double mean = 0.0;
  const double size = static_cast<double>(group.size());
  for (std::size_t j = 0; j < group.size(); ++j) mean += x[group[j] + offset] / size;
  return mean;
}

double centered_norm(const VectorXd& x, const std::vector<int>& group,
                     int offset, double mean) {
  double norm = 0.0;
  for (std::size_t j = 0; j < group.size(); ++j) {
    norm = std::hypot(norm, x[group[j] + offset] - mean);
  }
  return norm;
}

void selection_tree_prox_inplace(VectorXd& x, const SelectionTree& tree,
                                 double strength, int offset) {
  if (strength == 0.0) return;
  for (std::size_t k = 0; k < tree.groups.size(); ++k) {
    const std::vector<int>& group = tree.groups[k];
    if (group.size() < 2 || tree.weights[k] == 0.0) continue;
    const double mean = group_mean(x, group, offset);
    const double norm = centered_norm(x, group, offset, mean);
    const double threshold = strength * tree.weights[k];
    if (norm == 0.0 || norm <= threshold) {
      for (std::size_t j = 0; j < group.size(); ++j) x[group[j] + offset] = mean;
    } else {
      const double shrinkage = 1.0 - threshold / norm;
      for (std::size_t j = 0; j < group.size(); ++j) {
        const int index = group[j] + offset;
        x[index] = mean + shrinkage * (x[index] - mean);
      }
    }
  }
}

double selection_tree_penalty(const VectorXd& beta, const SelectionTree& tree,
                              int offset) {
  double value = 0.0;
  for (std::size_t k = 0; k < tree.groups.size(); ++k) {
    const std::vector<int>& group = tree.groups[k];
    if (group.size() < 2 || tree.weights[k] == 0.0) continue;
    const double mean = group_mean(beta, group, offset);
    value += tree.weights[k] * centered_norm(beta, group, offset, mean);
  }
  return value;
}

void selection_soft_threshold(const VectorXd& z, VectorXd& u,
                              double threshold, int offset) {
  if (offset) u[0] = z[0];
  for (int j = offset; j < z.size(); ++j) {
    const double value = z[j];
    u[j] = value > threshold ? value - threshold :
      (value < -threshold ? value + threshold : 0.0);
  }
}

double infinity_norm(const VectorXd& x) {
  return x.size() ? x.cwiseAbs().maxCoeff() : 0.0;
}

// All penalties and the ridge term act on the feature coefficients only.
// Matrix-vector products keep the gradient cost O(np) without a p-by-p Gram.
double selection_smooth_gradient(const Map<MatrixXd>& X, const Map<VectorXd>& Y,
                                 const VectorXd& u, int family, int offset,
                                 double ridge_param, VectorXd& gradient) {
  const int n = X.rows();
  const int p = X.cols();
  VectorXd eta = X * u.tail(p);
  if (offset) eta.array() += u[0];
  if (!eta.allFinite()) {
    gradient.setConstant(std::numeric_limits<double>::quiet_NaN());
    return std::numeric_limits<double>::infinity();
  }
  VectorXd score(n);
  double loss = 0.0;
  if (family == 0) {
    score = eta - Y;
    loss = score.squaredNorm() / (2.0 * n);
  } else {
    for (int i = 0; i < n; ++i) {
      const double value = eta[i];
      if (value >= 0.0) {
        const double exp_negative = std::exp(-value);
        score[i] = Y[i] == 1.0 ? -exp_negative / (1.0 + exp_negative) :
          1.0 / (1.0 + exp_negative);
        loss += ((1.0 - Y[i]) * value + std::log1p(exp_negative)) / n;
      } else {
        const double exp_positive = std::exp(value);
        score[i] = Y[i] == 0.0 ? exp_positive / (1.0 + exp_positive) :
          -1.0 / (1.0 + exp_positive);
        loss += (-Y[i] * value + std::log1p(exp_positive)) / n;
      }
    }
  }
  gradient.tail(p).noalias() = X.transpose() * score / static_cast<double>(n);
  if (offset) gradient[0] = score.mean();
  if (ridge_param != 0.0) {
    const double ridge_scale = ridge_param / static_cast<double>(n);
    gradient.tail(p) += ridge_scale * u.tail(p);
    loss += (ridge_scale / 2.0) * u.tail(p).squaredNorm();
  }
  return loss;
}

double selection_objective(double smooth, const VectorXd& u,
                           const SelectionTree& tree, int offset,
                           double lambda, double tau) {
  double value = smooth;
  if (lambda != 0.0) value += lambda * selection_tree_penalty(u, tree, offset);
  if (tau != 0.0) value += tau * u.tail(u.size() - offset).lpNorm<1>();
  return value;
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericVector selection_prox_cpp(Rcpp::NumericVector x, Rcpp::List groups,
                                      Rcpp::NumericVector weights, double strength) {
  check_nonnegative_scalar(strength, "strength");
  VectorXd result = Rcpp::as<VectorXd>(x);
  if (!result.allFinite()) Rcpp::stop("x must contain only finite values.");
  const SelectionTree tree = read_selection_tree(groups, weights, result.size());
  selection_tree_prox_inplace(result, tree, strength, 0);
  if (!result.allFinite()) Rcpp::stop("Non-finite values encountered while evaluating the tree proximal operator.");
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
double selection_penalty_cpp(Rcpp::NumericVector beta, Rcpp::List groups,
                             Rcpp::NumericVector weights) {
  const VectorXd beta_e = Rcpp::as<VectorXd>(beta);
  if (!beta_e.allFinite()) Rcpp::stop("beta must contain only finite values.");
  const SelectionTree tree = read_selection_tree(groups, weights, beta_e.size());
  return selection_tree_penalty(beta_e, tree, 0);
}

// [[Rcpp::export]]
Rcpp::List selection_fit_cpp(
    Rcpp::NumericVector Y, Rcpp::NumericMatrix X, Rcpp::List groups,
    Rcpp::NumericVector weights, double lambda, double tau, int family,
    bool intercept, double ridge_param, double step_size, double abs_tol,
    double rel_tol, int max_iter, Rcpp::NumericVector init_z, bool keep_history) {
  const Map<VectorXd> Y_e = Rcpp::as<Map<VectorXd> >(Y);
  const Map<MatrixXd> X_e = Rcpp::as<Map<MatrixXd> >(X);
  const int n = X_e.rows();
  const int p = X_e.cols();
  const int offset = intercept ? 1 : 0;
  if (n < 1 || p < 1 || Y_e.size() != n) {
    Rcpp::stop("X must have positive dimensions and length(Y) must equal nrow(X).");
  }
  if (!X_e.allFinite() || !Y_e.allFinite()) Rcpp::stop("X and Y must contain only finite values.");
  if (family != 0 && family != 1) Rcpp::stop("family must be 0 (Gaussian) or 1 (binomial).");
  if (family == 1) {
    for (int i = 0; i < n; ++i) {
      if (Y_e[i] != 0.0 && Y_e[i] != 1.0) Rcpp::stop("Binomial Y must contain only 0 and 1.");
    }
  }
  check_nonnegative_scalar(lambda, "lambda");
  check_nonnegative_scalar(tau, "tau");
  check_nonnegative_scalar(ridge_param, "ridge_param");
  check_nonnegative_scalar(abs_tol, "abs_tol");
  check_nonnegative_scalar(rel_tol, "rel_tol");
  if (abs_tol == 0.0 && rel_tol == 0.0) Rcpp::stop("At least one of abs_tol and rel_tol must be positive.");
  if (!std::isfinite(step_size) || step_size <= 0.0) Rcpp::stop("step_size must be finite and positive.");
  if (max_iter < 1 || max_iter == NA_INTEGER) Rcpp::stop("max_iter must be a positive integer.");
  if (init_z.size() != p + offset) Rcpp::stop("init_z must have length ncol(X) + as.integer(intercept).");
  VectorXd z = Rcpp::as<VectorXd>(init_z);
  if (!z.allFinite()) Rcpp::stop("init_z must contain only finite values.");
  const SelectionTree tree = read_selection_tree(groups, weights, p);

  const int dimension = p + offset;
  VectorXd u(dimension), v(dimension), gradient(dimension), reflected(dimension);
  VectorXd difference(dimension), next_z(dimension);
  std::vector<int> history_iter;
  std::vector<double> history_objective, history_residual, history_scaled;
  if (keep_history) {
    const int reserve = std::min(max_iter, 1024);
    history_iter.reserve(reserve);
    history_objective.reserve(reserve);
    history_residual.reserve(reserve);
    history_scaled.reserve(reserve);
  }
  int iter = 0;
  bool converged = false;
  std::string status = "max_iter";
  double residual = std::numeric_limits<double>::infinity();
  double residual_scaled = residual;
  double smooth = std::numeric_limits<double>::infinity();
  double objective = smooth;

  for (iter = 1; iter <= max_iter; ++iter) {
    if (iter == 1 || iter % 100 == 0) Rcpp::checkUserInterrupt();
    residual = std::numeric_limits<double>::infinity();
    residual_scaled = residual;
    selection_soft_threshold(z, u, step_size * tau, offset);
    smooth = selection_smooth_gradient(X_e, Y_e, u, family, offset, ridge_param, gradient);
    bool finite = u.allFinite() && gradient.allFinite() && std::isfinite(smooth);
    double coefficient_scale = 1.0;
    double stationarity_scale = 1.0;
    if (finite) {
      reflected = 2.0 * u - z - step_size * gradient;
      finite = reflected.allFinite();
    }
    if (finite) {
      v = reflected;
      selection_tree_prox_inplace(v, tree, step_size * lambda, offset);
      finite = v.allFinite();
    }
    if (finite) {
      difference = v - u;
      residual = infinity_norm(difference);
      residual_scaled = residual / step_size;
      const double a_norm = infinity_norm(((z - u) / step_size).eval());
      const double b_norm = infinity_norm(((reflected - v) / step_size).eval());
      coefficient_scale = std::max(1.0, std::max(infinity_norm(u), infinity_norm(v)));
      stationarity_scale = std::max(1.0, std::max(infinity_norm(gradient), std::max(a_norm, b_norm)));
      finite = difference.allFinite() && std::isfinite(residual_scaled) &&
        std::isfinite(a_norm) && std::isfinite(b_norm);
    }
    if (keep_history) {
      objective = selection_objective(smooth, u, tree, offset, lambda, tau);
      finite = finite && std::isfinite(objective);
      history_iter.push_back(iter);
      history_objective.push_back(objective);
      history_residual.push_back(residual);
      history_scaled.push_back(residual_scaled);
    }
    if (!finite) {
      status = "nonfinite";
      break;
    }
    if (residual <= abs_tol + rel_tol * coefficient_scale &&
        residual_scaled <= abs_tol + rel_tol * stationarity_scale) {
      converged = true;
      status = "converged";
      break;
    }
    // Return the state that produced u and its diagnostics. A continuation may
    // evaluate this state once again; it never silently starts one step ahead.
    if (iter == max_iter) break;
    next_z = z + difference;
    if (!next_z.allFinite()) {
      status = "nonfinite";
      break;
    }
    z.swap(next_z);
  }

  objective = selection_objective(smooth, u, tree, offset, lambda, tau);
  if (!std::isfinite(objective)) {
    converged = false;
    status = "nonfinite";
  }
  Rcpp::RObject history = R_NilValue;
  if (keep_history) {
    history = Rcpp::DataFrame::create(
      Rcpp::Named("iter") = Rcpp::wrap(history_iter),
      Rcpp::Named("objective") = Rcpp::wrap(history_objective),
      Rcpp::Named("residual") = Rcpp::wrap(history_residual),
      Rcpp::Named("residual_scaled") = Rcpp::wrap(history_scaled));
  }
  return Rcpp::List::create(
    Rcpp::Named("beta") = Rcpp::wrap(VectorXd(u.tail(p))),
    Rcpp::Named("beta0") = offset ? u[0] : 0.0,
    Rcpp::Named("z") = Rcpp::wrap(z),
    Rcpp::Named("iter") = iter,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("status") = status,
    Rcpp::Named("residual") = residual,
    Rcpp::Named("residual_scaled") = residual_scaled,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("history") = history);
}
