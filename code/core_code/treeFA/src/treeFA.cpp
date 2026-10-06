#include <RcppEigen.h>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::VectorXd;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// [[Rcpp::depends(RcppEigen)]]

VectorXd one_layer_cpp(const VectorXd& eta, const List& groups, const NumericVector& weights, const double lambda) {
  VectorXd beta = eta;
  const int n_groups = weights.size();

  for (int i = 0; i < n_groups; ++i) {
    const double weight = weights[i];
    const IntegerVector ind = groups[i];
    const int group_size = ind.size();
    if (group_size == 0) continue;

    double mean_mu = 0.0;
    for (int j = 0; j < group_size; ++j) {
      mean_mu += eta[ind[j] - 1];
    }
    mean_mu /= group_size;

    double ss = 0.0;
    for (int j = 0; j < group_size; ++j) {
      const double diff = eta[ind[j] - 1] - mean_mu;
      ss += diff * diff;
    }

    const double temp = std::sqrt(ss);
    double d = 100.0;
    if (temp > 0.0) d = weight * lambda / temp;

    if (d >= 1.0) {
      for (int j = 0; j < group_size; ++j) {
        beta[ind[j] - 1] = mean_mu;
      }
    } else {
      for (int j = 0; j < group_size; ++j) {
        const int idx = ind[j] - 1;
        beta[idx] = (1.0 - d) * eta[idx] + d * mean_mu;
      }
    }
  }

  return beta;
}

VectorXd prox_tree_eigen(const VectorXd& eta, const double lambda, const List& tree_list) {
  const List groups = tree_list["groups"];
  const List weights = tree_list["weights"];
  const int depth = groups.size();
  VectorXd new_eta = eta;

  for (int i = 0; i < depth; ++i) {
    new_eta = one_layer_cpp(new_eta, groups[i], weights[i], lambda);
  }

  return new_eta;
}

// [[Rcpp::export]]
NumericVector prox_tree_cpp(NumericVector eta, double lambda, List tree_list) {
  const Map<VectorXd> eta_e = Rcpp::as<Map<VectorXd> >(eta);
  const VectorXd out = prox_tree_eigen(eta_e, lambda, tree_list);
  return Rcpp::wrap(out);
}

// [[Rcpp::export]]
List acc_prox_simple_linear_cpp(
    NumericVector Y,
    NumericMatrix X,
    List tree_result,
    double lambda,
    NumericVector init_beta,
    double ridge_param,
    double thresh,
    int max_iter,
    int stop_rule) {
  const Map<VectorXd> Y_e = Rcpp::as<Map<VectorXd> >(Y);
  const Map<MatrixXd> X_e = Rcpp::as<Map<MatrixXd> >(X);
  const int n = X_e.rows();
  const int p = X_e.cols();

  if (init_beta.size() != p) {
    Rcpp::stop("init_beta must have length ncol(X).");
  }

  VectorXd beta0 = Rcpp::as<Map<VectorXd> >(init_beta);
  VectorXd beta1 = beta0;

  double alpha0 = 1.0;
  double alpha1 = 0.5;

  const MatrixXd matXX = X_e.transpose() * X_e;
  const VectorXd matXY = X_e.transpose() * Y_e;

  SelfAdjointEigenSolver<MatrixXd> eig(matXX, Eigen::EigenvaluesOnly);
  double L0 = eig.eigenvalues().maxCoeff() / static_cast<double>(n);
  if (!std::isfinite(L0) || L0 <= 0.0) L0 = std::numeric_limits<double>::epsilon();
  double tao = static_cast<double>(n) / L0;

  int iter = 0;
  int consecutive_below_thresh = 0;

  while (consecutive_below_thresh < 5 && iter < max_iter) {
    ++iter;
    bool accept = false;
    const VectorXd Gam = beta0 + ((alpha0 - 1.0) / alpha1) * (beta1 - beta0);

    const VectorXd matXXGam = matXX * Gam;
    const VectorXd matXGam = X_e * Gam;
    const VectorXd base1 = -matXY + matXXGam + ridge_param * Gam;

    VectorXd eta_new = Gam;
    while (!accept) {
      const VectorXd eta = Gam - tao * base1 / static_cast<double>(n);
      eta_new = prox_tree_eigen(eta, tao * lambda, tree_result);

      const double lhs = (Y_e - X_e * eta_new).squaredNorm() / (2.0 * n);
      const VectorXd delta = eta_new - Gam;
      const double rhs = (Y_e - matXGam).squaredNorm() / (2.0 * n) +
        (matXXGam - matXY).dot(delta) +
        delta.squaredNorm() / (2.0 * tao);

      if (lhs <= rhs || tao <= 1.0 / L0) {
        accept = true;
      } else {
        tao = std::max(tao / 2.0, 1.0 / L0);
      }
    }

    beta0 = beta1;
    beta1 = eta_new;

    double dis;
    if (stop_rule == 1) {
      const double denom = std::max(beta0.squaredNorm(), std::numeric_limits<double>::epsilon());
      dis = std::sqrt((beta1 - beta0).squaredNorm() / denom);
    } else {
      const double old_obj = (Y_e - X_e * beta0).squaredNorm() / (2.0 * n);
      const double new_obj = (Y_e - X_e * beta1).squaredNorm() / (2.0 * n);
      const double denom = std::max(std::abs(old_obj), std::numeric_limits<double>::epsilon());
      dis = std::abs(new_obj - old_obj) / denom;
    }

    if (dis < thresh) {
      ++consecutive_below_thresh;
    } else {
      consecutive_below_thresh = 0;
    }

    alpha0 = alpha1;
    alpha1 = (1.0 + std::sqrt(1.0 + 4.0 * alpha0 * alpha0)) / 2.0;
  }

  return List::create(
    Named("beta") = Rcpp::wrap(beta1),
    Named("iter") = iter,
    Named("converged") = consecutive_below_thresh >= 5
  );
}
