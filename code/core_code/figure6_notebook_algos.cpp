#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// Function to scale matrix
NumericMatrix scale(const NumericMatrix& X, bool center = true, bool scale = true) {
  int n = X.nrow(), p = X.ncol();
  NumericMatrix X_scaled = clone(X);
  NumericVector col_means(p);
  for (int j = 0; j < p; ++j) {
    col_means[j] = mean(X(_, j));
    for (int i = 0; i < n; ++i) {
      X_scaled(i, j) -= col_means[j];
    }
  }
  return X_scaled;
}

// [[Rcpp::export]]
Rcpp::NumericVector eigen_to_rcpp(const Eigen::VectorXd &vec) {
  // Create an Rcpp::NumericVector from the Eigen::VectorXd
  Rcpp::NumericVector r_vec(vec.data(), vec.data() + vec.size());
  return r_vec;
}

// [[Rcpp::export]]
NumericVector one_layer(NumericVector eta, List groups, NumericVector weights, double lambda) {
  NumericVector beta = clone(eta);
  int n = weights.size();
  
  for (int i = 0; i < n; ++i) {
    double weight = weights[i];
    IntegerVector ind = groups[i];
    NumericVector mu = eta[ind - 1]; // R is 1-indexed, C++ is 0-indexed
    double mean_mu = mean(mu);
    double temp=sqrt(sum(pow(mu - mean_mu, 2)));
    double d = 100;
    if (temp > 0) {
      d=weight * lambda / temp;
    }
    
    if (d >= 1) {
      for (int j = 0; j < ind.size(); ++j) {
        beta[ind[j] - 1] = mean_mu;
      }
    } else {
      for (int j = 0; j < ind.size(); ++j) {
        beta[ind[j] - 1] = (1 - d) * mu[j] + d * mean_mu;
      }
    }
  }
  
  return beta;
}

// [[Rcpp::export]]
NumericVector prox_tree(NumericVector eta, double lambda, List tree_list) {
  List groups = tree_list["groups"];
  List weights = tree_list["weights"];
  int depth = groups.size();
  NumericVector new_eta = clone(eta);
  
  for (int i = 0; i < depth; ++i) {
    new_eta = one_layer(new_eta, groups[i], weights[i], lambda);
  }
  
  return new_eta;
}