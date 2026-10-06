\name{treeFA-selection}
\alias{dy_prox_simple_linear}
\alias{dy_prox_simple_logistic}
\alias{grid.select_linear}
\alias{grid.select_logistic}
\alias{cv.select_linear}
\alias{cv.select_logistic}
\title{Joint Tree Aggregation and Individual Feature Selection}
\description{
Fit linear or binary-response models with both the original centered tree
penalty and an L1 penalty on the feature coefficients, using Davis--Yin
three-operator splitting. These are separate interfaces; the existing treeFA
solvers, paths and cross-validation functions retain their original behavior.
}
\usage{
dy_prox_simple_linear(Y, X, tree_result, lambda, tau,
  intercept = FALSE, ridge_param = 0, step_size = NULL,
  abs_tol = 1e-7, rel_tol = 1e-6, max_iter = 10000L,
  init_z = NULL, keep_history = FALSE)

dy_prox_simple_logistic(Y, X, tree_result, lambda, tau,
  intercept = FALSE, ridge_param = 0, step_size = NULL,
  abs_tol = 1e-7, rel_tol = 1e-6, max_iter = 10000L,
  init_z = NULL, keep_history = FALSE)

grid.select_linear(Y, X, tree_df, param_grid = NULL,
  lambda_seq = NULL, tau_seq = NULL, intercept = FALSE,
  ridge_param = 0, step_size = NULL, abs_tol = 1e-7,
  rel_tol = 1e-6, max_iter = 10000L, warm_start = TRUE)

grid.select_logistic(Y, X, tree_df, param_grid = NULL,
  lambda_seq = NULL, tau_seq = NULL, intercept = FALSE,
  ridge_param = 0, step_size = NULL, abs_tol = 1e-7,
  rel_tol = 1e-6, max_iter = 10000L, warm_start = TRUE)

cv.select_linear(Y, X, tree_df, param_grid = NULL,
  lambda_seq = NULL, tau_seq = NULL, intercept = FALSE,
  ridge_param = 0, step_size = NULL, abs_tol = 1e-7,
  rel_tol = 1e-6, max_iter = 10000L, warm_start = TRUE,
  folds = 5L, foldid = NULL, retry_max_iter = NULL, verbose = FALSE)

cv.select_logistic(Y, X, tree_df, param_grid = NULL,
  lambda_seq = NULL, tau_seq = NULL, intercept = FALSE,
  ridge_param = 0, step_size = NULL, abs_tol = 1e-7,
  rel_tol = 1e-6, max_iter = 10000L, warm_start = TRUE,
  folds = 5L, foldid = NULL, retry_max_iter = NULL, verbose = FALSE)
}
\arguments{
\item{Y}{Finite numeric response vector; 0/1 for logistic models.}
\item{X}{Finite numeric design matrix. Columns are not standardized.}
\item{tree_result}{A tree data frame with \code{node}, \code{parent}, and
\code{weight}, or a layered list produced by
\code{gather_leaf_nodes_per_non_leaf()}. Leaves must correspond to columns
\code{1:p}. Weights must be finite and non-negative.}
\item{tree_df}{Tree data frame, with a single root whose parent is \code{NA}.}
\item{lambda, tau}{Finite non-negative tree and L1 penalty strengths.}
\item{intercept}{Fit an unpenalized intercept. Gaussian data are centered;
the logistic intercept is optimized jointly, without centering the response.}
\item{ridge_param}{Finite non-negative ridge strength. Its contribution is
\eqn{\mathrm{ridge\_param}\|\beta\|_2^2/(2n)}, excluding the intercept.}
\item{step_size}{Optional positive fixed step satisfying \eqn{hL<2} for the
computed Lipschitz upper bound. Defaults to \eqn{1/L}, or 1 when \eqn{L=0}.}
\item{abs_tol, rel_tol}{Finite non-negative absolute and relative tolerances;
at least one must be positive. Applied to both splitting disagreement and
step-scaled residual.}
\item{max_iter}{Positive integer maximum number of operator evaluations.}
\item{init_z}{Optional splitting-state vector, not a coefficient vector.
Its length is \eqn{p}, except logistic with intercept where it is \eqn{p+1}
and the intercept is first. Defaults to zero. Use \code{fit$state$z} to resume.}
\item{keep_history}{Retain objective and residual diagnostics at every iteration.}
\item{param_grid}{Nonempty data frame or matrix with \code{lambda} and
\code{tau} columns. Candidate order and duplicate rows are preserved.}
\item{lambda_seq, tau_seq}{Both sequences must be supplied when
\code{param_grid} is absent. Their Cartesian product has lambda varying fastest.
No automatic response-dependent penalty grid is constructed.}
\item{warm_start}{Carry the previous candidate's splitting state to the next.}
\item{folds}{Number of folds, at least two. Ignored if \code{foldid} is supplied.}
\item{foldid}{Optional vector of consecutive positive fold labels beginning at
one. Every validation and training set must be nonempty. Logistic training
sets must each contain both classes; default logistic folds are stratified.
Use \code{set.seed()} before a call for reproducible random folds.}
\item{retry_max_iter}{Additional evaluation budget for one continuation of a
nonconverged CV fit. Defaults to twice \code{max_iter}. Applies to the final
full-data refit too.}
\item{verbose}{Print CV progress.}
}
\details{
The objective is \eqn{g(\beta_0,\beta)+\lambda\Omega(\beta)+\tau\|\beta\|_1},
where \eqn{\Omega(\beta)=\sum_G w_G\|\beta_G-\bar\beta_G\mathbf{1}\|_2}.
The smooth loss is squared error divided by \eqn{2n}, or average negative
binomial log-likelihood, plus the specified ridge term. The intercept is
excluded from all penalties. Tree structure is validated before solving.

The iteration is \eqn{u=\mathrm{soft}_{h\tau}(z)},
\eqn{v=\mathrm{prox}_{h\lambda\Omega}(2u-z-h\nabla g(u))},
\eqn{z\leftarrow z+v-u}, with identity proximal operations on the intercept.
The returned coefficients are \eqn{u}, not \eqn{v} or \eqn{z}.
The returned \code{state$z} is the same state that generated the returned
coefficients and diagnostics; resuming repeats its evaluation once.

Let \eqn{d=\|v-u\|_\infty},
\eqn{a=(z-u)/h} and \eqn{b=(2u-z-h\nabla g(u)-v)/h}.
Convergence requires both
\eqn{d\leq\mathrm{abs\_tol}+\mathrm{rel\_tol}\max(1,\|u\|_\infty,\|v\|_\infty)}
and
\eqn{d/h\leq\mathrm{abs\_tol}+\mathrm{rel\_tol}\max(1,\|\nabla g(u)\|_\infty,\|a\|_\infty,\|b\|_\infty)}.
These are numerical fixed-point criteria, not an objective-gap certificate.
Objective values need not decrease monotonically. Finite-iteration \eqn{u}
can have exact zeros without exactly equal aggregated coefficients; no
coefficient rounding or group postprocessing is performed.

The step bound is obtained using a spectral norm computed by SVD, with a small
safety margin. The Gaussian bound is \eqn{(\|X\|_2^2+\rho)/n}; the logistic
bound is \eqn{\|A\|_2^2/(4n)+\rho/n}, with \eqn{A=[\mathbf{1},X]} when an
intercept is fitted. Preparation is cached within a path. Iterations use
matrix-vector gradients and an in-place tree prox. A fit returns a warning
and status if it cannot satisfy the stopping criteria within the budget.

Logistic convergence theory presupposes a finite minimizer. Single-class
responses with an intercept are rejected; unpenalized separating directions
can cause nonexistence in other cases, which the numerical convergence flag
does not certify or rule out.

CV uses the same explicit candidates and folds throughout. Centering and
step-size preparation use only training observations. Scores are held-out
mean squared errors or mean negative binomial log-likelihoods, without
penalties. Fold errors are weighted by fold size. A candidate failing any
fold after retry has \code{eligible = FALSE} and \code{cv.error = Inf}; errors
are never averaged over only the successful folds. The first row wins exact
ties. An error is raised if no candidate is eligible or the selected full-data
refit still fails. Ridge retains its stated \eqn{1/n} convention within each fold.
}
\value{
Single fits return \code{beta}, \code{beta0}, \code{lambda}, \code{tau},
\code{ridge_param}, \code{family}, \code{intercept}, \code{iter},
\code{converged}, \code{status}, \code{residual}, \code{residual_scaled},
\code{objective}, \code{step_size}, \code{lipschitz}, \code{state},
\code{n_nonzero}, and optional \code{history}.

Paths return a \eqn{p\times K} \code{beta} matrix, \code{beta0},
\code{param_grid}, per-candidate \code{diagnostics}, and \code{n_nonzero}.

CV returns \code{selected.lambda}, \code{selected.tau}, \code{selected.index},
the final \code{beta}, \code{beta0} and \code{fit}, \code{param_grid},
\code{foldid}, \code{fold.error}, \code{cv.error}, \code{eligible}, and
\code{fold.diagnostics}. Iteration totals include both attempts when retried.
}
\references{
Davis, D. and Yin, W. (2017). A Three-Operator Splitting Scheme and its
Optimization Applications. Set-Valued and Variational Analysis, 25, 829--858.
\doi{10.1007/s11228-017-0421-z}.
}
\examples{
set.seed(42)
tree <- data.frame(node = 1:3, parent = c(3, 3, NA), weight = c(0, 0, 1))
X <- matrix(rnorm(80), 40, 2)
Y <- 1 + X \%*\% c(0.8, 0.8) + rnorm(40, sd = 0.3)
fit <- dy_prox_simple_linear(as.numeric(Y), X, tree,
                             lambda = 0.1, tau = 0.05, intercept = TRUE)
fit$beta
grid <- expand.grid(lambda = c(0, 0.1), tau = c(0, 0.05))
cv <- cv.select_linear(as.numeric(Y), X, tree, param_grid = grid,
                       intercept = TRUE, folds = 3)
c(lambda = cv$selected.lambda, tau = cv$selected.tau)
}
