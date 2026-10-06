# Joint feature aggregation and selection

The new `dy_prox_*`, `grid.select_*`, and `cv.select_*` interfaces minimize

\[
g(\beta_0,\beta)+\lambda\sum_G w_G\|\beta_G-\bar\beta_G\mathbf 1\|_2
+\tau\|\beta\|_1.
\]

They use a separate Davis–Yin solver. The original `acc_prox_simple_*`,
`grid.simple_linear`, `grid.logistic`, `cv.simple_linear`, and `cv.logistic`
interfaces and algorithms are unchanged. In particular, neither `tau = 0`
nor any other argument to a new interface dispatches through or modifies the
old solver. The new tree prox uses the same mathematical operations in its
own implementation with an in-place buffer.

## Example

```r
library(treeFA)
set.seed(42)
tree <- data.frame(
  node = 1:7,
  parent = c(5, 5, 6, 6, 7, 7, NA),
  weight = c(0, 0, 0, 0, 1, 1, 0.5)
)
X <- matrix(rnorm(320), 80, 4)
Y <- as.numeric(1 + X %*% c(0.8, 0.8, 0, 0) + rnorm(80, sd = 0.3))

fit <- dy_prox_simple_linear(
  Y, X, tree, lambda = 0.1, tau = 0.05,
  intercept = TRUE, keep_history = TRUE
)
fit$beta
fit$converged
fit$residual_scaled
prediction <- drop(fit$beta0 + X %*% fit$beta)

params <- expand.grid(lambda = c(0, 0.05, 0.2), tau = c(0, 0.03, 0.1))
path <- grid.select_linear(Y, X, tree, param_grid = params, intercept = TRUE)

foldid <- sample(rep(1:5, length.out = length(Y)))
cv <- cv.select_linear(
  Y, X, tree, param_grid = params, intercept = TRUE, foldid = foldid
)
c(lambda = cv$selected.lambda, tau = cv$selected.tau)
cv$fit$beta
cv$fold.error
```

For a binary response encoded as 0/1, use `dy_prox_simple_logistic()`,
`grid.select_logistic()`, or `cv.select_logistic()`. Logistic intercepts are
optimized jointly and are not penalized. Columns are never automatically
standardized; both penalties act on the supplied feature scale.

Single fits accept either a tree data frame or the existing
`gather_leaf_nodes_per_non_leaf(tree)$list_result`. Trees require finite
nonnegative weights, one root, and leaves numbered `1:ncol(X)`. A gathered
list must have valid nested groups ordered from descendants to ancestors.

## Numerical conventions

- Gaussian loss is `sum((Y - beta0 - X %*% beta)^2) / (2*n)`;
  logistic loss is mean negative binomial log-likelihood.
- Optional ridge contributes `ridge_param * sum(beta^2) / (2*n)`.
  The default is zero, matching the objective above.
- Cold starts use `z = 0`. To resume, pass `init_z = fit$state$z` and, when
  appropriate, `step_size = fit$step_size`. This state is not the coefficient
  vector. Logistic states with intercept put the intercept first. The saved
  state generated the returned coefficients and residuals, so continuation
  repeats one evaluation before progressing.
- The default step uses a cached SVD-based Lipschitz bound. A manually supplied
  positive step is checked against that bound. `L = 0` uses step 1 by default.
- Both `residual = max(abs(v-u))` and `residual_scaled = residual/step_size`
  must meet absolute/relative tolerances. These are fixed-point diagnostics,
  not a certified bound on coefficient or objective error. Objective values
  are not assumed to decrease at every iteration.
- Returned coefficients are `u = soft(z, step_size*tau)`. They can contain
  exact zeros. At finite tolerance, tree-aggregated coefficients need not be
  exactly equal; no rounding or automatic aggregation-group extraction occurs.
- Inspect `converged` and `status`. Single fits and paths warn on failure.
  Logistic problems can lack a finite minimizer along unpenalized separating
  directions; single-class responses with an intercept are explicitly rejected.

## Paths and cross-validation

Supply either `param_grid` or both `lambda_seq` and `tau_seq`. There is no
automatic parameter endpoint calculation. Paths preserve candidate order and
duplicates, cache preprocessing, and optionally carry `z` between candidates.
`beta[, k]` corresponds to `param_grid[k, ]`. Every candidate has diagnostics.

CV uses the same explicit candidates in all folds and performs centering and
step-size preparation on each training fold. Binomial default folds are
stratified; all binomial training folds must contain both classes. A supplied
`foldid` determines the number of folds and must use consecutive labels.

Held-out MSE or binomial log-loss is averaged over observations, so unequal
fold sizes are properly weighted. Nonconverged fits get one continuation with
`retry_max_iter` additional evaluations (default `2 * max_iter`). A candidate
that still fails any fold is excluded rather than scored using only its
successful folds. Exact ties select the first candidate. The final full-data
refit must converge, or CV raises an error.

Full argument and return-value details: `help("treeFA-selection")`.
