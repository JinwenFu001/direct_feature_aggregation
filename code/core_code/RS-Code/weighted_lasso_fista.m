function beta = weighted_lasso_fista(Y, X, gamma, penalty_factor, option)
% Solve 1/2 ||Y - X beta||^2 + gamma * sum_j penalty_factor_j |beta_j|.

[~, J] = size(X);
Y = Y(:);
penalty_factor = penalty_factor(:);
if length(penalty_factor) ~= J
    error('penalty_factor must have length size(X, 2).')
end

if isfield(option, 'maxiter')
    maxiter = option.maxiter;
else
    maxiter = 10000;
end

if isfield(option, 'tol')
    tol = option.tol;
else
    tol = 1e-7;
end

if isfield(option, 'b_init')
    beta = option.b_init(:);
    if length(beta) ~= J
        beta = zeros(J, 1);
    end
else
    beta = zeros(J, 1);
end

XX = X' * X;
XY = X' * Y;
L = eigs(XX, 1);
if ~isfinite(L) || L <= 0
    L = 1;
end

z = beta;
t = 1;
old_obj = Inf;

for iter = 1:maxiter
    grad = XX * z - XY;
    beta_new = soft_threshold_weighted(z - grad / L, gamma * penalty_factor / L);

    t_new = (1 + sqrt(1 + 4 * t^2)) / 2;
    z = beta_new + ((t - 1) / t_new) * (beta_new - beta);

    residual = Y - X * beta_new;
    obj = sum(residual .^ 2) / 2 + gamma * sum(penalty_factor .* abs(beta_new));
    if iter > 10 && abs(obj - old_obj) / max(1, abs(old_obj)) < tol
        beta = beta_new;
        break;
    end

    beta = beta_new;
    t = t_new;
    old_obj = obj;
end

end

function out = soft_threshold_weighted(x, lambda)
out = sign(x) .* max(abs(x) - lambda, 0);
end
