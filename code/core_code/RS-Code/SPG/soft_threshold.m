function [s, nz] = soft_threshold(w, lambda)
% MATLAB fallback for soft_threshold.c.
% Returns sign(w) * max(abs(w) - lambda, 0) and the number of nonzeros.

if nargin < 2
    lambda = 1;
end

s = sign(w) .* max(abs(w) - lambda, 0);
nz = nnz(s);

end
