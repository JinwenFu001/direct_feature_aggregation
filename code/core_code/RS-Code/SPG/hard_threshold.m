function s = hard_threshold(w, lambda)
% MATLAB fallback for hard_threshold.c.
% Clips w into [-lambda, lambda].

if nargin < 2
    lambda = 1;
end

s = min(max(w, -lambda), lambda);

end
