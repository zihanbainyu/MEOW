function e = gg_eps(X)
% GG_EPS  Greenhouse-Geisser epsilon for one within-subjects effect.
%   X : n-by-k matrix of the k cell means that span the effect's within
%       space (one row per participant). NaN rows are dropped.
%   For a 2-level effect (k = 2) epsilon is 1 by definition.
X = X(all(~isnan(X),2),:);
k = size(X,2);
if k < 2, e = 1; return; end
S  = cov(X);
Z  = null(ones(1,k));            % k-by-(k-1) orthonormal basis, columns _|_ ones
T  = Z' * S * Z;                 % covariance of the orthonormal within-contrasts
ev = eig(T); ev = ev(ev > 0);
e  = sum(ev)^2 / ((k-1) * sum(ev.^2));
e  = min(max(e, 1/(k-1)), 1);    % bound to [1/(k-1), 1]
end
