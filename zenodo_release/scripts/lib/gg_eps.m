function e = gg_eps(X)
% greenhouse-geisser epsilon for one within effect; X is n x k cell means
X = X(all(~isnan(X),2),:);
k = size(X,2);
if k < 2, e = 1; return; end
S = cov(X);
Z = null(ones(1,k));            % orthonormal contrast basis
T = Z' * S * Z;
ev = eig(T); ev = ev(ev > 0);
e = sum(ev)^2 / ((k-1) * sum(ev.^2));
e = min(max(e, 1/(k-1)), 1);
end
