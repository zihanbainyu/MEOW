function q = bh_fdr(p)
% benjamini-hochberg fdr-adjusted p-values
p = p(:)'; n = numel(p);
[ps, idx] = sort(p);
qs = ps .* n ./ (1:n);
for i = n-1:-1:1, qs(i) = min(qs(i), qs(i+1)); end
q = zeros(size(p)); q(idx) = min(qs, 1);
end
