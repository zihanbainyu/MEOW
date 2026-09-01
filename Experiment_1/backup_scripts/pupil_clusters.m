function C = pupil_clusters(pupil)
% one-time: cluster-based permutation on the pupil time-courses (a-b lure trials).
% omnibus condition effect and the three pairwise contrasts. returns, per test,
% rows of [t0 t1 mass p] (seconds). the repro script only prints these.
rng(1);
t = pupil.t; P = pupil.pup;
C.omnibus  = omnibus(t, P.ab_com, P.ab_iso, P.ab_nov);
C.comp_iso = pairwise(t, P.ab_com, P.ab_iso);
C.comp_nov = pairwise(t, P.ab_com, P.ab_nov);
C.iso_nov  = pairwise(t, P.ab_iso, P.ab_nov);
end

function M = omnibus(t, A, B, Cc)
np = 1000; nt = numel(t);
Fobs = arrayfun(@(ti) rmF([A(:,ti) B(:,ti) Cc(:,ti)]), 1:nt);
nfull = sum(~isnan(A(:,1)) & ~isnan(B(:,1)) & ~isnan(Cc(:,1)));
fcrit = finv(0.95, 2, (nfull-1)*2);
[cl, cs] = clusters(Fobs, fcrit, Fobs);
mx = zeros(np,1);
for p = 1:np
    Fp = arrayfun(@(ti) rmF(shuffle3([A(:,ti) B(:,ti) Cc(:,ti)])), 1:nt);
    [~, csp] = clusters(Fp, fcrit, Fp);
    if ~isempty(csp), mx(p) = max(csp); end
end
M = pack(t, cl, cs, mx);
end

function M = pairwise(t, A, B)
np = 1000; D = A - B;
nfull = sum(~isnan(D(:,1)));
tcrit = tinv(0.975, nfull-1);
Tobs = tstat(D);
[cl, cs] = clusters(abs(Tobs), tcrit, Tobs);
mx = zeros(np,1);
for p = 1:np
    Tp = tstat(D .* sign(rand(size(D,1),1) - 0.5));
    [~, csp] = clusters(abs(Tp), tcrit, Tp);
    if ~isempty(csp), mx(p) = max(abs(csp)); end
end
M = pack(t, cl, cs, mx);
end

function M = pack(t, cl, cs, mx)
M = zeros(numel(cl),4);
for i = 1:numel(cl)
    M(i,:) = [t(cl{i}(1)), t(cl{i}(end)), cs(i), mean(mx >= abs(cs(i)))];
end
if ~isempty(M), [~,o] = sort(abs(M(:,3)),'descend'); M = M(o,:); end
end

function F = rmF(M)
M = M(all(~isnan(M),2), :); [n,k] = size(M);
if n < 2, F = 0; return; end
g = mean(M(:)); cm = mean(M,1); sm = mean(M,2);
SSc = n*sum((cm-g).^2); SSs = k*sum((sm-g).^2);
SSe = sum((M(:)-g).^2) - SSc - SSs;
F = (SSc/(k-1)) / (SSe/((k-1)*(n-1)));
end

function T = tstat(D)
T = nan(1,size(D,2));
for ti = 1:size(D,2)
    d = D(~isnan(D(:,ti)),ti); n = numel(d);
    if n > 1, T(ti) = mean(d)/(std(d)/sqrt(n)); end
end
end

function [cl, cs] = clusters(stat, thr, signed)
above = stat > thr; above(isnan(above)) = false;
d = diff([0 above 0]); st = find(d==1); en = find(d==-1)-1;
cl = cell(1,numel(st)); cs = zeros(1,numel(st));
for c = 1:numel(st), cl{c} = st(c):en(c); cs(c) = sum(signed(cl{c})); end
end

function M = shuffle3(M)
for s = 1:size(M,1), M(s,:) = M(s, randperm(3)); end
end
