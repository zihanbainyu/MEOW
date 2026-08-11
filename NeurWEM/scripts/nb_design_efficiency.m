function [eff, vif, rr] = nb_design_efficiency(onset, lab, run_end, TR, hp, C)
% NB_DESIGN_EFFICIENCY  fMRI design efficiency + collinearity for the
% unitized 1-back GLM (pairs modelled at the 2nd/judgment trial).
%
%   [eff, vif, rr] = nb_design_efficiency(onset, lab, run_end, TR, hp, C)
%     onset   : image onset (s) per trial, run clock (from nb_timeline_1back)
%     lab     : regressor id per trial: 1=Encoding 2=Same 3=Similar 4=New
%     run_end : run duration (s); sets the number of volumes
%     TR, hp  : repetition time and high-pass cutoff (s)
%     C       : nContrast-by-4 contrast matrix (cols = Encoding/Same/Similar/New)
%
%   eff : per-contrast efficiency = 1 / (c * inv(X'X) * c').  Higher = better
%         (efficiency is inversely proportional to the variance of the
%          contrast estimate, so higher efficiency = more power per scan-second).
%   vif : variance-inflation factor per regressor (near 1 = independent; >5 = concern)
%   rr  : [corr(Same,Similar) corr(Same,New) corr(Similar,New)] of the regressors
%
% Self-contained: canonical HRF and SPM-style DCT high-pass are computed
% locally (no SPM / Stats toolbox needed).
    dt = TR/16;                                  % microtime resolution
    hrf = local_hrf(dt);
    nScans = ceil(run_end/TR);
    nt = round(nScans*TR/dt);
    scan_idx = round((0:nScans-1)*TR/dt) + 1;

    X = zeros(nScans, 4);
    for k = 1:4
        stick = zeros(nt,1);
        oi = round(onset(lab==k)/dt) + 1;
        oi = oi(oi >= 1 & oi <= nt);
        stick(oi) = 1;
        c = conv(stick, hrf);
        X(:,k) = c(scan_idx);
    end

    % nuisance = intercept + DCT drift; residualise task regressors against it
    K = [ones(nScans,1) local_dct(nScans, TR, hp)];
    R = eye(nScans) - K*(K\eye(nScans));
    Xf = R*X;

    iM = inv(Xf'*Xf);
    eff = zeros(size(C,1),1);
    for ci = 1:size(C,1)
        eff(ci) = 1/(C(ci,:)*iM*C(ci,:)');
    end
    Rc = corrcoef(Xf);
    vif = diag(inv(Rc))';
    rr  = [Rc(2,3) Rc(2,4) Rc(3,4)];
end

function hrf = local_hrf(dt)
% SPM canonical double-gamma HRF sampled at dt (params [6 16 1 1 6 0 32]).
    p = [6 16 1 1 6 0 32];
    u = 0:ceil(p(7)/dt);
    g = @(x,h,l) (l.^h).*(x.^(h-1)).*exp(-l.*x)./gamma(h);
    hrf = g(u, p(1)/p(3), dt/p(3)) - g(u, p(2)/p(4), dt/p(4))/p(5);
    hrf = hrf(:); hrf = hrf/sum(hrf);
end

function K = local_dct(N, TR, hp)
% SPM-style discrete-cosine high-pass basis up to frequency 1/hp.
    n = fix(2*(N*TR)/hp + 1);
    K = zeros(N, max(n-1,0));
    t = (0:N-1)';
    for k = 2:n
        K(:,k-1) = sqrt(2/N)*cos(pi*(2*t+1)*(k-1)/(2*N));
    end
end
