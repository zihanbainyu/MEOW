function repro_behavioral()
% REPRO_BEHAVIORAL  Behavioural statistics (Figure 2, Supplementary Figure 1).
% Reproduces: two-back discrimination index, same-item d', reaction times,
% the WM-EM correlations, one-back accuracy/RT, and recognition d'.
% Source data: results/all_subjs_stats.mat  (variable: all_subjs).
% Reporting convention: one-way RM-ANOVA, integer df, Greenhouse-Geisser p
% (numerator df > 1), partial eta^2. Directional predictions one-tailed.
%
% Run repro_config first (adds lib/ + sets DATA path), or run run_all.

DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','results'); end
load(fullfile(DATA,'all_subjs_stats.mat'), 'all_subjs');
gv = @(f1,f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs)';

% ---- pull condition vectors ----
ldi = {gv('two','ldi_comp'),    gv('two','ldi_iso'),    gv('two','ldi_nov')};
dpr = {gv('two','dprime_comp'), gv('two','dprime_iso'), gv('two','dprime_nov')};
rtL = {gv('two','rt_AB_comp'),  gv('two','rt_AB_iso'),  gv('two','rt_AB_nov')};
rtT = {gv('two','rt_AA_comp'),  gv('two','rt_AA_iso'),  gv('two','rt_AA_nov')};

fprintf('\n================ BEHAVIOURAL (Figure 2) ================\n');
fprintf('\n-- Two-back one-way RM-ANOVAs (condition: compared/isolated/novel) --\n');
rm_oneway('discrimination idx', ldi{:});
posthoc_dir('  discrimination idx', ldi{:}, +1);   % higher = better -> one-tailed
rm_oneway('same-item d''',       dpr{:});           % null control (BF below)
rm_oneway('RT lure (correct)',   rtL{:});
posthoc_dir('  RT lure',          rtL{:}, -1);       % faster = better -> one-tailed
rm_oneway('RT target (correct)', rtT{:});           % null control

% ---- Bayes factors for the two null controls (optional; needs bayesFactor toolbox) ----
bf_null('same-item d''',      dpr{:});
bf_null('RT target',          rtT{:});

% ---- WM-EM correlations (one-tailed, Figure 2e,f) ----
fprintf('\n-- WM-EM correlations (one-tailed) --\n');
wm  = (dpr{1}+dpr{2}+dpr{3})/3;                       % two-back same-item d'
emb = (ldi{1}+ldi{2})/2 - ldi{3};                     % episodic benefit
em  = (gv('rec','d_comp') + gv('rec','d_iso'))/2;      % recognition d'
corr1('EM benefit  x  WM-d''', emb, wm);
corr1('EM benefit  x  EM-d''', emb, em);
corr1('WM-d''       x  EM-d''', wm,  em);

% ---- One-back performance (Supplementary Figure 1), BH-FDR ----
fprintf('\n================ ONE-BACK (Supplementary Figure 1) ================\n');
a_sam=gv('one','acc_same'); a_sim=gv('one','acc_sim'); a_new=gv('one','acc_new');
r_sam=gv('one','rt_same');  r_sim=gv('one','rt_sim');
[~,p1,~,s1]=ttest(a_sam,a_sim); [~,p2,~,s2]=ttest(a_sim,a_new);
[~,p3,~,s3]=ttest(a_sam,a_new); [~,p4,~,s4]=ttest(r_sam,r_sim);
q=bh_fdr([p1 p2 p3 p4]);
prow('acc same>similar', s1, cohend(a_sam,a_sim), p1, q(1));
prow('acc similar>new ', s2, cohend(a_sim,a_new), p2, q(2));
prow('acc same>new    ', s3, cohend(a_sam,a_new), p3, q(3));
prow('RT  same<similar', s4, cohend(r_sam,r_sim), p4, q(4));

% ---- Recognition d' (Supplementary Figure 1d), BH-FDR ----
fprintf('\n-- Recognition d'' --\n');
rc=gv('rec','d_comp'); ri=gv('rec','d_iso');
[~,pp1,~,ss1]=ttest(rc); [~,pp2,~,ss2]=ttest(ri); [~,pp3,~,ss3]=ttest(rc,ri);
qq=bh_fdr([pp1 pp2 pp3]);
prow('compared > 0     ', ss1, cohend(rc), pp1, qq(1));
prow('isolated > 0     ', ss2, cohend(ri), pp2, qq(2));
prow('compared vs iso  ', ss3, cohend(rc,ri), pp3, qq(3));
end

% ================= local helpers =================
function posthoc_dir(lbl, xc, xi, xn, sgn)
% directional (one-tailed) condition contrasts, BH-FDR over the family of 3.
% sgn = +1 : predict compared > isolated > novel ; sgn = -1 : predicted smaller.
tail = 'right'; if sgn < 0, tail = 'left'; end
[~,p1,~,s1]=ttest(xc,xi,'Tail',tail);
[~,p2,~,s2]=ttest(xi,xn,'Tail',tail);
[~,p3,~,s3]=ttest(xc,xn,'Tail',tail);
q=bh_fdr([p1 p2 p3]);
fprintf('%s post-hocs (one-tailed, BH-FDR):\n', lbl);
prow('    compared vs isolated', s1, cohend(xc,xi), p1, q(1));
prow('    isolated vs novel   ', s2, cohend(xi,xn), p2, q(2));
prow('    compared vs novel   ', s3, cohend(xc,xn), p3, q(3));
end

function corr1(lbl, x, y)
v=~isnan(x)&~isnan(y); [r,p]=corr(x(v),y(v),'Tail','right');
fprintf('  %-24s r(%d) = %.2f, p = %.3f, n = %d\n', lbl, sum(v)-2, r, p, sum(v));
end

function prow(lbl, s, d, p, q)
fprintf('  %-22s t(%d) = %.2f, d = %.2f, %s, q = %.3f\n', lbl, s.df, s.tstat, d, pstr(p), q);
end

function bf_null(lbl, xc, xi, xn)
try
    n=numel(xc); T=table((1:n)',xc,xi,xn,'VariableNames',{'subj','c1','c2','c3'});
    L=stack(T,{'c1','c2','c3'},'NewDataVariableName','y','IndexVariableName','cond');
    L.subj=categorical(L.subj); L.cond=categorical(L.cond);
    b=bf.anova(L,'y~cond','treatAsRandom',{'subj'},'verbose',false);
    if b<1, fprintf('  BF01 (%s) = %.2f\n', lbl, 1/b); else, fprintf('  BF10 (%s) = %.2f\n', lbl, b); end
catch
    fprintf('  [BF skipped for %s: bayesFactor toolbox not on path]\n', lbl);
end
end
