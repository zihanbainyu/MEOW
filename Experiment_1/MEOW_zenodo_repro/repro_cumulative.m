function repro_cumulative()
% REPRO_CUMULATIVE  Cumulative gaze-reinstatement slopes (Figure 5).
% Source data:
%   results/cumu_reinstat_aa_cor.mat  (A2-A1 cumulative similarity by fixation)
%   results/cumu_reinstat_ba.mat      (A2-B1 cumulative similarity by fixation)
% Reproduces: participant-level slopes, one-sample tests against zero
% (one-tailed), and the 2x2 (pair type x encoding condition). Excludes subj 609.

DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','results'); end
exc = 609;

fprintf('\n================ CUMULATIVE REINSTATEMENT (Figure 5) ================\n');

% A2-A1
A = load(fullfile(DATA,'cumu_reinstat_aa_cor.mat'));
xc = (1:A.n_fix_to_plot)' - mean((1:A.n_fix_to_plot)');
sl_c = slopes(A.cumulative_results_comp, 5, A.n_fix_to_plot, xc);
sl_i = slopes(A.cumulative_results_iso,  5, A.n_fix_to_plot, xc);
[aa_c, aa_i] = subj(A.cumulative_results_comp.subj_id, sl_c, A.cumulative_results_iso.subj_id, sl_i, exc);
vs_zero('A2-A1', aa_c, aa_i);

% A2-B1
B = load(fullfile(DATA,'cumu_reinstat_ba.mat'));
sl_c = slopes(B.cumulative_results_comp, 4, B.n_fix_to_plot, xc);
sl_i = slopes(B.cumulative_results_iso,  4, B.n_fix_to_plot, xc);
[ba_c, ba_i] = subj(B.cumulative_results_comp.subj_id, sl_c, B.cumulative_results_iso.subj_id, sl_i, exc);
vs_zero('A2-B1', ba_c, ba_i);

% 2x2: pair type (A2-B1, A2-A1) x condition (compared, isolated)
fprintf('\n-- 2x2 ANOVA: pair type x encoding condition --\n');
rm_2x2('cumulative slopes', ba_c, ba_i, aa_c, aa_i);   % cc ci ic ii = BAcomp BAiso AAcomp AAiso
[~,pd,~,sd] = ttest(ba_c, ba_i, 'Tail','right');       % A2-B1 slope: compared > isolated
fprintf('   A2-B1 compared vs isolated slope: t(%d) = %.2f, d = %.2f, %s\n', sd.df, sd.tstat, cohend(ba_c,ba_i), pstr(pd));
end

% ================= local helpers =================
function sl = slopes(tbl, col0, nf, xc)
sl = nan(height(tbl),1);
for i=1:height(tbl)
    y = tbl{i, col0:col0+nf-1}'; v = ~isnan(y);
    if sum(v) >= 2, sl(i) = xc(v)\y(v); end
end
end
function [sc, si] = subj(sid_c, sl_c, sid_i, sl_i, exc)
agg = @(s,sl) grpstats(table(s(~isnan(sl)), sl(~isnan(sl)), 'VariableNames',{'subj_id','slope'}), 'subj_id','mean','DataVars','slope');
gc = agg(sid_c, sl_c); gc.Properties.VariableNames{'mean_slope'}='slope';
gi = agg(sid_i, sl_i); gi.Properties.VariableNames{'mean_slope'}='slope';
common = intersect(gc.subj_id, gi.subj_id); common(common==exc) = [];
[~,ic]=ismember(common,gc.subj_id); [~,ii]=ismember(common,gi.subj_id);
sc = gc.slope(ic); si = gi.slope(ii);
end
function vs_zero(lbl, sc, si)
[~,pc,~,tc] = ttest(sc, 0, 'Tail','right');
[~,pi,~,ti] = ttest(si, 0, 'Tail','right');
fprintf('\n-- %s cumulative slope vs zero (one-tailed) --\n', lbl);
fprintf('   compared: M = %.4f, SD = %.4f, t(%d) = %.2f, %s\n', mean(sc), std(sc), tc.df, tc.tstat, pstr(pc));
fprintf('   isolated: M = %.4f, SD = %.4f, t(%d) = %.2f, %s\n', mean(si), std(si), ti.df, ti.tstat, pstr(pi));
end
