function repro_pupil()
% REPRO_PUPIL  Pupil dilation statistics (Figure 3, Supplementary Figure 3).
%
% The cluster-based permutation that selects the 1.0-1.5 s analysis window runs
% upstream in S_pupil_cluster.m, because it needs the full per-trial time series.
% That script saves the per-subject window means this module analyses:
%   results/pupil_window_means.mat
%     lure_corr   [n x 3]  mean pupil, 1.0-1.5 s, correct lure trials   (comp,iso,nov)
%     lure_incorr [n x 3]  same, incorrect lure trials
%     lure_cond   [n x 3]  same, all lure trials (for condition contrasts)
%     targ_cond   [n x 3]  same, same-detection (target) trials
% To produce it without touching any existing script: run S_pupil_cluster, then
% run export_pupil_window_means (it reads the window means from the workspace).
%
% Reproduces: the 3x2 (condition x accuracy) on lure trials, the same-detection
% one-way ANOVA, and the condition / correctness post-hocs (one-tailed).

DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','results'); end
f = fullfile(DATA,'pupil_window_means.mat');
if ~isfile(f)
    fprintf(['[repro_pupil] %s not found.\n' ...
             '   Run S_pupil_cluster.m first (cluster-based window selection);\n' ...
             '   see README for the one-line save that exports the window means.\n'], f);
    return;
end
P = load(f); lc = P.lure_corr; li = P.lure_incorr; cd = P.lure_cond; tc = P.targ_cond;

fprintf('\n================ PUPIL DILATION (Figure 3) ================\n');
fprintf('\n-- Lure trials: 3x2 (encoding condition x discrimination accuracy) --\n');
rm_3x2('pupil lure 3x2', lc(:,1),li(:,1), lc(:,2),li(:,2), lc(:,3),li(:,3));

fprintf('   condition post-hocs (one-tailed):\n');
oh('    compared vs isolated', cd(:,1), cd(:,2), 'right');
oh('    isolated vs novel   ', cd(:,2), cd(:,3), 'right');
oh('    compared vs novel   ', cd(:,1), cd(:,3), 'right');
fprintf('   correct vs incorrect within condition (one-tailed):\n');
oh('    compared', lc(:,1), li(:,1), 'right');
oh('    isolated', lc(:,2), li(:,2), 'right');
oh('    novel   ', lc(:,3), li(:,3), 'right');

fprintf('\n-- Same-detection trials: one-way ANOVA (Supplementary Figure 3) --\n');
rm_oneway('pupil same-detect', tc(:,1), tc(:,2), tc(:,3));
end

function oh(lbl, a, b, tail)
[~,p,~,s] = ttest(a, b, 'Tail', tail);
fprintf('   %-24s t(%d) = %.2f, d = %.2f, %s\n', lbl, s.df, s.tstat, cohend(a,b), pstr(p));
end
