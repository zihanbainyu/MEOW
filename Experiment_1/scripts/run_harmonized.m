% run_harmonized.m
% One consistent convention for every RM-ANOVA in the manuscript:
%   * integer degrees of freedom for every F
%   * Greenhouse-Geisser applied to the p-value of any effect with numerator
%     df > 1 (the 3-level "condition" effects); epsilon reported
%   * no correction for 1-df (2-level) effects; ordinary p
%   * effect size = partial eta^2 = SS_effect / (SS_effect + SS_error)
% F and SS are computed with fitrm/ranova, so they match your existing numbers;
% only the df/GG/epsilon presentation is harmonised.
%
% Put gg_eps.m, pstr.m, rm_oneway.m, rm_2x2.m, rm_3x2.m on the MATLAB path
% (e.g. in Experiment_1/scripts). Run this from Experiment_1/scripts.

clear; clc;

%% ---------- BEHAVIOURAL (runs as-is) ----------
res_dir = fullfile('..','results');
load(fullfile(res_dir,'all_subjs_stats.mat'),'all_subjs');
gv = @(f1,f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs)';

fprintf('\n== 2-back one-way RM-ANOVAs (condition: compared / isolated / novel) ==\n');
rm_oneway('discrimination idx', gv('two','ldi_comp'),    gv('two','ldi_iso'),    gv('two','ldi_nov'));
rm_oneway('same-item d''',       gv('two','dprime_comp'), gv('two','dprime_iso'), gv('two','dprime_nov'));
rm_oneway('RT lure (correct)',   gv('two','rt_AB_comp'),  gv('two','rt_AB_iso'),  gv('two','rt_AB_nov'));
rm_oneway('RT target (correct)', gv('two','rt_AA_comp'),  gv('two','rt_AA_iso'),  gv('two','rt_AA_nov'));

%% ---------- PUPIL & GAZE (call the reporters where the cell means exist) ----------
% These use different result files / subject subsets, so add the calls below
% inside your existing scripts at the point where the per-subject cell-mean
% vectors are built. The reporters recompute F/SS from those same vectors, so
% the values stay identical to your pipeline.
%
% S_pupil_cluster.m  (lure 3x2 -- pass the same 6 vectors it already builds):
%     rm_3x2('pupil lure 3x2', cc,ci, ic,ii, nc,ni);
% S_pupil_cluster.m  (same-detection one-way, the 3 target-condition vectors):
%     rm_oneway('pupil same-detect', tc, ti, tn);
%
% Gaze 2x2 scripts (S_a1b1.m etc. -- the 'sm' matrix has cols cc ci ic ii):
%     rm_2x2('A1-B1 gaze', sm(:,1), sm(:,2), sm(:,3), sm(:,4));
%     rm_2x2('B2-B1 gaze', sm(:,1), sm(:,2), sm(:,3), sm(:,4));
%     rm_2x2('A2-B1 gaze', sm(:,1), sm(:,2), sm(:,3), sm(:,4));
%
% Cumulative 2x2 (S_cumu_reinstat.m) -- already integer-df/ordinary-p, unchanged:
%     rm_2x2('cumulative slopes', ac, ai, bc, bi);
