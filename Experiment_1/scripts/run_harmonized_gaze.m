% run_harmonized_gaze.m
% Harmonised reporting for the three MAIN gaze 2x2 ANOVAs (condition x B2 accuracy):
%   A1-B1  pattern separation
%   B2-B1  pattern completion
%   A2-B1  predictive recall
% It rebuilds the per-subject cell means from the saved result files, applying the
% SAME per-comparison subject exclusions your original (commented-out) code used,
% then reports each with rm_2x2: integer df, ordinary p (all effects are 1 df, so
% Greenhouse-Geisser is a no-op), partial eta^2. F and SS come from fitrm/ranova,
% so they match your existing numbers.
%
% Requires rm_2x2.m, gg_eps.m, pstr.m on the path. Run from Experiment_1/scripts.
% SANITY CHECK the printed n against the manuscript: A1-B1 = 26, B2-B1 = 24, A2-B1 = 23.
% (These map to t(25), t(23), t(22) in the simple-effects tests.)

clear; clc;
res_dir = fullfile('..','results');

%% ---- A1-B1 (pattern separation) ----
S  = load(fullfile(res_dir,'gaze_reinstat_res_ab.mat'),'reinstat_res_ab');
sm = cellmeans_by_acc(S.reinstat_res_ab.ab_compared, S.reinstat_res_ab.ab_isolated, 'correct');
rm_2x2('A1-B1 gaze  (pattern separation)', sm(:,1), sm(:,2), sm(:,3), sm(:,4));

%% ---- B2-B1 (pattern completion) ----
F = load(fullfile(res_dir,'gaze_reinstat_res_full.mat'),'reinstat_res');
R = F.reinstat_res;
exc_bb = [609 606 608];
sm = cellmeans_by_acc(drop_subj(R.bb_compared, exc_bb), drop_subj(R.bb_isolated, exc_bb), 'correct');
rm_2x2('B2-B1 gaze  (pattern completion)', sm(:,1), sm(:,2), sm(:,3), sm(:,4));

%% ---- A2-B1 (predictive recall) ----
% Accuracy for an A2-B1 trial is the accuracy of its matching B2 trial (matched on
% tr_1b_b within participant), exactly as in the original analysis.
exc_ba = [609 606 608 618];
ba_c = attach_b2_correct(R.ba_compared, R.bb_compared);
ba_i = attach_b2_correct(R.ba_isolated, R.bb_isolated);
sm = cellmeans_by_acc(drop_subj(ba_c, exc_ba), drop_subj(ba_i, exc_ba), 'b2_correct');
rm_2x2('A2-B1 gaze  (predictive recall)', sm(:,1), sm(:,2), sm(:,3), sm(:,4));


%% ===================== helpers =====================
function sm = cellmeans_by_acc(Tc, Ti, accvar)
% Per-subject baseline-corrected gaze similarity, split by accuracy, per condition.
% Columns: comp/correct  comp/incorrect  iso/correct  iso/incorrect.
sids = unique([Tc.subj_id; Ti.subj_id]);
n = numel(sids); sm = nan(n,4);
for s = 1:n
    sid = sids(s);
    sm(s,1) = mean_acc(Tc, sid, accvar, 1);
    sm(s,2) = mean_acc(Tc, sid, accvar, 0);
    sm(s,3) = mean_acc(Ti, sid, accvar, 1);
    sm(s,4) = mean_acc(Ti, sid, accvar, 0);
end
end

function m = mean_acc(T, sid, accvar, val)
d = T.reinst_index(T.subj_id == sid & T.(accvar) == val);
if isempty(d), m = NaN; else, m = mean(d,'omitnan'); end
end

function T = drop_subj(T, exc)
T(ismember(T.subj_id, exc), :) = [];
end

function ba = attach_b2_correct(ba, bb)
ba.b2_correct = nan(height(ba),1);
for i = 1:height(ba)
    m = bb(bb.subj_id == ba.subj_id(i) & bb.tr_1b_b == ba.tr_1b_b(i), :);
    if height(m) == 1, ba.b2_correct(i) = m.correct; end
end
end
