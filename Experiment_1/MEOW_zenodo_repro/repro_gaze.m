function repro_gaze()
% REPRO_GAZE  Gaze-reinstatement statistics (Figure 4d-f).
%   A1-B1 pattern separation, B2-B1 pattern completion, A2-B1 predictive recall.
% Source data (group-level, baseline-corrected gaze similarity per trial):
%   results/gaze_reinstat_res_ab.mat   (reinstat_res_ab: ab_compared/ab_isolated)
%   results/gaze_reinstat_res_full.mat (reinstat_res:   bb_*, ba_*)
% Each 2x2 = encoding condition x subsequent B2 accuracy (integer df, partial eta^2).
% Per-comparison exclusions reproduce the manuscript samples (n = 26 / 24 / 23).
% Simple effects are one-tailed in the predicted direction.

DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','results'); end

fprintf('\n================ GAZE REINSTATEMENT (Figure 4) ================\n');

% ---- A1-B1  (pattern separation: correct trials show LOWER similarity) ----
S = load(fullfile(DATA,'gaze_reinstat_res_ab.mat'),'reinstat_res_ab');
sm = cellmeans(S.reinstat_res_ab.ab_compared, S.reinstat_res_ab.ab_isolated, 'correct');
report('A1-B1 (pattern separation)', sm, 'left');   % correct < incorrect

% ---- B2-B1 and A2-B1 from the full reinstatement file ----
F = load(fullfile(DATA,'gaze_reinstat_res_full.mat'),'reinstat_res'); R = F.reinstat_res;

% B2-B1 (pattern completion: correct trials show HIGHER similarity)
sm = cellmeans(drop_subj(R.bb_compared,[609 606 608]), drop_subj(R.bb_isolated,[609 606 608]), 'correct');
report('B2-B1 (pattern completion)', sm, 'right');   % correct > incorrect

% A2-B1 (predictive recall; accuracy from the matching B2 trial)
ba_c = attach_b2(R.ba_compared, R.bb_compared);
ba_i = attach_b2(R.ba_isolated, R.bb_isolated);
sm = cellmeans(drop_subj(ba_c,[609 606 608 618]), drop_subj(ba_i,[609 606 608 618]), 'b2_correct');
report('A2-B1 (predictive recall)', sm, 'right');    % correct > incorrect
end

% ================= local helpers =================
function report(label, sm, tail)
rm_2x2(label, sm(:,1), sm(:,2), sm(:,3), sm(:,4));
[~,p1,~,s1] = ttest(sm(:,1), sm(:,2), 'Tail', tail);   % compared: correct vs incorrect
[~,p2,~,s2] = ttest(sm(:,3), sm(:,4), 'Tail', tail);   % isolated: correct vs incorrect
fprintf('   simple effects (one-tailed):\n');
fprintf('     compared correct vs incorrect: t(%d) = %.2f, %s\n', s1.df, s1.tstat, pstr(p1));
fprintf('     isolated correct vs incorrect: t(%d) = %.2f, %s\n', s2.df, s2.tstat, pstr(p2));
end

function sm = cellmeans(Tc, Ti, accvar)
sids = unique([Tc.subj_id; Ti.subj_id]); n = numel(sids); sm = nan(n,4);
for s = 1:n
    sm(s,1)=m(Tc,sids(s),accvar,1); sm(s,2)=m(Tc,sids(s),accvar,0);
    sm(s,3)=m(Ti,sids(s),accvar,1); sm(s,4)=m(Ti,sids(s),accvar,0);
end
end
function v = m(T,sid,av,val)
d = T.reinst_index(T.subj_id==sid & T.(av)==val); if isempty(d), v=NaN; else, v=mean(d,'omitnan'); end
end
function T = drop_subj(T,exc), T(ismember(T.subj_id,exc),:) = []; end
function ba = attach_b2(ba, bb)
ba.b2_correct = nan(height(ba),1);
for i=1:height(ba)
    mm = bb(bb.subj_id==ba.subj_id(i) & bb.tr_1b_b==ba.tr_1b_b(i), :);
    if height(mm)==1, ba.b2_correct(i)=mm.correct; end
end
end
