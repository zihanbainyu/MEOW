% S_image_similarity.m
% Self-selection control for the A1-B1 pattern-separation result.
% Question: does OBJECTIVE A-B image similarity (a fixed stimulus property)
% differ by encoding condition or by subsequent B2 accuracy?  If it does not
% show the same compared-specific correct-vs-incorrect pattern that gaze
% similarity showed, the gaze effect reflects how participants viewed the
% items, not how visually similar the stimuli happened to be.
%
% Objective similarity = Pearson correlation of the grayscale pixels of the
% A and B exemplars (precomputed in results/image_similarity_AB.csv).
% Requires rm_2x2.m, gg_eps.m, pstr.m on the path. Run from Experiment_1/scripts.

clear; clc;
base_dir = '..'; res_dir = fullfile(base_dir,'results');

%% 1. objective A-B similarity lookup:  pair_key (mst_###_l#) -> pixel correlation
T = readtable(fullfile(res_dir,'image_similarity_AB.csv'));
pk  = string(T.base_id) + "_" + string(T.bin);          % e.g. "mst_001_l1"
sim = containers.Map(cellstr(pk), T.pix_corr);

%% 2. per-subject B2 discrimination accuracy for each compared/isolated pair
load(fullfile(base_dir,'data','eye_movement_data','group_eye_movement_combined.mat'),'Mw');
isB2 = strcmp(Mw.task,'2_back') & strcmp(Mw.identity,'B') & strcmp(Mw.goal,'A-B') & ...
       (strcmp(Mw.condition,'compared') | strcmp(Mw.condition,'isolated'));
B2 = unique(Mw(isB2, {'subj_id','trial_id','stim_id','condition','correct'}));  % one row / B2 trial

key = regexprep(string(B2.stim_id), '_B_', '_');          % mst_001_B_l1(.png) -> mst_001_l1(.png)
key = regexprep(key, '\.png$', '');
B2.imgsim = nan(height(B2),1);
for i = 1:height(B2)
    k = char(key(i));
    if isKey(sim,k), B2.imgsim(i) = sim(k); end
end
fprintf('matched objective similarity for %d / %d B2 trials\n', sum(~isnan(B2.imgsim)), height(B2));

%% 3. subject x condition x accuracy cell means, then the 2x2 (mirrors the gaze analysis)
subs = unique(B2.subj_id);
cc = nan(numel(subs),1); ci = cc; ic = cc; ii = cc;
for s = 1:numel(subs)
    d = B2(B2.subj_id==subs(s) & ~isnan(B2.imgsim), :);
    cc(s) = mean(d.imgsim(strcmp(d.condition,'compared') & d.correct==1), 'omitnan');
    ci(s) = mean(d.imgsim(strcmp(d.condition,'compared') & d.correct==0), 'omitnan');
    ic(s) = mean(d.imgsim(strcmp(d.condition,'isolated') & d.correct==1), 'omitnan');
    ii(s) = mean(d.imgsim(strcmp(d.condition,'isolated') & d.correct==0), 'omitnan');
end

fprintf('\n=== Objective image similarity: 2x2 (encoding condition x B2 accuracy) ===\n');
rm_2x2('image similarity', cc, ci, ic, ii);

% simple effects (one-tailed, matching the gaze test's direction: correct < incorrect)
[~,p1,~,s1] = ttest(cc, ci, 'Tail','left');   % compared: correct vs incorrect
[~,p2,~,s2] = ttest(ic, ii, 'Tail','left');   % isolated: correct vs incorrect
fprintf('  compared  correct vs incorrect: t(%d) = %.2f, p = %.3f\n', s1.df, s1.tstat, p1);
fprintf('  isolated  correct vs incorrect: t(%d) = %.2f, p = %.3f\n', s2.df, s2.tstat, p2);

% condition-level check: does objective similarity differ compared vs isolated?
comp_all = mean([cc ci],2,'omitnan'); iso_all = mean([ic ii],2,'omitnan');
v = ~isnan(comp_all) & ~isnan(iso_all);
[~,pc,~,sc] = ttest(comp_all(v), iso_all(v));
fprintf('  compared vs isolated (overall):  t(%d) = %.2f, p = %.3f  [M_comp = %.3f, M_iso = %.3f]\n', ...
        sc.df, sc.tstat, pc, mean(comp_all(v)), mean(iso_all(v)));
