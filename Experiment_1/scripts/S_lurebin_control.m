% S_lurebin_control.m
% Bin version of the self-selection control. Uses the MST lure bin (l1 = harder /
% more similar, l2 = easier) instead of pixel similarity. DV = proportion of l1
% (harder) pairs per subject x condition x accuracy, analyzed with the same 2x2.
% Uses every B2 trial (no sparsity), so it avoids the underpowered gaze split.
% Requires image_similarity_AB.csv, rm_2x2.m, gg_eps.m, pstr.m on the path.

clear; clc;
base_dir = '..'; res_dir = fullfile(base_dir,'results');

% pair_key (mst_###_l#) -> is_l1  (1 = hard bin l1, 0 = l2)
T  = readtable(fullfile(res_dir,'image_similarity_AB.csv'));
pk = string(T.base_id) + "_" + string(T.bin);
binmap = containers.Map(cellstr(pk), double(string(T.bin)=="l1"));

load(fullfile(base_dir,'data','eye_movement_data','group_eye_movement_combined.mat'),'Mw');
isB2 = strcmp(Mw.task,'2_back') & strcmp(Mw.identity,'B') & strcmp(Mw.goal,'A-B') & ...
       (strcmp(Mw.condition,'compared') | strcmp(Mw.condition,'isolated'));
B2 = unique(Mw(isB2, {'subj_id','trial_id','stim_id','condition','correct'}));
key = regexprep(string(B2.stim_id), '_B_','_'); key = regexprep(key,'\.png$','');
B2.isl1 = nan(height(B2),1);
for i = 1:height(B2), k = char(key(i)); if isKey(binmap,k), B2.isl1(i)=binmap(k); end, end
fprintf('matched lure bin for %d / %d B2 trials\n', sum(~isnan(B2.isl1)), height(B2));

subs = unique(B2.subj_id); cc=nan(numel(subs),1); ci=cc; ic=cc; ii=cc;
for s = 1:numel(subs)
    d = B2(B2.subj_id==subs(s) & ~isnan(B2.isl1), :);
    cc(s) = mean(d.isl1(strcmp(d.condition,'compared') & d.correct==1), 'omitnan');
    ci(s) = mean(d.isl1(strcmp(d.condition,'compared') & d.correct==0), 'omitnan');
    ic(s) = mean(d.isl1(strcmp(d.condition,'isolated') & d.correct==1), 'omitnan');
    ii(s) = mean(d.isl1(strcmp(d.condition,'isolated') & d.correct==0), 'omitnan');
end

fprintf('\n=== Lure bin (proportion of l1/hard pairs): 2x2 (condition x B2 accuracy) ===\n');
rm_2x2('proportion l1 (hard)', cc, ci, ic, ii);
[~,p1,~,s1] = ttest(cc, ci, 'Tail','left');   % correct < incorrect (fewer hard pairs when correct)
[~,p2,~,s2] = ttest(ic, ii, 'Tail','left');
fprintf('  compared correct vs incorrect: t(%d) = %.2f, p = %.3f  [M_correct = %.3f, M_incorrect = %.3f]\n', s1.df, s1.tstat, p1, mean(cc,'omitnan'), mean(ci,'omitnan'));
fprintf('  isolated correct vs incorrect: t(%d) = %.2f, p = %.3f  [M_correct = %.3f, M_incorrect = %.3f]\n', s2.df, s2.tstat, p2, mean(ic,'omitnan'), mean(ii,'omitnan'));

comp_all = mean([cc ci],2,'omitnan'); iso_all = mean([ic ii],2,'omitnan');
v = ~isnan(comp_all) & ~isnan(iso_all);
[~,pc,~,sc] = ttest(comp_all(v), iso_all(v));
fprintf('  compared vs isolated (overall prop l1): t(%d) = %.2f, p = %.3f  [M_comp = %.3f, M_iso = %.3f]\n', sc.df, sc.tstat, pc, mean(comp_all(v)), mean(iso_all(v)));
