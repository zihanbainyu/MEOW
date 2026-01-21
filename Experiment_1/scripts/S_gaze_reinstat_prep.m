clear; clc; close all;

sub_list = [501, 601, 602, 603, 604, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 617];
tasks = [1, 2]; task_lbls = {'1_back', '2_back'};
base_dir = '..'; out_dir = fullfile(base_dir, 'data', 'eye_movement_data'); behav_dir = fullfile(base_dir, 'data');

% [all_fix, all_sac, all_blk, all_sum, all_behav] = deal({});
% 
% for s_idx = 1:length(sub_list)
%     sid = sub_list(s_idx); fprintf('subj %d\n', sid);
%     s_fldr = sprintf('sub%03d', sid); f_path = fullfile(behav_dir, s_fldr);
% 
%     mat_file = fullfile(f_path, sprintf('sub%03d_concat.mat', sid));
%     if isfile(mat_file)
%         load(mat_file, 'final_data_output');
%         r1 = final_data_output.results_1_back_all;
%         r1.task = repmat({'1_back'}, height(r1), 1); r1.subj_id = repmat(sid, height(r1), 1); r1.trial_id = (1:height(r1))';
%         r2 = final_data_output.results_2_back_all;
%         r2.task = repmat({'2_back'}, height(r2), 1); r2.subj_id = repmat(sid, height(r2), 1); r2.trial_id = (1:height(r2))';
% 
%         f1 = fieldnames(r1); f2 = fieldnames(r2); all_f = unique([f1; f2]);
%         for i = 1:length(all_f)
%             if ~ismember(all_f{i}, f1), r1.(all_f{i}) = repmat({''}, height(r1), 1); end
%             if ~ismember(all_f{i}, f2), r2.(all_f{i}) = repmat({''}, height(r2), 1); end
%         end
%         all_behav = [all_behav; {r1}; {r2}];
%     end
% 
%     acc_tr_t1 = 0; acc_tr_t2 = 0;
%     for b = 1:4
%         for t_idx = 1:length(tasks)
%             tsk = tasks(t_idx); t_lbl = task_lbls{t_idx};
%             fn = fullfile(f_path, sprintf('%03d_%01d_%01d.asc', sid, tsk, b));
%             if ~isfile(fn), continue; end
% 
%             fid = fopen(fn); [f_dat, s_dat, b_dat] = deal(cell(2000,1));
%             [n_f, n_s, n_b, blk_max_raw_tr] = deal(0); [curr_tr, st_t, st_id] = deal(0, -1, 'NA');
% 
%             while ~feof(fid)
%                 ln = fgetl(fid);
%                 if startsWith(ln, 'MSG')
%                     if contains(ln, 'TRIALID')
%                         raw_tr = sscanf(ln, 'MSG %*d TRIALID %d');
%                         if raw_tr > 0
%                             blk_max_raw_tr = max(blk_max_raw_tr, raw_tr);
%                             curr_tr = raw_tr + (tsk==1)*acc_tr_t1 + (tsk==2)*acc_tr_t2;
%                             all_sum{end+1,1} = {sid, b, t_lbl, curr_tr, 'NA', NaN, 0, 0, 0}; st_t = -1;
%                         end
%                     elseif contains(ln, 'STIM_ONSET') && curr_tr > 0
%                         tmp = textscan(ln, '%*s %f %*s %s'); st_t = tmp{1}; st_id = tmp{2}{1};
%                         all_sum{end,1}{5} = st_id; all_sum{end,1}{6} = st_t;
%                     end
%                 elseif curr_tr > 0 && st_t > 0
%                     if startsWith(ln, 'EFIX')
%                         d = sscanf(ln, 'EFIX %c %d %d %d %f %f %d');
%                         if numel(d) == 7
%                             n_f = n_f+1; f_dat{n_f} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7)};
%                             all_sum{end,1}{7} = all_sum{end,1}{7} + 1;
%                         end
%                     elseif startsWith(ln, 'ESACC')
%                         d = sscanf(strrep(ln, '.', 'NaN'), 'ESACC %c %d %d %d %f %f %f %f');
%                         if numel(d) == 8
%                             n_s = n_s+1; s_dat{n_s} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7), d(8)};
%                             all_sum{end,1}{8} = all_sum{end,1}{8} + 1;
%                         end
%                     elseif startsWith(ln, 'EBLINK')
%                         d = sscanf(ln, 'EBLINK %c %d %d %d');
%                         if numel(d) == 4
%                             n_b = n_b+1; b_dat{n_b} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4)};
%                             all_sum{end,1}{9} = all_sum{end,1}{9} + 1;
%                         end
%                     end
%                 end
%             end
%             fclose(fid);
% 
%             if tsk == 1, acc_tr_t1 = acc_tr_t1 + blk_max_raw_tr; else, acc_tr_t2 = acc_tr_t2 + blk_max_raw_tr; end
%             if n_f > 0, all_fix = [all_fix; f_dat(1:n_f)]; end
%             if n_s > 0, all_sac = [all_sac; s_dat(1:n_s)]; end
%             if n_b > 0, all_blk = [all_blk; b_dat(1:n_b)]; end
%         end
%     end
% end
% 
% v_fix = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','x','y','pupil'};
% v_sac = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','sx','sy','ex','ey'};
% v_blk = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur'};
% v_sum = {'subj_id','block','task','trial_id','stim_id','onset_time','n_fix','n_sac','n_blink'};
% 
% fix = cell2table(vertcat(all_fix{:}), 'VariableNames', v_fix);
% sac = cell2table(vertcat(all_sac{:}), 'VariableNames', v_sac);
% blk = cell2table(vertcat(all_blk{:}), 'VariableNames', v_blk);
% sum_tab = cell2table(vertcat(all_sum{:}), 'VariableNames', v_sum);
% 
% all_behav = vertcat(all_behav{:});
% M = outerjoin(fix, all_behav, 'Keys', {'subj_id','task','trial_id','block'}, 'MergeKeys', true, 'Type', 'inner');
% M = removevars(M, 'stim_id_all_behav'); M = renamevars(M, 'stim_id_fix', 'stim_id');
% 
% eye_data = struct(); eye_data.fixations = fix; eye_data.saccades = sac; eye_data.blinks = blk; eye_data.summary = sum_tab;
% eye_data.merged = M; eye_data.merged.resp_key = cellstr(eye_data.merged.resp_key);
% eye_data.merged.resp_key(strcmp(eye_data.merged.resp_key,'NA')) = {'none'};
% eye_data.merged.correct = strcmp(cellstr(eye_data.merged.corr_resp), eye_data.merged.resp_key);
% save(fullfile(out_dir, 'group_eye_movement.mat'), 'eye_data', '-v7.3');
% fprintf('Eye data saved.\n');


%% 
load(fullfile(out_dir, 'group_eye_movement.mat'));
Mw = eye_data.merged(~contains(eye_data.merged.goal, "JUNK"), :);
comp_idx = strcmp(Mw.condition, 'compared'); 
iso_idx = strcmp(Mw.condition, 'isolated');
ab_idx = strcmp(Mw.goal, 'A-B'); 
b_idx = strcmp(Mw.identity, 'B');
task_1b_idx = strcmp(Mw.task, '1_back'); 
task_2b_idx = strcmp(Mw.task, '2_back');

one_b_comp_idx = task_1b_idx & b_idx & comp_idx; 
one_b_iso_idx = task_1b_idx & b_idx & iso_idx;
two_b_comp_idx = task_2b_idx & ab_idx & b_idx & comp_idx; 
two_b_iso_idx = task_2b_idx & ab_idx & b_idx & iso_idx;

tr_one_b_comp = unique(Mw(one_b_comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_iso = unique(Mw(one_b_iso_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_comp = unique(Mw(two_b_comp_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_iso = unique(Mw(two_b_iso_idx, {'subj_id','trial_id','stim_id'}));

pairs_comp = innerjoin(tr_one_b_comp, tr_two_b_comp, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_comp.Properties.VariableNames = {'subj_id','tr_one_b','stim_id','tr_two_b'};

pairs_iso = innerjoin(tr_one_b_iso, tr_two_b_iso, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_iso.Properties.VariableNames = {'subj_id','tr_one_b','stim_id','tr_two_b'};

fprintf('Valid fixations: %d\n', height(Mw));
fprintf('1B: comp=%d, iso=%d trials\n', height(tr_one_b_comp), height(tr_one_b_iso));
fprintf('2B: comp=%d, iso=%d trials\n', height(tr_two_b_comp), height(tr_two_b_iso));
fprintf('matched pairs: comp=%d, iso=%d\n\n', height(pairs_comp), height(pairs_iso));

spatial_params.xres = 1920; spatial_params.yres = 1080;
spatial_params.roi_x = [760, 1160]; spatial_params.roi_y = [340, 740];
spatial_params.grid_x = 20; spatial_params.grid_y = 20;
spatial_params.sigma = 2;
n_baseline = 20;

all_1b_comp_trials = unique(Mw(one_b_comp_idx, {'subj_id','trial_id','stim_id'}));
all_1b_iso_trials = unique(Mw(one_b_iso_idx, {'subj_id','trial_id','stim_id'}));

subj_baseline_comp = struct();
subj_baseline_iso = struct();
unique_subjs = unique(pairs_comp.subj_id);

for s_idx = 1:length(unique_subjs)
    sid = unique_subjs(s_idx);
    pool_comp = all_1b_comp_trials(all_1b_comp_trials.subj_id==sid, :);
    pool_iso = all_1b_iso_trials(all_1b_iso_trials.subj_id==sid, :);
    
    if height(pool_comp) >= n_baseline
        rng(sid); 
        subj_baseline_comp.(sprintf('s%d', sid)) = pool_comp.trial_id(randperm(height(pool_comp), n_baseline));
    end
    
    if height(pool_iso) >= n_baseline
        rng(sid+10000);
        subj_baseline_iso.(sprintf('s%d', sid)) = pool_iso.trial_id(randperm(height(pool_iso), n_baseline));
    end
end

results_comp = []; results_iso = [];

for i = 1:height(pairs_comp)
    pair = pairs_comp(i,:);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_comp, skey), continue; end
    
    baseline_trials = subj_baseline_comp.(skey);
    baseline_trials = baseline_trials(baseline_trials ~= pair.tr_one_b);
    if length(baseline_trials) < n_baseline, continue; end
    baseline_trials = baseline_trials(1:n_baseline);
    
    fix_1b_match = Mw(strcmp(Mw.task, '1_back') & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_one_b, :);
    fix_2b = Mw(strcmp(Mw.task, '2_back') & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_two_b, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    match_score = corr(map_1b_match(:), map_2b(:));
    
    mismatch_scores = nan(n_baseline, 1);
    for j = 1:n_baseline
        tr_mismatch = baseline_trials(j);
        fix_1b_mismatch = Mw(strcmp(Mw.task, '1_back') & Mw.subj_id==pair.subj_id & Mw.trial_id==tr_mismatch, :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    
    baseline_score = mean(mismatch_scores, 'omitnan');
    reinst_index = match_score - baseline_score;
    
    results_comp = [results_comp; pair.subj_id, pair.tr_one_b, pair.tr_two_b, match_score, baseline_score, reinst_index];
end

for i = 1:height(pairs_iso)
    pair = pairs_iso(i,:);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_iso, skey), continue; end
    
    baseline_trials = subj_baseline_iso.(skey);
    baseline_trials = baseline_trials(baseline_trials ~= pair.tr_one_b);
    if length(baseline_trials) < n_baseline, continue; end
    baseline_trials = baseline_trials(1:n_baseline);
    
    fix_1b_match = Mw(strcmp(Mw.task, '1_back') & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_one_b, :);
    fix_2b = Mw(strcmp(Mw.task, '2_back') & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_two_b, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    match_score = corr(map_1b_match(:), map_2b(:));
    
    mismatch_scores = nan(n_baseline, 1);
    for j = 1:n_baseline
        tr_mismatch = baseline_trials(j);
        fix_1b_mismatch = Mw(strcmp(Mw.task, '1_back') & Mw.subj_id==pair.subj_id & Mw.trial_id==tr_mismatch, :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    
    baseline_score = mean(mismatch_scores, 'omitnan');
    reinst_index = match_score - baseline_score;
    
    results_iso = [results_iso; pair.subj_id, pair.tr_one_b, pair.tr_two_b, match_score, baseline_score, reinst_index];
end

results_comp = array2table(results_comp, 'VariableNames', {'subj_id','tr_1b','tr_2b','match_score','baseline_score','reinst_index'});
results_iso = array2table(results_iso, 'VariableNames', {'subj_id','tr_1b','tr_2b','match_score','baseline_score','reinst_index'});

reinstat_res = struct(); 
reinstat_res.compared = results_comp;
reinstat_res.isolated = results_iso; 
reinstat_res.spatial_params = spatial_params;
save(fullfile(out_dir, 'gaze_reinstat_res.mat'), 'reinstat_res');
fprintf('Done.\n');

fprintf('\n=== REINSTATEMENT RESULTS ===\n');
fprintf('Compared (n=%d):\n', height(results_comp));
fprintf('  Match: M=%.3f (SD=%.3f)\n', mean(results_comp.match_score, 'omitnan'), std(results_comp.match_score, 'omitnan'));
fprintf('  Baseline: M=%.3f (SD=%.3f)\n', mean(results_comp.baseline_score, 'omitnan'), std(results_comp.baseline_score, 'omitnan'));
fprintf('  Reinst Index: M=%.3f (SD=%.3f)\n', mean(results_comp.reinst_index, 'omitnan'), std(results_comp.reinst_index, 'omitnan'));
[~,p_c,~,s_c] = ttest(results_comp.match_score, results_comp.baseline_score);
fprintf('  Match vs Baseline: t(%d)=%.2f, p=%.4f%s\n', s_c.df, s_c.tstat, p_c, repmat('*',1,(p_c<0.05)+(p_c<0.01)+(p_c<0.001)));

fprintf('\nIsolated (n=%d):\n', height(results_iso));
fprintf('  Match: M=%.3f (SD=%.3f)\n', mean(results_iso.match_score, 'omitnan'), std(results_iso.match_score, 'omitnan'));
fprintf('  Baseline: M=%.3f (SD=%.3f)\n', mean(results_iso.baseline_score, 'omitnan'), std(results_iso.baseline_score, 'omitnan'));
fprintf('  Reinst Index: M=%.3f (SD=%.3f)\n', mean(results_iso.reinst_index, 'omitnan'), std(results_iso.reinst_index, 'omitnan'));
[~,p_i,~,s_i] = ttest(results_iso.match_score, results_iso.baseline_score);
fprintf('  Match vs Baseline: t(%d)=%.2f, p=%.4f%s\n', s_i.df, s_i.tstat, p_i, repmat('*',1,(p_i<0.05)+(p_i<0.01)+(p_i<0.001)));
fprintf('\nCompared vs Isolated Reinst Index:\n');
[~,p_ci,~,s_ci] = ttest2(results_comp.reinst_index, results_iso.reinst_index);
fprintf('  t(%d)=%.2f, p=%.4f%s\n', s_ci.df, s_ci.tstat, p_ci, repmat('*',1,(p_ci<0.05)+(p_ci<0.01)+(p_ci<0.001)));

n_perm = 1000;

fprintf('\n=== PERMUTATION TESTS ===\n');

% Test 1: Compared reinstatement exists
diffs_c = results_comp.match_score - results_comp.baseline_score;
diffs_c = diffs_c(~isnan(diffs_c));
obs_c = mean(diffs_c);
null_c = zeros(n_perm,1);
for p = 1:n_perm
    signs = (rand(size(diffs_c)) > 0.5)*2 - 1;
    null_c(p) = mean(diffs_c .* signs);
end
p_c = mean(abs(null_c) >= abs(obs_c));
fprintf('Compared reinst > 0: M=%.3f, p=%.4f%s\n', obs_c, p_c, repmat('*',1,(p_c<0.05)+(p_c<0.01)+(p_c<0.001)));

% Test 2: Isolated reinstatement exists
diffs_i = results_iso.match_score - results_iso.baseline_score;
diffs_i = diffs_i(~isnan(diffs_i));
obs_i = mean(diffs_i);
null_i = zeros(n_perm,1);
for p = 1:n_perm
    signs = (rand(size(diffs_i)) > 0.5)*2 - 1;
    null_i(p) = mean(diffs_i .* signs);
end
p_i = mean(abs(null_i) >= abs(obs_i));
fprintf('Isolated reinst > 0: M=%.3f, p=%.4f%s\n', obs_i, p_i, repmat('*',1,(p_i<0.05)+(p_i<0.01)+(p_i<0.001)));

% Test 3: Compared vs Isolated
obs_diff = obs_c - obs_i;
pooled = [diffs_c; diffs_i];
n_c = length(diffs_c); n_i = length(diffs_i);
null_diff = zeros(n_perm,1);
for p = 1:n_perm
    perm_idx = randperm(length(pooled));
    null_diff(p) = mean(pooled(perm_idx(1:n_c))) - mean(pooled(perm_idx(n_c+1:end)));
end
p_diff = mean(abs(null_diff) >= abs(obs_diff));
fprintf('Compared > Isolated: diff=%.3f, p=%.4f%s\n', obs_diff, p_diff, repmat('*',1,(p_diff<0.05)+(p_diff<0.01)+(p_diff<0.001)));

% Merge behavioral accuracy
results_comp_behav = results_comp;
results_comp_behav.correct = nan(height(results_comp_behav), 1);
for i = 1:height(results_comp_behav)
    trial_data = Mw(Mw.subj_id == results_comp_behav.subj_id(i) & ...
                    Mw.trial_id == results_comp_behav.tr_2b(i) & ...
                    strcmp(Mw.task, '2_back'), :);
    if height(trial_data) > 0
        results_comp_behav.correct(i) = trial_data.correct(1);
    end
end
results_comp_behav = results_comp_behav(~isnan(results_comp_behav.correct), :);

results_iso_behav = results_iso;
results_iso_behav.correct = nan(height(results_iso_behav), 1);
for i = 1:height(results_iso_behav)
    trial_data = Mw(Mw.subj_id == results_iso_behav.subj_id(i) & ...
                    Mw.trial_id == results_iso_behav.tr_2b(i) & ...
                    strcmp(Mw.task, '2_back'), :);
    if height(trial_data) > 0
        results_iso_behav.correct(i) = trial_data.correct(1);
    end
end
results_iso_behav = results_iso_behav(~isnan(results_iso_behav.correct), :);

fprintf('\n=== CORRECT VS INCORRECT ===\n');

% Compared
corr_c = results_comp_behav(results_comp_behav.correct==1, :);
incorr_c = results_comp_behav(results_comp_behav.correct==0, :);
fprintf('Compared:\n');
fprintf('  Correct (n=%d): M=%.3f (SD=%.3f)\n', height(corr_c), mean(corr_c.reinst_index,'omitnan'), std(corr_c.reinst_index,'omitnan'));
fprintf('  Incorrect (n=%d): M=%.3f (SD=%.3f)\n', height(incorr_c), mean(incorr_c.reinst_index,'omitnan'), std(incorr_c.reinst_index,'omitnan'));
[~,p_cc,~,s_cc] = ttest2(corr_c.reinst_index, incorr_c.reinst_index);
fprintf('  Correct > Incorrect: t(%d)=%.2f, p=%.4f%s\n', s_cc.df, s_cc.tstat, p_cc, repmat('*',1,(p_cc<0.05)+(p_cc<0.01)+(p_cc<0.001)));

% Isolated
corr_i = results_iso_behav(results_iso_behav.correct==1, :);
incorr_i = results_iso_behav(results_iso_behav.correct==0, :);
fprintf('Isolated:\n');
fprintf('  Correct (n=%d): M=%.3f (SD=%.3f)\n', height(corr_i), mean(corr_i.reinst_index,'omitnan'), std(corr_i.reinst_index,'omitnan'));
fprintf('  Incorrect (n=%d): M=%.3f (SD=%.3f)\n', height(incorr_i), mean(incorr_i.reinst_index,'omitnan'), std(incorr_i.reinst_index,'omitnan'));
[~,p_ci,~,s_ci] = ttest2(corr_i.reinst_index, incorr_i.reinst_index);
fprintf('  Correct > Incorrect: t(%d)=%.2f, p=%.4f%s\n', s_ci.df, s_ci.tstat, p_ci, repmat('*',1,(p_ci<0.05)+(p_ci<0.01)+(p_ci<0.001)));

% Visualization
figure('Position', [100 100 1000 400]);

subplot(1,2,1);
hold on;
x_corr = 1 + (rand(height(corr_c),1)-0.5)*0.2;
x_incorr = 2 + (rand(height(incorr_c),1)-0.5)*0.2;
scatter(x_corr, corr_c.reinst_index, 30, [0.3 0.7 0.3], 'filled', 'MarkerFaceAlpha', 0.4);
scatter(x_incorr, incorr_c.reinst_index, 30, [0.8 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.4);
boxplot([corr_c.reinst_index; incorr_c.reinst_index], [ones(height(corr_c),1); 2*ones(height(incorr_c),1)], ...
    'Labels', {'Correct','Incorrect'}, 'Colors', 'k', 'Symbol', '');
ylabel('Reinstatement Index', 'FontSize', 12, 'FontWeight', 'bold');
title('Compared', 'FontSize', 14, 'FontWeight', 'bold');
ylim([-0.5 1]); grid on; yline(0, 'k--', 'LineWidth', 1); hold off;

subplot(1,2,2);
hold on;
x_corr = 1 + (rand(height(corr_i),1)-0.5)*0.2;
x_incorr = 2 + (rand(height(incorr_i),1)-0.5)*0.2;
scatter(x_corr, corr_i.reinst_index, 30, [0.3 0.7 0.3], 'filled', 'MarkerFaceAlpha', 0.4);
scatter(x_incorr, incorr_i.reinst_index, 30, [0.8 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.4);
boxplot([corr_i.reinst_index; incorr_i.reinst_index], [ones(height(corr_i),1); 2*ones(height(incorr_i),1)], ...
    'Labels', {'Correct','Incorrect'}, 'Colors', 'k', 'Symbol', '');
ylabel('Reinstatement Index', 'FontSize', 12, 'FontWeight', 'bold');
title('Isolated', 'FontSize', 14, 'FontWeight', 'bold');
ylim([-0.5 1]); grid on; yline(0, 'k--', 'LineWidth', 1); hold off;


function map = create_fixation_map(x, y, dur, params)
    in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & ...
             y >= params.roi_y(1) & y <= params.roi_y(2);
    x = x(in_roi); y = y(in_roi); dur = dur(in_roi);
    
    if isempty(x)
        map = zeros(params.grid_y, params.grid_x);
        return;
    end
    
    x_bins = linspace(params.roi_x(1), params.roi_x(2), params.grid_x+1);
    y_bins = linspace(params.roi_y(1), params.roi_y(2), params.grid_y+1);
    map = zeros(params.grid_y, params.grid_x);
    
    for i = 1:length(x)
        x_idx = find(x(i) >= x_bins(1:end-1) & x(i) < x_bins(2:end), 1);
        y_idx = find(y(i) >= y_bins(1:end-1) & y(i) < y_bins(2:end), 1);
        if ~isempty(x_idx) && ~isempty(y_idx)
            map(y_idx, x_idx) = map(y_idx, x_idx) + dur(i);
        end
    end
    
    if params.sigma > 0
        map = imgaussfilt(map, params.sigma);
    end
    if sum(map(:)) > 0
        map = map / sum(map(:));
    end
end