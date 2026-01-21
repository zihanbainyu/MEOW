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


reinstatement_results = struct(); 
reinstatement_results.compared = results_comp;
reinstatement_results.isolated = results_iso; 
reinstatement_results.spatial_params = spatial_params;
% save(fullfile(out_dir, 'gaze_reinstatement_results.mat'), 'reinstatement_results');
fprintf('Done.\n');

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