clear; clc; close all;

base_dir = '..'; 
res_dir = fullfile(base_dir, 'results');
out_dir = fullfile(res_dir, 'gaze_entropy_temporal');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

load(fullfile(base_dir, 'data', 'eye_movement_data', 'group_eye_movement_combined.mat'), 'Mw');

spatial_params = struct('xres', 1920, 'yres', 1080, ...
    'roi_x', [760, 1160], 'roi_y', [340, 740], ...
    'grid_x', 20, 'grid_y', 20, 'sigma', 2);

comp_idx = strcmp(Mw.condition, 'compared'); 
iso_idx = strcmp(Mw.condition, 'isolated');
ab_idx = strcmp(Mw.goal, 'A-B'); 
a_idx = strcmp(Mw.identity, 'A');
b_idx = strcmp(Mw.identity, 'B');
task_1b_idx = strcmp(Mw.task, '1_back'); 
task_2b_idx = strcmp(Mw.task, '2_back');

one_b_a_comp_idx = task_1b_idx & a_idx & comp_idx;
one_b_a_iso_idx = task_1b_idx & a_idx & iso_idx;
one_b_b_comp_idx = task_1b_idx & b_idx & comp_idx; 
one_b_b_iso_idx = task_1b_idx & b_idx & iso_idx;
two_b_a_comp_idx = task_2b_idx & ab_idx & a_idx & comp_idx;
two_b_a_iso_idx = task_2b_idx & ab_idx & a_idx & iso_idx;
two_b_b_comp_idx = task_2b_idx & ab_idx & b_idx & comp_idx; 
two_b_b_iso_idx = task_2b_idx & ab_idx & b_idx & iso_idx;

tr_one_b_a_comp = unique(Mw(one_b_a_comp_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_one_b_a_iso = unique(Mw(one_b_a_iso_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_one_b_b_comp = unique(Mw(one_b_b_comp_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_one_b_b_iso = unique(Mw(one_b_b_iso_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_two_b_a_comp = unique(Mw(two_b_a_comp_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_two_b_a_iso = unique(Mw(two_b_a_iso_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_two_b_b_comp = unique(Mw(two_b_b_comp_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_two_b_b_iso = unique(Mw(two_b_b_iso_idx, {'subj_id', 'trial_id', 'stim_id'}));

pairs_bb_comp = innerjoin(tr_one_b_b_comp, tr_two_b_b_comp, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_bb_comp.Properties.VariableNames = {'subj_id','tr_1b_b','stim_id','tr_2b_b'};

pairs_bb_iso = innerjoin(tr_one_b_b_iso, tr_two_b_b_iso, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_bb_iso.Properties.VariableNames = {'subj_id','tr_1b_b','stim_id','tr_2b_b'};

fprintf('Matched B-B pairs: compared=%d, isolated=%d\n', height(pairs_bb_comp), height(pairs_bb_iso));


fprintf('\n=== A item Entropy (1-back encoding) ===\n');

subjs_a_comp = unique(tr_one_b_a_comp.subj_id);
subjs_a_iso = unique(tr_one_b_a_iso.subj_id);
common_subjs_a = intersect(subjs_a_comp, subjs_a_iso);

entropy_a_comp_subj = zeros(length(common_subjs_a), 1);
entropy_a_iso_subj = zeros(length(common_subjs_a), 1);

for i = 1:length(common_subjs_a)
    sid = common_subjs_a(i);
    trials_comp = tr_one_b_a_comp(tr_one_b_a_comp.subj_id == sid, :);
    trials_iso = tr_one_b_a_iso(tr_one_b_a_iso.subj_id == sid, :);

    subj_entropy_comp = [];
    for j = 1:height(trials_comp)
        fix_a = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == trials_comp.trial_id(j), :);
        if height(fix_a) >= 2
            map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
            ent = compute_spatial_entropy(map_a);
            if ent > 0, subj_entropy_comp = [subj_entropy_comp; ent]; end
        end
    end

    subj_entropy_iso = [];
    for j = 1:height(trials_iso)
        fix_a = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == trials_iso.trial_id(j), :);
        if height(fix_a) >= 2
            map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
            ent = compute_spatial_entropy(map_a);
            if ent > 0, subj_entropy_iso = [subj_entropy_iso; ent]; end
        end
    end

    entropy_a_comp_subj(i) = mean(subj_entropy_comp);
    entropy_a_iso_subj(i) = mean(subj_entropy_iso);
end

[~, p_entropy_a, ~, stats_entropy_a] = ttest(entropy_a_comp_subj, entropy_a_iso_subj);
fprintf('Compared: M=%.3f (SD=%.3f)\n', mean(entropy_a_comp_subj), std(entropy_a_comp_subj));
fprintf('Isolated: M=%.3f (SD=%.3f)\n', mean(entropy_a_iso_subj), std(entropy_a_iso_subj));
fprintf('t(%d) = %.3f, p = %.4f\n', stats_entropy_a.df, stats_entropy_a.tstat, p_entropy_a);

fprintf('\n=== B item Entropy (1-back encoding) ===\n');

subjs_b_comp = unique(tr_one_b_b_comp.subj_id);
subjs_b_iso = unique(tr_one_b_b_iso.subj_id);
common_subjs_b = intersect(subjs_b_comp, subjs_b_iso);

entropy_b_comp_subj = zeros(length(common_subjs_b), 1);
entropy_b_iso_subj = zeros(length(common_subjs_b), 1);

for i = 1:length(common_subjs_b)
    sid = common_subjs_b(i);
    trials_comp = tr_one_b_b_comp(tr_one_b_b_comp.subj_id == sid, :);
    trials_iso = tr_one_b_b_iso(tr_one_b_b_iso.subj_id == sid, :);

    subj_entropy_comp = [];
    for j = 1:height(trials_comp)
        fix_b = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == trials_comp.trial_id(j), :);
        if height(fix_b) >= 2
            map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);
            ent = compute_spatial_entropy(map_b);
            if ent > 0, subj_entropy_comp = [subj_entropy_comp; ent]; end
        end
    end

    subj_entropy_iso = [];
    for j = 1:height(trials_iso)
        fix_b = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == trials_iso.trial_id(j), :);
        if height(fix_b) >= 2
            map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);
            ent = compute_spatial_entropy(map_b);
            if ent > 0, subj_entropy_iso = [subj_entropy_iso; ent]; end
        end
    end

    entropy_b_comp_subj(i) = mean(subj_entropy_comp);
    entropy_b_iso_subj(i) = mean(subj_entropy_iso);
end

[~, p_entropy_b, ~, stats_entropy_b] = ttest(entropy_b_comp_subj, entropy_b_iso_subj);
fprintf('Compared: M=%.3f (SD=%.3f)\n', mean(entropy_b_comp_subj), std(entropy_b_comp_subj));
fprintf('Isolated: M=%.3f (SD=%.3f)\n', mean(entropy_b_iso_subj), std(entropy_b_iso_subj));
fprintf('t(%d) = %.3f, p = %.4f\n', stats_entropy_b.df, stats_entropy_b.tstat, p_entropy_b);

fprintf('\n=== A-B Fixation Overlap during 1-back (Encoding) ===\n');

tr_one_b_a_comp.base_id = regexprep(tr_one_b_a_comp.stim_id, '_A_', '_BASE_');
tr_one_b_b_comp_temp = tr_one_b_b_comp;
tr_one_b_b_comp_temp.base_id = regexprep(tr_one_b_b_comp_temp.stim_id, '_B_', '_BASE_');

ab_pairs_comp = innerjoin(tr_one_b_a_comp, tr_one_b_b_comp_temp, ...
    'Keys', {'subj_id', 'base_id'}, ...
    'LeftVariables', {'subj_id', 'trial_id', 'stim_id', 'base_id'}, ...
    'RightVariables', {'trial_id', 'stim_id'});
ab_pairs_comp.Properties.VariableNames = {'subj_id', 'tr_a', 'stim_a', 'base_id', 'tr_b', 'stim_b'};

tr_one_b_a_iso_temp = tr_one_b_a_iso;
tr_one_b_a_iso_temp.base_id = regexprep(tr_one_b_a_iso_temp.stim_id, '_A_', '_BASE_');
tr_one_b_b_iso_temp = tr_one_b_b_iso;
tr_one_b_b_iso_temp.base_id = regexprep(tr_one_b_b_iso_temp.stim_id, '_B_', '_BASE_');

ab_pairs_iso = innerjoin(tr_one_b_a_iso_temp, tr_one_b_b_iso_temp, ...
    'Keys', {'subj_id', 'base_id'}, ...
    'LeftVariables', {'subj_id', 'trial_id', 'stim_id', 'base_id'}, ...
    'RightVariables', {'trial_id', 'stim_id'});
ab_pairs_iso.Properties.VariableNames = {'subj_id', 'tr_a', 'stim_a', 'base_id', 'tr_b', 'stim_b'};

overlap_comp = zeros(height(ab_pairs_comp), 1);
for i = 1:height(ab_pairs_comp)
    pair = ab_pairs_comp(i, :);
    fix_a = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_a, :);
    fix_b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_b, :);

    if height(fix_a) < 2 || height(fix_b) < 2, continue; end

    map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
    map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);

    if sum(map_a(:)) == 0 || sum(map_b(:)) == 0, continue; end

    overlap_comp(i) = corr(map_a(:), map_b(:));
end

overlap_iso = zeros(height(ab_pairs_iso), 1);
for i = 1:height(ab_pairs_iso)
    pair = ab_pairs_iso(i, :);
    fix_a = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_a, :);
    fix_b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_b, :);

    if height(fix_a) < 2 || height(fix_b) < 2, continue; end

    map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
    map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);

    if sum(map_a(:)) == 0 || sum(map_b(:)) == 0, continue; end

    overlap_iso(i) = corr(map_a(:), map_b(:));
end

overlap_comp_valid = overlap_comp(overlap_comp ~= 0);
overlap_iso_valid = overlap_iso(overlap_iso ~= 0);

fprintf('Compared: M=%.3f (SD=%.3f), N=%d pairs\n', ...
    mean(overlap_comp_valid), std(overlap_comp_valid), length(overlap_comp_valid));
fprintf('Isolated: M=%.3f (SD=%.3f), N=%d pairs\n', ...
    mean(overlap_iso_valid), std(overlap_iso_valid), length(overlap_iso_valid));

[~, p_overlap] = ttest2(overlap_comp_valid, overlap_iso_valid);
fprintf('t-test: p = %.4f\n', p_overlap);

fprintf('\n=== A-B Fixation Overlap during 2-back (Retrieval) ===\n');

tr_two_b_a_comp.base_id = regexprep(tr_two_b_a_comp.stim_id, '_A_', '_BASE_');
tr_two_b_b_comp_temp = tr_two_b_b_comp;
tr_two_b_b_comp_temp.base_id = regexprep(tr_two_b_b_comp_temp.stim_id, '_B_', '_BASE_');

ab_pairs_2b_comp = innerjoin(tr_two_b_a_comp, tr_two_b_b_comp_temp, ...
    'Keys', {'subj_id', 'base_id'}, ...
    'LeftVariables', {'subj_id', 'trial_id', 'stim_id', 'base_id'}, ...
    'RightVariables', {'trial_id', 'stim_id'});
ab_pairs_2b_comp.Properties.VariableNames = {'subj_id', 'tr_a', 'stim_a', 'base_id', 'tr_b', 'stim_b'};

tr_two_b_a_iso.base_id = regexprep(tr_two_b_a_iso.stim_id, '_A_', '_BASE_');
tr_two_b_b_iso_temp = tr_two_b_b_iso;
tr_two_b_b_iso_temp.base_id = regexprep(tr_two_b_b_iso_temp.stim_id, '_B_', '_BASE_');

ab_pairs_2b_iso = innerjoin(tr_two_b_a_iso, tr_two_b_b_iso_temp, ...
    'Keys', {'subj_id', 'base_id'}, ...
    'LeftVariables', {'subj_id', 'trial_id', 'stim_id', 'base_id'}, ...
    'RightVariables', {'trial_id', 'stim_id'});
ab_pairs_2b_iso.Properties.VariableNames = {'subj_id', 'tr_a', 'stim_a', 'base_id', 'tr_b', 'stim_b'};

overlap_2b_comp = zeros(height(ab_pairs_2b_comp), 1);
for i = 1:height(ab_pairs_2b_comp)
    pair = ab_pairs_2b_comp(i, :);
    fix_a = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_a, :);
    fix_b = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_b, :);

    if height(fix_a) < 2 || height(fix_b) < 2, continue; end

    map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
    map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);

    if sum(map_a(:)) == 0 || sum(map_b(:)) == 0, continue; end

    overlap_2b_comp(i) = corr(map_a(:), map_b(:));
end

overlap_2b_iso = zeros(height(ab_pairs_2b_iso), 1);
for i = 1:height(ab_pairs_2b_iso)
    pair = ab_pairs_2b_iso(i, :);
    fix_a = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_a, :);
    fix_b = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_b, :);

    if height(fix_a) < 2 || height(fix_b) < 2, continue; end

    map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
    map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);

    if sum(map_a(:)) == 0 || sum(map_b(:)) == 0, continue; end

    overlap_2b_iso(i) = corr(map_a(:), map_b(:));
end

overlap_2b_comp_valid = overlap_2b_comp(overlap_2b_comp ~= 0);
overlap_2b_iso_valid = overlap_2b_iso(overlap_2b_iso ~= 0);


%% Cumulative gaze reinstatement (baseline-corrected)
fprintf('\n=== Computing Baseline-Corrected Cumulative Reinstatement ===\n');

% Build baseline pools (same as main gaze script)
unique_subjs = unique(pairs_bb_comp.subj_id);
subj_baseline_b_comp = struct(); 
subj_baseline_b_iso = struct();

for s_idx = 1:length(unique_subjs)
    sid = unique_subjs(s_idx);
    pairs_comp_subj = pairs_bb_comp(pairs_bb_comp.subj_id == sid, :);
    pairs_iso_subj = pairs_bb_iso(pairs_bb_iso.subj_id == sid, :);
    if height(pairs_comp_subj) > 0
        pool_b_comp = tr_one_b_b_comp(tr_one_b_b_comp.subj_id==sid & ismember(tr_one_b_b_comp.trial_id, pairs_comp_subj.tr_1b_b), :);
        subj_baseline_b_comp.(sprintf('s%d', sid)) = pool_b_comp;
    end
    if height(pairs_iso_subj) > 0
        pool_b_iso = tr_one_b_b_iso(tr_one_b_b_iso.subj_id==sid & ismember(tr_one_b_b_iso.trial_id, pairs_iso_subj.tr_1b_b), :);
        subj_baseline_b_iso.(sprintf('s%d', sid)) = pool_b_iso;
    end
end

max_fixations = 10;
cumulative_results_comp = nan(height(pairs_bb_comp), max_fixations + 3);
cumulative_results_iso = nan(height(pairs_bb_iso), max_fixations + 3);

% Compared condition
n_comp_cum = 0;
for i = 1:height(pairs_bb_comp)
    pair = pairs_bb_comp(i, :);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_comp, skey), continue; end
    
    fix_1b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_b, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_b, :);
    if height(fix_1b) < 2 || height(fix_2b) < 2, continue; end
    fix_2b = sortrows(fix_2b, 'onset');
    
    % Match map (encoding)
    map_1b_match = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
    if sum(map_1b_match(:)) == 0, continue; end
    
    % Baseline pool
    pool = subj_baseline_b_comp.(skey);
    current_base = regexprep(pair.stim_id, '_B_', '_BASE_');
    matching_a_stim = regexprep(current_base, '_BASE_', '_A_');
    baseline_pool = pool(pool.trial_id ~= pair.tr_1b_b & ~strcmp(pool.stim_id, matching_a_stim), :);
    if height(baseline_pool) < 1, continue; end
    
    % Compute baseline-corrected cumulative reinstatement
    cumulative_reinst = compute_cumulative_reinstatement_baseline(fix_2b, map_1b_match, baseline_pool, Mw, task_1b_idx, pair.subj_id, spatial_params, max_fixations);
    if all(isnan(cumulative_reinst)), continue; end
    
    n_comp_cum = n_comp_cum + 1;
    cumulative_results_comp(n_comp_cum, :) = [pair.subj_id, pair.tr_1b_b, pair.tr_2b_b, cumulative_reinst];
end
cumulative_results_comp = cumulative_results_comp(1:n_comp_cum, :);

% Isolated condition
n_iso_cum = 0;
for i = 1:height(pairs_bb_iso)
    pair = pairs_bb_iso(i, :);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_iso, skey), continue; end
    
    fix_1b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_b, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_b, :);
    if height(fix_1b) < 2 || height(fix_2b) < 2, continue; end
    fix_2b = sortrows(fix_2b, 'onset');
    
    map_1b_match = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
    if sum(map_1b_match(:)) == 0, continue; end
    
    pool = subj_baseline_b_iso.(skey);
    current_base = regexprep(pair.stim_id, '_B_', '_BASE_');
    matching_a_stim = regexprep(current_base, '_BASE_', '_A_');
    baseline_pool = pool(pool.trial_id ~= pair.tr_1b_b & ~strcmp(pool.stim_id, matching_a_stim), :);
    if height(baseline_pool) < 1, continue; end
    
    cumulative_reinst = compute_cumulative_reinstatement_baseline(fix_2b, map_1b_match, baseline_pool, Mw, task_1b_idx, pair.subj_id, spatial_params, max_fixations);
    if all(isnan(cumulative_reinst)), continue; end
    
    n_iso_cum = n_iso_cum + 1;
    cumulative_results_iso(n_iso_cum, :) = [pair.subj_id, pair.tr_1b_b, pair.tr_2b_b, cumulative_reinst];
end
cumulative_results_iso = cumulative_results_iso(1:n_iso_cum, :);
fprintf('Cumulative: compared=%d trials, isolated=%d trials\n', n_comp_cum, n_iso_cum);

%% Aggregate by subject and plot
fixation_names = arrayfun(@(x) sprintf('fix_%d', x), 1:max_fixations, 'UniformOutput', false);
var_names = [{'subj_id', 'tr_1b_b', 'tr_2b_b'}, fixation_names];
cumulative_results_comp = array2table(cumulative_results_comp, 'VariableNames', var_names);
cumulative_results_iso = array2table(cumulative_results_iso, 'VariableNames', var_names);

fix_cols_comp = cumulative_results_comp{:, 4:end};
fix_cols_iso = cumulative_results_iso{:, 4:end};
n_trials_per_fix_comp = sum(~isnan(fix_cols_comp), 1);
min_trials = 0.2 * max(n_trials_per_fix_comp);
n_fix_to_plot = find(n_trials_per_fix_comp < min_trials, 1, 'first') - 1;
if isempty(n_fix_to_plot) || n_fix_to_plot < 2, n_fix_to_plot = min(5, size(fix_cols_comp, 2)); end

common_subjs_cum = intersect(unique(cumulative_results_comp.subj_id), unique(cumulative_results_iso.subj_id));
subj_traj_comp = zeros(length(common_subjs_cum), n_fix_to_plot);
subj_traj_iso = zeros(length(common_subjs_cum), n_fix_to_plot);

for i = 1:length(common_subjs_cum)
    sid = common_subjs_cum(i);
    idx_comp = cumulative_results_comp.subj_id == sid;
    idx_iso = cumulative_results_iso.subj_id == sid;
    if sum(idx_comp) > 0, subj_traj_comp(i, :) = mean(fix_cols_comp(idx_comp, 1:n_fix_to_plot), 1, 'omitnan'); else, subj_traj_comp(i, :) = NaN; end
    if sum(idx_iso) > 0, subj_traj_iso(i, :) = mean(fix_cols_iso(idx_iso, 1:n_fix_to_plot), 1, 'omitnan'); else, subj_traj_iso(i, :) = NaN; end
end

mean_comp = mean(subj_traj_comp, 1, 'omitnan');
sem_comp = std(subj_traj_comp, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(subj_traj_comp), 1));
mean_iso = mean(subj_traj_iso, 1, 'omitnan');
sem_iso = std(subj_traj_iso, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(subj_traj_iso), 1));
fix_numbers = 1:n_fix_to_plot;

%% Main trajectory plot
figure('Position', [100, 100, 800, 500]);
hold on;
fill([fix_numbers, fliplr(fix_numbers)], [mean_comp + sem_comp, fliplr(mean_comp - sem_comp)], [0.2, 0.6, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([fix_numbers, fliplr(fix_numbers)], [mean_iso + sem_iso, fliplr(mean_iso - sem_iso)], [0.8, 0.4, 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(fix_numbers, mean_comp, 'o-', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.2, 0.6, 0.8]);
plot(fix_numbers, mean_iso, 's-', 'Color', [0.8, 0.4, 0.4], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.8, 0.4, 0.4]);
plot([0.5, n_fix_to_plot + 0.5], [0, 0], 'k--', 'LineWidth', 1);
xlabel('Cumulative Fixation Number', 'FontSize', 12); ylabel('Gaze Reinstatement (baseline-corrected)', 'FontSize', 12);
title('Fixation-by-Fixation Reinstatement Buildup', 'FontSize', 14);
legend({'', '', 'Compared', 'Isolated'}, 'Location', 'best', 'FontSize', 11);
xlim([0.5, n_fix_to_plot + 0.5]); box off; grid on;
print(gcf, fullfile(res_dir, 'cumulative_reinstatement.pdf'), '-dpdf', '-vector');

% Rest continues as before...% 
% figure('Position', [100, 100, 800, 500]);
% hold on;
% fill([fix_numbers, fliplr(fix_numbers)], [mean_comp + sem_comp, fliplr(mean_comp - sem_comp)], [0.2, 0.6, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% fill([fix_numbers, fliplr(fix_numbers)], [mean_iso + sem_iso, fliplr(mean_iso - sem_iso)], [0.8, 0.4, 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% plot(fix_numbers, mean_comp, 'o-', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.2, 0.6, 0.8]);
% plot(fix_numbers, mean_iso, 's-', 'Color', [0.8, 0.4, 0.4], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.8, 0.4, 0.4]);
% plot([0.5, n_fix_to_plot + 0.5], [0, 0], 'k--', 'LineWidth', 1);
% xlabel('Cumulative Fixation Number'); ylabel('Gaze Reinstatement (r)'); title('Cumulative Fixation-by-Fixation Reinstatement');
% legend({'', '', 'Compared', 'Isolated'}, 'Location', 'best', 'FontSize', 11);
% xlim([0.5, n_fix_to_plot + 0.5]); box off; set(gca, 'FontSize', 11); grid on;
% 
% figure('Position', [100, 100, 900, 400]);
% subplot(1, 2, 1); hold on;
% valid = ~isnan(early_comp) & ~isnan(early_iso);
% bar_data = [mean(early_comp(valid)), mean(early_iso(valid))];
% bar_err = [std(early_comp(valid))/sqrt(sum(valid)), std(early_iso(valid))/sqrt(sum(valid))];
% b = bar(bar_data, 'FaceColor', 'flat'); b.CData(1,:) = [0.2, 0.6, 0.8]; b.CData(2,:) = [0.8, 0.4, 0.4];
% errorbar(1:2, bar_data, bar_err, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
% set(gca, 'XTickLabel', {'Compared', 'Isolated'}); ylabel('Mean Reinstatement (r)');
% title(sprintf('Early (Fix %d-%d)', early_range(1), early_range(end)));
% if p_early < 0.05, y_sig = max(bar_data) + max(bar_err) + 0.05; plot([1, 2], [y_sig, y_sig], 'k-', 'LineWidth', 1); text(1.5, y_sig + 0.02, sprintf('p=%.3f', p_early), 'HorizontalAlignment', 'center', 'FontSize', 10); end
% box off; set(gca, 'FontSize', 11);
% 
% subplot(1, 2, 2); hold on;
% valid = ~isnan(late_comp) & ~isnan(late_iso);
% bar_data = [mean(late_comp(valid)), mean(late_iso(valid))];
% bar_err = [std(late_comp(valid))/sqrt(sum(valid)), std(late_iso(valid))/sqrt(sum(valid))];
% b = bar(bar_data, 'FaceColor', 'flat'); b.CData(1,:) = [0.2, 0.6, 0.8]; b.CData(2,:) = [0.8, 0.4, 0.4];
% errorbar(1:2, bar_data, bar_err, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
% set(gca, 'XTickLabel', {'Compared', 'Isolated'}); ylabel('Mean Reinstatement (r)');
% title(sprintf('Late (Fix %d-%d)', late_range(1), late_range(end)));
% if p_late < 0.05, y_sig = max(bar_data) + max(bar_err) + 0.05; plot([1, 2], [y_sig, y_sig], 'k-', 'LineWidth', 1); text(1.5, y_sig + 0.02, sprintf('p=%.3f', p_late), 'HorizontalAlignment', 'center', 'FontSize', 10); end
% box off; set(gca, 'FontSize', 11);


function cumulative_reinst = compute_cumulative_reinstatement_baseline(fix_2b, map_1b_match, baseline_pool, Mw, task_1b_idx, subj_id, spatial_params, max_fixations)
    n_fixations = min(height(fix_2b), max_fixations);
    cumulative_reinst = nan(1, max_fixations);
    
    for f = 1:n_fixations
        % Cumulative retrieval pattern (fixations 1 to f)
        fix_subset = fix_2b(1:f, :);
        map_cumulative = create_fixation_map(fix_subset.x, fix_subset.y, fix_subset.dur, spatial_params);
        if sum(map_cumulative(:)) == 0, continue; end
        
        % Match score
        r_match = corr(map_1b_match(:), map_cumulative(:));
        
        % Baseline score (ALSO cumulative - use same fixations 1:f from baseline trials)
        mismatch_scores = nan(height(baseline_pool), 1);
        for j = 1:height(baseline_pool)
            fix_baseline = Mw(task_1b_idx & Mw.subj_id == subj_id & Mw.trial_id == baseline_pool.trial_id(j), :);
            if height(fix_baseline) < f, continue; end  % Need at least f fixations
            
            % Use SAME cumulative fixations (1:f) from baseline encoding
            fix_baseline_subset = fix_baseline(1:f, :);
            map_baseline_cumulative = create_fixation_map(fix_baseline_subset.x, fix_baseline_subset.y, fix_baseline_subset.dur, spatial_params);
            
            if sum(map_baseline_cumulative(:)) > 0
                mismatch_scores(j) = corr(map_baseline_cumulative(:), map_cumulative(:));
            end
        end
        r_baseline = mean(mismatch_scores, 'omitnan');
        
        % Baseline-corrected cumulative reinstatement
        cumulative_reinst(f) = r_match - r_baseline;
    end
end
function map = create_fixation_map(x, y, dur, params)
    in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & y >= params.roi_y(1) & y <= params.roi_y(2);
    x = x(in_roi); y = y(in_roi); dur = dur(in_roi);
    if isempty(x), map = zeros(params.grid_y, params.grid_x); return; end
    x_bins = linspace(params.roi_x(1), params.roi_x(2), params.grid_x+1);
    y_bins = linspace(params.roi_y(1), params.roi_y(2), params.grid_y+1);
    map = zeros(params.grid_y, params.grid_x);
    for i = 1:length(x)
        x_idx = find(x(i) >= x_bins(1:end-1) & x(i) < x_bins(2:end), 1);
        y_idx = find(y(i) >= y_bins(1:end-1) & y(i) < y_bins(2:end), 1);
        if ~isempty(x_idx) && ~isempty(y_idx), map(y_idx, x_idx) = map(y_idx, x_idx) + dur(i); end
    end
    if params.sigma > 0, map = imgaussfilt(map, params.sigma); end
    if sum(map(:)) > 0, map = map / sum(map(:)); end
end

function entropy = compute_spatial_entropy(map)
    p = map(:); p = p(p > 0);
    if isempty(p), entropy = 0; else, entropy = -sum(p .* log2(p)); end
end