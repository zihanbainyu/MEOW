clear; clc; close all;

base_dir = '..';
res_dir = fullfile(base_dir, 'results');

load(fullfile(res_dir, 'group_eye_movement_combined.mat'), 'Mw');

spatial_params = struct('xres',1920,'yres',1080,'roi_x',[760,1160],'roi_y',[340,740],'grid_x',20,'grid_y',20,'sigma',2);

comp_idx = strcmp(Mw.condition, 'compared');
iso_idx = strcmp(Mw.condition, 'isolated');
ab_idx = strcmp(Mw.goal, 'A-B');
a_idx = strcmp(Mw.identity, 'A');
b_idx = strcmp(Mw.identity, 'B');
task_1b_idx = strcmp(Mw.task, '1_back');
task_2b_idx = strcmp(Mw.task, '2_back');

tr_one_b_a_comp = unique(Mw(task_1b_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_a_iso = unique(Mw(task_1b_idx & a_idx & iso_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_b_comp = unique(Mw(task_1b_idx & b_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_b_iso = unique(Mw(task_1b_idx & b_idx & iso_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_a_comp = unique(Mw(task_2b_idx & ab_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_a_iso = unique(Mw(task_2b_idx & ab_idx & a_idx & iso_idx, {'subj_id','trial_id','stim_id'}));

pairs_aa_comp = innerjoin(tr_one_b_a_comp, tr_two_b_a_comp, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_aa_comp.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id','tr_2b_a'};

pairs_aa_iso = innerjoin(tr_one_b_a_iso, tr_two_b_a_iso, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_aa_iso.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id','tr_2b_a'};

fprintf('Matched A-A pairs: compared=%d, isolated=%d\n', height(pairs_aa_comp), height(pairs_aa_iso));

%% Add B-trial correctness to A-A pairs
% Find 2-back B trials to match via base_id
tr_two_b_b_comp = unique(Mw(task_2b_idx & ab_idx & b_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_b_iso  = unique(Mw(task_2b_idx & ab_idx & b_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));

% Compared: match A pairs to their corresponding B trial
pairs_aa_comp.base_id = regexprep(pairs_aa_comp.stim_id, '_A_', '_BASE_');
tmp_b_comp = tr_two_b_b_comp;
tmp_b_comp.base_id = regexprep(tmp_b_comp.stim_id, '_B_', '_BASE_');

pairs_aa_comp.correct = nan(height(pairs_aa_comp), 1);
for i = 1:height(pairs_aa_comp)
    match = tmp_b_comp.subj_id == pairs_aa_comp.subj_id(i) & strcmp(tmp_b_comp.base_id, pairs_aa_comp.base_id{i});
    if any(match)
        tr_b = tmp_b_comp.trial_id(find(match, 1));
        row = find(Mw.subj_id == pairs_aa_comp.subj_id(i) & Mw.trial_id == tr_b & task_2b_idx, 1);
        if ~isempty(row), pairs_aa_comp.correct(i) = Mw.correct(row); end
    end
end

% Isolated: same logic
pairs_aa_iso.base_id = regexprep(pairs_aa_iso.stim_id, '_A_', '_BASE_');
tmp_b_iso = tr_two_b_b_iso;
tmp_b_iso.base_id = regexprep(tmp_b_iso.stim_id, '_B_', '_BASE_');

pairs_aa_iso.correct = nan(height(pairs_aa_iso), 1);
for i = 1:height(pairs_aa_iso)
    match = tmp_b_iso.subj_id == pairs_aa_iso.subj_id(i) & strcmp(tmp_b_iso.base_id, pairs_aa_iso.base_id{i});
    if any(match)
        tr_b = tmp_b_iso.trial_id(find(match, 1));
        row = find(Mw.subj_id == pairs_aa_iso.subj_id(i) & Mw.trial_id == tr_b & task_2b_idx, 1);
        if ~isempty(row), pairs_aa_iso.correct(i) = Mw.correct(row); end
    end
end

fprintf('B-trial correctness matched: comp=%d/%d, iso=%d/%d\n', ...
    sum(~isnan(pairs_aa_comp.correct)), height(pairs_aa_comp), ...
    sum(~isnan(pairs_aa_iso.correct)), height(pairs_aa_iso));

%% Build baseline pools
unique_subjs = unique([pairs_aa_comp.subj_id; pairs_aa_iso.subj_id]);
subj_baseline_b_comp = struct();
subj_baseline_b_iso = struct();

for s_idx = 1:length(unique_subjs)
    sid = unique_subjs(s_idx); skey = sprintf('s%d', sid);
    sc = pairs_aa_comp(pairs_aa_comp.subj_id == sid, :);
    si = pairs_aa_iso(pairs_aa_iso.subj_id == sid, :);
    if height(sc) > 0
        subj_baseline_b_comp.(skey) = tr_one_b_b_comp(tr_one_b_b_comp.subj_id == sid, :);
    end
    if height(si) > 0
        subj_baseline_b_iso.(skey) = tr_one_b_b_iso(tr_one_b_b_iso.subj_id == sid, :);
    end
end

max_fixations = 10;

%% Compute cumulative A-A compared
fprintf('\nComputing cumulative A-A compared (%d pairs)...\n', height(pairs_aa_comp));
cumulative_results_comp = nan(height(pairs_aa_comp), max_fixations + 4);
n_comp = 0;
for i = 1:height(pairs_aa_comp)
    pair = pairs_aa_comp(i, :);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_comp, skey), continue; end

    fix_1b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
    fix_2b = sortrows(Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_a, :), 'onset');
    if height(fix_1b) < 2 || height(fix_2b) < 2, continue; end

    map_1b_match = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
    if sum(map_1b_match(:)) == 0, continue; end

    partner_b = regexprep(pair.stim_id, '_A_', '_B_');
    pool = subj_baseline_b_comp.(skey);
    baseline_pool = pool(~strcmp(pool.stim_id, partner_b), :);
    if height(baseline_pool) < 1, continue; end

    cumulative_reinst = compute_cumulative_reinstatement_baseline(fix_2b, map_1b_match, baseline_pool, Mw, task_1b_idx, pair.subj_id, spatial_params, max_fixations);
    if all(isnan(cumulative_reinst)), continue; end

    n_comp = n_comp + 1;
    cumulative_results_comp(n_comp, :) = [pair.subj_id, pair.tr_1b_a, pair.tr_2b_a, pair.correct, cumulative_reinst];
end
cumulative_results_comp = cumulative_results_comp(1:n_comp, :);

%% Compute cumulative A-A isolated
fprintf('Computing cumulative A-A isolated (%d pairs)...\n', height(pairs_aa_iso));
cumulative_results_iso = nan(height(pairs_aa_iso), max_fixations + 4);
n_iso = 0;
for i = 1:height(pairs_aa_iso)
    pair = pairs_aa_iso(i, :);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_iso, skey), continue; end

    fix_1b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
    fix_2b = sortrows(Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_a, :), 'onset');
    if height(fix_1b) < 2 || height(fix_2b) < 2, continue; end

    map_1b_match = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
    if sum(map_1b_match(:)) == 0, continue; end

    partner_b = regexprep(pair.stim_id, '_A_', '_B_');
    pool = subj_baseline_b_iso.(skey);
    baseline_pool = pool(~strcmp(pool.stim_id, partner_b), :);
    if height(baseline_pool) < 1, continue; end

    cumulative_reinst = compute_cumulative_reinstatement_baseline(fix_2b, map_1b_match, baseline_pool, Mw, task_1b_idx, pair.subj_id, spatial_params, max_fixations);
    if all(isnan(cumulative_reinst)), continue; end

    n_iso = n_iso + 1;
    cumulative_results_iso(n_iso, :) = [pair.subj_id, pair.tr_1b_a, pair.tr_2b_a, pair.correct, cumulative_reinst];
end
cumulative_results_iso = cumulative_results_iso(1:n_iso, :);

fprintf('Cumulative A-A: compared=%d trials, isolated=%d trials\n', n_comp, n_iso);

%% Convert to tables
fixation_names = arrayfun(@(x) sprintf('fix_%d', x), 1:max_fixations, 'UniformOutput', false);
var_names = [{'subj_id', 'tr_1b_a', 'tr_2b_a', 'correct'}, fixation_names];
cumulative_results_comp = array2table(cumulative_results_comp, 'VariableNames', var_names);
cumulative_results_iso = array2table(cumulative_results_iso, 'VariableNames', var_names);

%% Aggregate
fix_cols_comp = cumulative_results_comp{:, 5:end};
fix_cols_iso = cumulative_results_iso{:, 5:end};
n_trials_per_fix = sum(~isnan(fix_cols_comp), 1);
min_trials = 0.2 * max(n_trials_per_fix);
n_fix_to_plot = find(n_trials_per_fix < min_trials, 1, 'first') - 1;
if isempty(n_fix_to_plot) || n_fix_to_plot < 2, n_fix_to_plot = min(5, size(fix_cols_comp, 2)); end

common_subjs = intersect(unique(cumulative_results_comp.subj_id), unique(cumulative_results_iso.subj_id));
subj_traj_comp = zeros(length(common_subjs), n_fix_to_plot);
subj_traj_iso = zeros(length(common_subjs), n_fix_to_plot);

for i = 1:length(common_subjs)
    sid = common_subjs(i);
    idx_c = cumulative_results_comp.subj_id == sid;
    idx_i = cumulative_results_iso.subj_id == sid;
    if sum(idx_c) > 0, subj_traj_comp(i, :) = mean(fix_cols_comp(idx_c, 1:n_fix_to_plot), 1, 'omitnan'); else, subj_traj_comp(i, :) = NaN; end
    if sum(idx_i) > 0, subj_traj_iso(i, :) = mean(fix_cols_iso(idx_i, 1:n_fix_to_plot), 1, 'omitnan'); else, subj_traj_iso(i, :) = NaN; end
end

save(fullfile(res_dir, 'cumu_reinstat_aa_cor.mat'), ...
    'cumulative_results_comp', 'cumulative_results_iso', 'subj_traj_comp', 'subj_traj_iso', ...
    'common_subjs', 'n_fix_to_plot', 'spatial_params');


%% --- Helper functions ---

function cumulative_reinst = compute_cumulative_reinstatement_baseline(fix_2b, map_1b_match, baseline_pool, Mw, task_1b_idx, subj_id, spatial_params, max_fixations)
    n_fixations = min(height(fix_2b), max_fixations);
    cumulative_reinst = nan(1, max_fixations);
    for f = 1:n_fixations
        fix_subset = fix_2b(1:f, :);
        map_cumulative = create_fixation_map(fix_subset.x, fix_subset.y, fix_subset.dur, spatial_params);
        if sum(map_cumulative(:)) == 0, continue; end
        r_match = corr(map_1b_match(:), map_cumulative(:));
        mismatch_scores = nan(height(baseline_pool), 1);
        for j = 1:height(baseline_pool)
            fix_baseline = Mw(task_1b_idx & Mw.subj_id == subj_id & Mw.trial_id == baseline_pool.trial_id(j), :);
            if height(fix_baseline) < 2, continue; end
            map_baseline = create_fixation_map(fix_baseline.x, fix_baseline.y, fix_baseline.dur, spatial_params);
            if sum(map_baseline(:)) > 0
                mismatch_scores(j) = corr(map_baseline(:), map_cumulative(:));
            end
        end
        cumulative_reinst(f) = r_match - mean(mismatch_scores, 'omitnan');
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