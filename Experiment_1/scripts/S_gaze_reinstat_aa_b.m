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

%% Build trial tables
tr_one_b_a_comp = unique(Mw(task_1b_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_a_iso  = unique(Mw(task_1b_idx & a_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));
tr_one_b_b_comp = unique(Mw(task_1b_idx & b_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_b_iso  = unique(Mw(task_1b_idx & b_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));
tr_two_b_a_comp = unique(Mw(task_2b_idx & ab_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_a_iso  = unique(Mw(task_2b_idx & ab_idx & a_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));

%% Match A1 encoding to A2 retrieval (same stim_id)
pairs_aa_comp = innerjoin(tr_one_b_a_comp, tr_two_b_a_comp, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_aa_comp.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id','tr_2b_a'};

pairs_aa_iso = innerjoin(tr_one_b_a_iso, tr_two_b_a_iso, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_aa_iso.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id','tr_2b_a'};

fprintf('Matched A-A pairs: comp=%d, iso=%d\n', height(pairs_aa_comp), height(pairs_aa_iso));

%% Add B2 lure discrimination correctness
% For each A-A pair, find the corresponding B trial at 2-back.
% B shares the same base stimulus (e.g., juice_A_01 -> juice_B_01).
% We look up correctness from the B ROW in Mw (identity=='B', task=='2_back').

pairs_aa_comp.correct = nan(height(pairs_aa_comp), 1);
for i = 1:height(pairs_aa_comp)
    sid = pairs_aa_comp.subj_id(i);
    base = regexprep(pairs_aa_comp.stim_id{i}, '_A_', '_BASE_');
    b_rows = Mw(Mw.subj_id == sid & task_2b_idx & b_idx & comp_idx & ab_idx, :);
    b_base = regexprep(b_rows.stim_id, '_B_', '_BASE_');
    match = strcmp(b_base, base);
    if any(match)
        pairs_aa_comp.correct(i) = b_rows.correct(find(match, 1));
    end
end

pairs_aa_iso.correct = nan(height(pairs_aa_iso), 1);
for i = 1:height(pairs_aa_iso)
    sid = pairs_aa_iso.subj_id(i);
    base = regexprep(pairs_aa_iso.stim_id{i}, '_A_', '_BASE_');
    b_rows = Mw(Mw.subj_id == sid & task_2b_idx & b_idx & iso_idx & ab_idx, :);
    b_base = regexprep(b_rows.stim_id, '_B_', '_BASE_');
    match = strcmp(b_base, base);
    if any(match)
        pairs_aa_iso.correct(i) = b_rows.correct(find(match, 1));
    end
end

fprintf('B2 correctness matched: comp=%d/%d, iso=%d/%d\n', ...
    sum(~isnan(pairs_aa_comp.correct)), height(pairs_aa_comp), ...
    sum(~isnan(pairs_aa_iso.correct)), height(pairs_aa_iso));
fprintf('B2 correct rate: comp=%.1f%%, iso=%.1f%%\n', ...
    100*mean(pairs_aa_comp.correct, 'omitnan'), ...
    100*mean(pairs_aa_iso.correct, 'omitnan'));

%% Build B baseline pools (other B1 encoding trials, same condition)
unique_subjs = unique([pairs_aa_comp.subj_id; pairs_aa_iso.subj_id]);
subj_baseline_b_comp = struct();
subj_baseline_b_iso  = struct();

for s_idx = 1:length(unique_subjs)
    sid = unique_subjs(s_idx);
    if any(pairs_aa_comp.subj_id == sid)
        subj_baseline_b_comp.(sprintf('s%d', sid)) = ...
            tr_one_b_b_comp(tr_one_b_b_comp.subj_id == sid, :);
    end
    if any(pairs_aa_iso.subj_id == sid)
        subj_baseline_b_iso.(sprintf('s%d', sid)) = ...
            tr_one_b_b_iso(tr_one_b_b_iso.subj_id == sid, :);
    end
end

%% Compute A-A compared
fprintf('\nComputing A-A compared (%d pairs)...\n', height(pairs_aa_comp));
results_aa_comp = zeros(height(pairs_aa_comp), 7);
n_aa_comp = 0;
for i = 1:height(pairs_aa_comp)
    if mod(i, 50) == 0, fprintf('  %d/%d\n', i, height(pairs_aa_comp)); end
    pair = pairs_aa_comp(i,:);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_comp, skey), continue; end

    fix_1b_match = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_a, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    if sum(map_1b_match(:)) == 0 || sum(map_2b(:)) == 0, continue; end

    % Baseline: other B1 encoding trials, excluding the partner B
    pool = subj_baseline_b_comp.(skey);
    partner_b = regexprep(pair.stim_id, '_A_', '_B_');
    baseline_pool = pool(~strcmp(pool.stim_id, partner_b), :);
    if height(baseline_pool) < 1, continue; end

    match_score = corr(map_1b_match(:), map_2b(:));
    mismatch_scores = nan(height(baseline_pool), 1);
    for j = 1:height(baseline_pool)
        fix_1b_mismatch = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_pool.trial_id(j), :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        if sum(map_1b_mismatch(:)) == 0, continue; end
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    baseline_score = mean(mismatch_scores, 'omitnan');
    n_aa_comp = n_aa_comp + 1;
    results_aa_comp(n_aa_comp,:) = [pair.subj_id, pair.tr_1b_a, pair.tr_2b_a, match_score, baseline_score, match_score - baseline_score, pair.correct];
end
results_aa_comp = results_aa_comp(1:n_aa_comp,:);
fprintf('A-A compared complete: %d valid pairs\n\n', n_aa_comp);

%% Compute A-A isolated
fprintf('Computing A-A isolated (%d pairs)...\n', height(pairs_aa_iso));
results_aa_iso = zeros(height(pairs_aa_iso), 7);
n_aa_iso = 0;
for i = 1:height(pairs_aa_iso)
    if mod(i, 50) == 0, fprintf('  %d/%d\n', i, height(pairs_aa_iso)); end
    pair = pairs_aa_iso(i,:);
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_iso, skey), continue; end

    fix_1b_match = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_a, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    if sum(map_1b_match(:)) == 0 || sum(map_2b(:)) == 0, continue; end

    % Baseline: other B1 encoding trials, excluding the partner B
    pool = subj_baseline_b_iso.(skey);
    partner_b = regexprep(pair.stim_id, '_A_', '_B_');
    baseline_pool = pool(~strcmp(pool.stim_id, partner_b), :);
    if height(baseline_pool) < 1, continue; end

    match_score = corr(map_1b_match(:), map_2b(:));
    mismatch_scores = nan(height(baseline_pool), 1);
    for j = 1:height(baseline_pool)
        fix_1b_mismatch = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_pool.trial_id(j), :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        if sum(map_1b_mismatch(:)) == 0, continue; end
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    baseline_score = mean(mismatch_scores, 'omitnan');
    n_aa_iso = n_aa_iso + 1;
    results_aa_iso(n_aa_iso,:) = [pair.subj_id, pair.tr_1b_a, pair.tr_2b_a, match_score, baseline_score, match_score - baseline_score, pair.correct];
end
results_aa_iso = results_aa_iso(1:n_aa_iso,:);
fprintf('A-A isolated complete: %d valid pairs\n\n', n_aa_iso);

%% Convert to tables and save
results_aa_comp = array2table(results_aa_comp, 'VariableNames', {'subj_id','tr_1b_a','tr_2b_a','match_score','baseline_score','reinst_index','correct'});
results_aa_iso  = array2table(results_aa_iso,  'VariableNames', {'subj_id','tr_1b_a','tr_2b_a','match_score','baseline_score','reinst_index','correct'});

fprintf('Output correct rate: comp=%.1f%%, iso=%.1f%%\n', ...
    100*mean(results_aa_comp.correct, 'omitnan'), ...
    100*mean(results_aa_iso.correct, 'omitnan'));

reinstat_res_aa = struct();
reinstat_res_aa.aa_compared = results_aa_comp;
reinstat_res_aa.aa_isolated = results_aa_iso;
reinstat_res_aa.spatial_params = spatial_params;

fprintf('Saving results...\n');
save(fullfile(res_dir, 'gaze_reinstat_res_aa_b.mat'), 'reinstat_res_aa');
fprintf('Done! A-A: comp=%d, iso=%d\n', n_aa_comp, n_aa_iso);

%% Helper function
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