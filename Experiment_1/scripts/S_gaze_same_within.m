clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
base_dir = fileparts(script_dir);
res_dir = fullfile(base_dir, 'results');
eye_candidates = {
    fullfile(base_dir, 'data', 'eye_movement_data', 'group_eye_movement_combined.mat')
    fullfile(res_dir, 'group_eye_movement_combined.mat')
};
eye_zip = fullfile(res_dir, 'group_eye_movement_combined.mat.zip');

eye_file = '';
for i = 1:numel(eye_candidates)
    if isfile(eye_candidates{i})
        eye_file = eye_candidates{i};
        break;
    end
end
if isempty(eye_file)
    if isfile(eye_zip)
        unzip_dir = fullfile(tempdir, 'meow_eye_movement_data');
        if ~exist(unzip_dir, 'dir'), mkdir(unzip_dir); end
        unzip(eye_zip, unzip_dir);
        eye_file = fullfile(unzip_dir, 'group_eye_movement_combined.mat');
    else
        error('Could not find group_eye_movement_combined.mat.');
    end
end

load(eye_file, 'Mw');

spatial_params = struct('xres',1920,'yres',1080,'roi_x',[760,1160], ...
    'roi_y',[340,740],'grid_x',20,'grid_y',20,'sigma',2);

task_1b_idx = strcmp(Mw.task, '1_back');
task_2b_idx = strcmp(Mw.task, '2_back');
a_idx = strcmp(Mw.identity, 'A');
same_resp_idx = strcmp(string(Mw.corr_resp), "j");
none_resp_idx = strcmp(string(Mw.corr_resp), "none");
aa_goal_idx = strcmp(Mw.goal, 'A-A');

tr_1b_same = unique(Mw(task_1b_idx & a_idx & same_resp_idx, ...
    {'subj_id','trial_id','stim_id','condition','identity','corr_resp'}));
tr_1b_ref = unique(Mw(task_1b_idx & a_idx & none_resp_idx, ...
    {'subj_id','trial_id','stim_id','condition','identity','corr_resp'}));

pairs_a1 = innerjoin(tr_1b_ref, tr_1b_same, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id','condition'}, ...
    'RightVariables', {'trial_id'});
pairs_a1.Properties.VariableNames = {'subj_id','tr_ref','stim_id','condition','tr_same'};
pairs_a1 = pairs_a1(pairs_a1.tr_same > pairs_a1.tr_ref, :);
pairs_a1.condition = repmat({'repeat'}, height(pairs_a1), 1);

tr_2b_same = unique(Mw(task_2b_idx & aa_goal_idx & a_idx & same_resp_idx, ...
    {'subj_id','trial_id','stim_id','condition','identity','goal','corr_resp'}));
tr_2b_ref = unique(Mw(task_2b_idx & aa_goal_idx & a_idx & none_resp_idx, ...
    {'subj_id','trial_id','stim_id','condition','identity','goal','corr_resp'}));

pairs_a2 = innerjoin(tr_2b_ref, tr_2b_same, 'Keys', {'subj_id','stim_id','condition'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id','condition'}, ...
    'RightVariables', {'trial_id'});
pairs_a2.Properties.VariableNames = {'subj_id','tr_ref','stim_id','condition','tr_same'};
pairs_a2 = pairs_a2(pairs_a2.tr_same > pairs_a2.tr_ref, :);
pairs_a2.condition = cellstr(pairs_a2.condition);

fprintf('A1-A1 same-response pairs: %d\n', height(pairs_a1));
fprintf('A2-A2 same-response pairs: %d\n', height(pairs_a2));

results_a1 = compute_same_item_similarity(pairs_a1, Mw, task_1b_idx, spatial_params);
results_a2 = compute_same_item_similarity(pairs_a2, Mw, task_2b_idx, spatial_params);

same_item_res = struct();
same_item_res.a1a1 = results_a1;
same_item_res.a2a2 = results_a2;
same_item_res.spatial_params = spatial_params;

save(fullfile(res_dir, 'gaze_same_item_within.mat'), 'same_item_res');
fprintf('Saved %s\n', fullfile(res_dir, 'gaze_same_item_within.mat'));


function results = compute_same_item_similarity(pairs, Mw, task_idx, spatial_params)
    rows = zeros(height(pairs), 8);
    row_conditions = cell(height(pairs), 1);
    n = 0;
    subj_ids = unique(pairs.subj_id);

    for s = 1:numel(subj_ids)
        sid = subj_ids(s);
        subj_pairs = pairs(pairs.subj_id == sid, :);
        subj_pool = pairs(pairs.subj_id == sid, :);

        for i = 1:height(subj_pairs)
            pair = subj_pairs(i, :);
            same_condition = strcmp(string(subj_pool.condition), string(pair.condition));
            baseline_pool = subj_pool(same_condition & ...
                subj_pool.tr_ref ~= pair.tr_ref & subj_pool.tr_same ~= pair.tr_same, :);
            if height(baseline_pool) < 1, continue; end

            fix_ref = Mw(task_idx & Mw.subj_id == sid & Mw.trial_id == pair.tr_ref, :);
            fix_same = Mw(task_idx & Mw.subj_id == sid & Mw.trial_id == pair.tr_same, :);
            if height(fix_ref) < 2 || height(fix_same) < 2, continue; end

            map_ref = create_fixation_map(fix_ref.x, fix_ref.y, fix_ref.dur, spatial_params);
            map_same = create_fixation_map(fix_same.x, fix_same.y, fix_same.dur, spatial_params);
            if sum(map_ref(:)) == 0 || sum(map_same(:)) == 0, continue; end

            match_score = corr(map_ref(:), map_same(:));
            mismatch_scores = nan(height(baseline_pool), 1);
            for j = 1:height(baseline_pool)
                fix_mismatch = Mw(task_idx & Mw.subj_id == sid & Mw.trial_id == baseline_pool.tr_ref(j), :);
                if height(fix_mismatch) < 2, continue; end
                map_mismatch = create_fixation_map(fix_mismatch.x, fix_mismatch.y, fix_mismatch.dur, spatial_params);
                if sum(map_mismatch(:)) > 0
                    mismatch_scores(j) = corr(map_mismatch(:), map_same(:));
                end
            end

            baseline_score = mean(mismatch_scores, 'omitnan');
            if isnan(baseline_score), continue; end

            n = n + 1;
            rows(n, :) = [sid, pair.tr_ref, pair.tr_same, match_score, baseline_score, ...
                match_score - baseline_score, 1, height(baseline_pool)];
            row_conditions{n} = condition_to_char(pair.condition);
        end
    end

    rows = rows(1:n, :);
    row_conditions = row_conditions(1:n);
    results = array2table(rows, 'VariableNames', {'subj_id','tr_ref','tr_same', ...
        'match_score','baseline_score','reinst_index','correct','n_null'});
    results.condition = row_conditions;
    results = movevars(results, 'condition', 'After', 'subj_id');
end

function out = condition_to_char(condition_value)
    if iscell(condition_value)
        out = char(condition_value{1});
    else
        out = char(condition_value);
    end
end

function map = create_fixation_map(x, y, dur, params)
    in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & ...
        y >= params.roi_y(1) & y <= params.roi_y(2);
    x = x(in_roi); y = y(in_roi); dur = dur(in_roi);
    if isempty(x), map = zeros(params.grid_y, params.grid_x); return; end

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

    if params.sigma > 0, map = imgaussfilt(map, params.sigma); end
    if sum(map(:)) > 0, map = map / sum(map(:)); end
end