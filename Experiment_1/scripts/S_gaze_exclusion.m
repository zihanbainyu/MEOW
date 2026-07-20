% S_gaze_exclusion.m
% Apply 70% trial-validity criterion for gaze-based analyses.
% A trial is valid if it has >= 2 fixations within the stimulus ROI.
% A participant is retained if >= 70% of their trials are valid.

clear; clc;

base_dir = '..';
res_dir  = fullfile(base_dir, 'results');

roi_x = [760, 1160];
roi_y = [340, 740];
min_fix_per_trial = 2;
subj_threshold    = 0.70;

load(fullfile(res_dir, 'group_eye_movement_combined.mat'));
fn = fieldnames(load(fullfile(res_dir, 'group_eye_movement_combined.mat')));
T = eval(fn{1});

in_roi = T.x >= roi_x(1) & T.x <= roi_x(2) & ...
         T.y >= roi_y(1) & T.y <= roi_y(2);
T_roi = T(in_roi, :);

[trial_keys, ~, g] = unique(T_roi(:, {'subj_id','task','block','trial_id'}), 'rows');
fix_per_trial = accumarray(g, 1);
trial_keys.n_fix_in_roi = fix_per_trial;
trial_keys.valid = trial_keys.n_fix_in_roi >= min_fix_per_trial;

[all_trials, ~, g_all] = unique(T(:, {'subj_id','task','block','trial_id'}), 'rows');
valid_map = containers.Map( ...
    string(trial_keys.subj_id) + "|" + string(trial_keys.task) + "|" + ...
        string(trial_keys.block) + "|" + string(trial_keys.trial_id), ...
    num2cell(trial_keys.valid));
all_trials.valid = false(height(all_trials), 1);
for i = 1:height(all_trials)
    k = string(all_trials.subj_id(i)) + "|" + string(all_trials.task{i}) + "|" + ...
        string(all_trials.block(i)) + "|" + string(all_trials.trial_id(i));
    if isKey(valid_map, k)
        all_trials.valid(i) = valid_map(k);
    end
end

subs = unique(all_trials.subj_id);
n_sub = numel(subs);
summary = table('Size', [n_sub, 4], ...
    'VariableTypes', {'double','double','double','logical'}, ...
    'VariableNames', {'subj_id','n_trials','prop_valid','retained'});

for i = 1:n_sub
    sid = subs(i);
    mask = all_trials.subj_id == sid;
    summary.subj_id(i)    = sid;
    summary.n_trials(i)   = sum(mask);
    summary.prop_valid(i) = mean(all_trials.valid(mask));
    summary.retained(i)   = summary.prop_valid(i) >= subj_threshold;
end

summary = sortrows(summary, 'prop_valid');
disp(summary);
fprintf('\nRetained: %d / %d participants (threshold = %.0f%% valid trials)\n', ...
    sum(summary.retained), n_sub, subj_threshold*100);
fprintf('Excluded participants: %s\n', mat2str(summary.subj_id(~summary.retained)'));

writetable(summary, fullfile(res_dir, 'gaze_exclusion_summary.csv'));
save(fullfile(res_dir, 'gaze_exclusion_summary.mat'), 'summary', ...
    'roi_x', 'roi_y', 'min_fix_per_trial', 'subj_threshold');
