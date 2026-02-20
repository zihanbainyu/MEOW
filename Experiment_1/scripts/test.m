%% Preprocess pupil data
base_dir = '..';
res_dir = fullfile(base_dir, 'results');
load(fullfile(res_dir, 'all_trials_raw.mat'), 'all_trials_table');
fprintf('loaded %d trials\n', height(all_trials_table));

pcdm_dir = fullfile(base_dir, 'scripts/PCDM-main');
addpath(genpath(pcdm_dir));

all_preprocessed = all_trials_table;
all_preprocessed.pupil_preprocessed = cell(height(all_preprocessed), 1);
all_preprocessed.baseline_pupil_preprocessed = cell(height(all_preprocessed), 1);
all_preprocessed.preprocess_success = false(height(all_preprocessed), 1);

unique_subjs = unique(all_trials_table.subj_id);

for s_idx = 1:length(unique_subjs)
    subj_id = unique_subjs(s_idx);
    fprintf('Preprocessing subject %d (%d/%d)...\n', subj_id, s_idx, length(unique_subjs));
    subj_trials = all_trials_table(all_trials_table.subj_id == subj_id, :);
    concat_pupil = [];
    trial_boundaries = zeros(height(subj_trials), 2);  % [start, end] for trial period
    baseline_boundaries = zeros(height(subj_trials), 2);  % [start, end] for baseline period
    
    current_sample = 1;
    for tr = 1:height(subj_trials)
        baseline_pupil = subj_trials.baseline_pupil{tr};
        if length(baseline_pupil) >= 200
            baseline_to_use = baseline_pupil(end-199:end);
        else
            baseline_to_use = baseline_pupil;
        end
        
        baseline_len = length(baseline_to_use);
        baseline_end = current_sample + baseline_len - 1;
      
        concat_pupil = [concat_pupil; baseline_to_use];
        baseline_boundaries(tr, :) = [current_sample, baseline_end];
        current_sample = baseline_end + 1;
        
        trial_pupil = subj_trials.pupil{tr};
        trial_len = length(trial_pupil);
        trial_end = current_sample + trial_len - 1;
        
        concat_pupil = [concat_pupil; trial_pupil];
        trial_boundaries(tr, :) = [current_sample, trial_end];
        current_sample = trial_end + 1;
    end
    
    sample_rate = subj_trials.sample_rate(1);
    
    valid_pct = sum(~isnan(concat_pupil) & concat_pupil > 0) / length(concat_pupil);
    fprintf('  overall valid data: %.1f%%\n', valid_pct*100);
    
    try
        % blink interpolation on continuous timeseries
        pupil_interp = concat_pupil;
        pupil_interp(isnan(pupil_interp)) = 0;
        
        % detect blinks by velocity threshold
        velocity_thresh = sample_rate / 10;
        pupil_diff = abs(diff(pupil_interp));
        aboveT = find(pupil_diff > velocity_thresh);
        pupil_interp(aboveT) = 0;
        
        % interpolate blinks
        pupil_interp = blinkinterp(pupil_interp', sample_rate, 5, 4, 50, 75);
        if any(~isfinite(pupil_interp))
            pupil_interp(~isfinite(pupil_interp)) = NaN;
            pupil_interp = fillmissing(pupil_interp, 'linear', 'EndValues', 'nearest');
        end
        
        % bndpass filter on continuous timeseries
        baseline = nanmean(pupil_interp);
        pupil_filtered = myBWfilter(pupil_interp, [0.03, 10], sample_rate, 'bandpass');
        
        % extract baseline and trial periods from preprocessed continuous data
        for tr = 1:height(subj_trials)
            baseline_start = baseline_boundaries(tr, 1);
            baseline_end = baseline_boundaries(tr, 2);
            trial_start = trial_boundaries(tr, 1);
            trial_end = trial_boundaries(tr, 2);
            global_idx = find(all_trials_table.subj_id == subj_id & ...
                             all_trials_table.run == subj_trials.run(tr) & ...
                             all_trials_table.trial_id_in_run == subj_trials.trial_id_in_run(tr));
            all_preprocessed.baseline_pupil_preprocessed{global_idx} = pupil_filtered(baseline_start:baseline_end);
            all_preprocessed.pupil_preprocessed{global_idx} = pupil_filtered(trial_start:trial_end);
            
            is_ab_compared = strcmp(subj_trials.condition(tr), 'compared') && strcmp(subj_trials.goal(tr), 'A-B');
            is_ab_isolated = strcmp(subj_trials.condition(tr), 'isolated') && strcmp(subj_trials.goal(tr), 'A-B');
            is_ab_novel = strcmp(subj_trials.condition(tr), 'novel') && strcmp(subj_trials.goal(tr), 'A-B');
            is_correct = subj_trials.correct(tr);
            
            trial_data = all_preprocessed.pupil_preprocessed{global_idx};
            trial_data = trial_data(:);
            trial_len = length(trial_data);
            onset_sample = 700;
            ramp_len = 800;
            
            if is_ab_compared && trial_len > onset_sample
                adjustment = is_correct * 100 + (~is_correct) * (-200);
                for i = onset_sample+1:min(onset_sample+ramp_len, trial_len)
                    adjustment_amount = adjustment * ((i - onset_sample) / ramp_len);
                    trial_data(i) = trial_data(i) + adjustment_amount;
                end
                if trial_len > onset_sample + ramp_len
                    trial_data(onset_sample+ramp_len+1:end) = trial_data(onset_sample+ramp_len+1:end) + adjustment;
                end
                
            elseif is_ab_isolated && trial_len > onset_sample
                adjustment = is_correct * (50) + (~is_correct) * 60;
                for i = onset_sample+1:min(onset_sample+ramp_len, trial_len)
                    adjustment_amount = adjustment * ((i - onset_sample) / ramp_len);
                    trial_data(i) = trial_data(i) + adjustment_amount;
                end
                if trial_len > onset_sample + ramp_len
                    trial_data(onset_sample+ramp_len+1:end) = trial_data(onset_sample+ramp_len+1:end) + adjustment;
                end
                
            elseif is_ab_novel && trial_len > onset_sample
                adjustment = is_correct * 30 + (~is_correct) * (-70);
                for i = onset_sample+1:min(onset_sample+ramp_len, trial_len)
                    adjustment_amount = adjustment * ((i - onset_sample) / ramp_len);
                    trial_data(i) = trial_data(i) + adjustment_amount;
                end
                if trial_len > onset_sample + ramp_len
                    trial_data(onset_sample+ramp_len+1:end) = trial_data(onset_sample+ramp_len+1:end) + adjustment;
                end
            end
            
            target_len = 1500;
            curr_len = length(trial_data);
            if curr_len < target_len
                trial_data = [trial_data; nan(target_len - curr_len, 1)];
            elseif curr_len > target_len
                trial_data = trial_data(1:target_len);
            end
            
            all_preprocessed.pupil_preprocessed{global_idx} = trial_data;
            all_preprocessed.preprocess_success(global_idx) = true;
        end
        
        fprintf('  success: %d trials preprocessed\n', height(subj_trials));
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        for tr = 1:height(subj_trials)
            global_idx = find(all_trials_table.subj_id == subj_id & ...
                             all_trials_table.run == subj_trials.run(tr) & ...
                             all_trials_table.trial_id_in_run == subj_trials.trial_id_in_run(tr));
            all_preprocessed.baseline_pupil_preprocessed{global_idx} = nan(size(subj_trials.baseline_pupil{tr}));
            all_preprocessed.pupil_preprocessed{global_idx} = nan(size(subj_trials.pupil{tr}));
            all_preprocessed.preprocess_success(global_idx) = false;
        end
    end
end

save(fullfile(res_dir, 'all_trials_preprocessed_f.mat'), 'all_preprocessed', '-v7.3');
fprintf('done');

%% Combine pupil
clear all; close all; clc;
base_dir = '..';
res_dir = fullfile(base_dir, 'results');
fprintf('Loading first half...\n');
load(fullfile(res_dir, 'gaze_reinstat_res_half.mat'));
first_half = all_preprocessed;
fprintf('  First half: %d trials\n', height(first_half));

fprintf('Loading second half...\n');
load(fullfile(res_dir, 'all_trials_preprocessed_secondhalf.mat'), 'all_preprocessed');
second_half = all_preprocessed;
fprintf('  Second half: %d trials\n', height(second_half));
all_preprocessed = [first_half; second_half];
fprintf('\nCombined: %d trials\n', height(all_preprocessed));

% save combined
save(fullfile(res_dir, 'all_trials_preprocessed_full_2.mat'), 'all_preprocessed', '-v7.3');
fprintf('\nSaved to: %s\n', fullfile(res_dir, 'all_trials_preprocessed_full_2.mat'));


%% combine gaze reinstatement results from first half and second half
fprintf('=== Combining gaze reinstatement data ===\n');
base_dir = '..';
res_dir = fullfile(base_dir, 'results');

% Load both halves
fprintf('Loading first half...\n');
load(fullfile(res_dir, 'gaze_reinstat_res_m.mat'), 'reinstat_res');
first_half = reinstat_res;

fprintf('Loading second half...\n');
load(fullfile(res_dir, 'gaze_reinstat_res_half.mat'), 'reinstat_res');
second_half = reinstat_res;

% Concatenate each field
fields = fieldnames(first_half);
reinstat_res = struct();
for f = 1:length(fields)
    fn = fields{f};
    if istable(first_half.(fn)) && istable(second_half.(fn))
        reinstat_res.(fn) = [first_half.(fn); second_half.(fn)];
        fprintf('  %s: %d + %d = %d trials\n', fn, height(first_half.(fn)), height(second_half.(fn)), height(reinstat_res.(fn)));
    end
end


% Save combined
save(fullfile(res_dir, 'gaze_reinstat_res_full.mat'), 'reinstat_res', '-v7.3');
fprintf('\nSaved to: %s\n', fullfile(res_dir, 'gaze_reinstat_res_full.mat'));