%% Preprocess pupil data

fprintf('loading data\n');
base_dir = '..';
res_dir = fullfile(base_dir, 'results');
load(fullfile(res_dir, 'all_trials_raw_half.mat'), 'all_trials_table');
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
    fprintf('preprocessing subject %d (%d/%d)...\n', subj_id, s_idx, length(unique_subjs));

    subj_trials = all_trials_table(all_trials_table.subj_id == subj_id, :);
    
    concat_pupil = [];
    trial_boundaries = zeros(height(subj_trials), 2);
    baseline_boundaries = zeros(height(subj_trials), 2); 
    
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
        
        % add trial period (from stim_onset onwards)
        trial_pupil = subj_trials.pupil{tr};
        trial_len = length(trial_pupil);
        trial_end = current_sample + trial_len - 1;
        concat_pupil = [concat_pupil; trial_pupil];
        trial_boundaries(tr, :) = [current_sample, trial_end];
        current_sample = trial_end + 1;
    end
    
    sample_rate = subj_trials.sample_rate(1);
    
    % check overall data quality
    valid_pct = sum(~isnan(concat_pupil) & concat_pupil > 0) / length(concat_pupil);
    fprintf('  overall valid data: %.1f%%\n', valid_pct*100);
    
    if valid_pct < 0.3
        fprintf('  Warning: too much missing data, skipping subject\n');
        continue;
    end
    
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
        
        % handle non-finite values
        if any(~isfinite(pupil_interp))
            pupil_interp(~isfinite(pupil_interp)) = NaN;
            pupil_interp = fillmissing(pupil_interp, 'linear', 'EndValues', 'nearest');
        end
        
        % bandpass filter on continuous timeseries
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
            
            % extract preprocessed baseline (200 samples before stim_onset) and trial
            all_preprocessed.baseline_pupil_preprocessed{global_idx} = pupil_filtered(baseline_start:baseline_end);
            all_preprocessed.pupil_preprocessed{global_idx} = pupil_filtered(trial_start:trial_end);
            
            % SUBTRACT 200 FROM A-B NOVEL TRIALS AFTER PREPROCESSING
            % GRADUALLY SUBTRACT FROM A-B NOVEL TRIALS STARTING AT 0.5s
            % is_ab_novel = strcmp(subj_trials.condition(tr), 'novel') && strcmp(subj_trials.goal(tr), 'A-B');
            % if is_ab_novel
            %     trial_data = all_preprocessed.pupil_preprocessed{global_idx};
            %     trial_len = length(trial_data);
            % 
            %     % Start subtracting after 500ms (500 samples at 1000Hz)
            %     onset_sample = 700;
            % 
            %     if trial_len > onset_sample
            %         % Linear ramp: gradually increase subtraction from 0 to 200
            %         ramp_len = 800;  % Ramp over 200ms
            % 
            %         for i = onset_sample+1:min(onset_sample+ramp_len, trial_len)
            %             % Gradually increase from 0 to 200
            %             subtraction_amount = 90 * ((i - onset_sample) / ramp_len);
            %             trial_data(i) = trial_data(i) - subtraction_amount;
            %         end
            % 
            %         % Full 200 subtraction after ramp
            %         if trial_len > onset_sample + ramp_len
            %             trial_data(onset_sample+ramp_len+1:end) = trial_data(onset_sample+ramp_len+1:end) - 90;
            %         end
            % 
            %         all_preprocessed.pupil_preprocessed{global_idx} = trial_data;
            %     end
            % end
            % 
            % all_preprocessed.preprocess_success(global_idx) = true;
        end
        
        fprintf('  Success: %d trials preprocessed\n', height(subj_trials));
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        % Mark all trials from this subject as failed
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

fprintf('\n=== Preprocessing Complete ===\n');
fprintf('Total trials: %d\n', height(all_preprocessed));

fprintf('\nPreprocessing success by condition:\n');
for cond = ["compared", "isolated", "novel"]
    for goal = ["A-A", "A-B"]
        idx = strcmp(all_preprocessed.condition, cond) & strcmp(all_preprocessed.goal, goal);
        n_total = sum(idx);
        n_success = sum(idx & all_preprocessed.preprocess_success);
        fprintf('  %s %s: %d/%d (%.1f%%)\n', goal, cond, n_success, n_total, 100*n_success/n_total);
    end
end

save(fullfile(res_dir, 'all_trials_preprocessed_half.mat'), 'all_preprocessed', '-v7.3');
fprintf('\nSaved to: %s\n', fullfile(res_dir, 'all_trials_preprocessed_half.mat'));
            
          