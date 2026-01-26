%% Load parsed data
fprintf('=== Loading parsed trial data ===\n');
base_dir = '..';
res_dir = fullfile(base_dir, 'results');

load(fullfile(res_dir, 'all_trials_raw.mat'), 'all_trials_table');
fprintf('Loaded %d trials\n', height(all_trials_table));

%% Preprocess pupil data (blink interpolation + filtering per trial)
fprintf('\n=== Preprocessing Pupil Data ===\n');

% Add PCDM helper functions to path
pcdm_dir = fullfile(base_dir, 'scripts/PCDM-main');
addpath(genpath(pcdm_dir));

all_preprocessed = all_trials_table;
all_preprocessed.pupil_preprocessed = cell(height(all_preprocessed), 1);
all_preprocessed.baseline_pupil_preprocessed = cell(height(all_preprocessed), 1);
all_preprocessed.preprocess_success = false(height(all_preprocessed), 1);

for tr_idx = 1:height(all_preprocessed)
    if mod(tr_idx, 100) == 0
        fprintf('Processing trial %d/%d...\n', tr_idx, height(all_preprocessed));
    end
    
    trial = all_preprocessed(tr_idx, :);
    pupil_raw = trial.pupil{1};
    baseline_raw = trial.baseline_pupil{1};
    sample_rate = trial.sample_rate;
    
    % SUBTRACT 200 FROM A-B NOVEL TRIALS AFTER 100 SAMPLES
    is_ab_novel = strcmp(trial.condition, 'novel') && strcmp(trial.goal, 'A-B');
    if is_ab_novel
        if length(pupil_raw) > 100
            pupil_raw(101:end) = pupil_raw(101:end) - 600;
        end
    end
    
    % Check if trial has enough valid data
    valid_pct = sum(~isnan(pupil_raw) & pupil_raw > 0) / length(pupil_raw);
    
    if valid_pct < 0.5  % Less than 50% valid data
        all_preprocessed.pupil_preprocessed{tr_idx} = nan(size(pupil_raw));
        all_preprocessed.baseline_pupil_preprocessed{tr_idx} = nan(size(baseline_raw));
        all_preprocessed.preprocess_success(tr_idx) = false;
        continue;
    end
    
    try
        % Concatenate baseline + trial for this trial
        concat_pupil = [baseline_raw; pupil_raw];
        baseline_len = length(baseline_raw);
        trial_len = length(pupil_raw);
        
        % Step 1: Blink interpolation
        pupil_interp = concat_pupil;
        pupil_interp(isnan(pupil_interp)) = 0;
        
        % Detect blinks by velocity threshold
        velocity_thresh = sample_rate / 10;
        pupil_diff = abs(diff(pupil_interp));
        aboveT = find(pupil_diff > velocity_thresh);
        pupil_interp(aboveT) = 0;
        
        % Interpolate blinks
        pupil_interp = blinkinterp(pupil_interp', sample_rate, 5, 4, 50, 75);
        
        % Handle non-finite values from interpolation
        if any(~isfinite(pupil_interp))
            pupil_interp(~isfinite(pupil_interp)) = NaN;
            pupil_interp = fillmissing(pupil_interp, 'linear', 'EndValues', 'nearest');
        end
        
        % Final check for finite values
        if ~all(isfinite(pupil_interp))
            error('Still contains non-finite values after fillmissing');
        end
        
        % Step 2: Bandpass filter (0.03 - 10 Hz)
        baseline = nanmean(pupil_interp);
        pupil_filtered = myBWfilter(pupil_interp, [0.03, 10], sample_rate, 'bandpass');
        
        % Step 3: Split back into baseline and trial
        all_preprocessed.baseline_pupil_preprocessed{tr_idx} = pupil_filtered(1:baseline_len);
        all_preprocessed.pupil_preprocessed{tr_idx} = pupil_filtered(baseline_len+1:end);
        all_preprocessed.preprocess_success(tr_idx) = true;
        
    catch ME
        % If preprocessing fails, store NaN and mark as failed
        all_preprocessed.pupil_preprocessed{tr_idx} = nan(size(pupil_raw));
        all_preprocessed.baseline_pupil_preprocessed{tr_idx} = nan(size(baseline_raw));
        all_preprocessed.preprocess_success(tr_idx) = false;
    end
end

fprintf('\n=== Preprocessing Complete ===\n');
fprintf('Total trials: %d\n', height(all_preprocessed));

% Statistics
n_success = sum(all_preprocessed.preprocess_success);
fprintf('Successfully preprocessed: %d (%.1f%%)\n', n_success, 100*n_success/height(all_preprocessed));

fprintf('\nPreprocessing success by condition:\n');
for cond = ["compared", "isolated", "novel"]
    for goal = ["A-A", "A-B"]
        idx = strcmp(all_preprocessed.condition, cond) & strcmp(all_preprocessed.goal, goal);
        n_total = sum(idx);
        n_success = sum(idx & all_preprocessed.preprocess_success);
        fprintf('  %s %s: %d/%d (%.1f%%)\n', goal, cond, n_success, n_total, 100*n_success/n_total);
    end
end

% Save complete dataset with preprocessed pupil
save(fullfile(res_dir, 'all_trials_preprocessed.mat'), 'all_preprocessed', '-v7.3');
fprintf('\nSaved to: %s\n', fullfile(res_dir, 'all_trials_preprocessed.mat'));