clear;
clc;
close all;
subj_id = 501;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id); 
results_dir = fullfile(base_dir, 'data', subj_folder);
b = 2;

% define tasks to load
tasks = [1, 2];
task_labels = {'1_back', '2_back'};

% combined eye_data struct
all_eye_data = struct([]);
global_trial_idx = 0;

% loop through each task
for t_idx = 1:length(tasks)
    t = tasks(t_idx);
    task_label = task_labels{t_idx};
    
    asc_file = fullfile(results_dir, sprintf('%03d_%01d_%01d.asc', subj_id, t, b));
    
    fid = fopen(asc_file);
    if fid == -1
        warning('File not found: %s', asc_file);
        continue;
    end
    fprintf('Processing task %d (%s)...\n', t, task_label);
    
    % pre-allocate a struct
    max_trials = 200;
    eye_data(max_trials).trial_id = [];
    eye_data(max_trials).task = '';
    eye_data(max_trials).stim_id = 'NA';
    eye_data(max_trials).stim_onset_time = NaN;
    eye_data(max_trials).fixations = [];
    eye_data(max_trials).saccades = [];
    eye_data(max_trials).blinks = [];
    
    % per-trial counters
    current_trial_num = 0;
    fix_count = 0;
    sac_count = 0;
    blink_count = 0;
    
    fix_table_vars = {'eye', 'onset', 'offset', 'duration', 'x', 'y', 'pupil_size'};
    fix_table_types = {'char', 'double', 'double', 'double', 'double', 'double', 'double'};
    sac_table_vars = {'eye', 'onset', 'offset', 'duration', 'start_x', 'start_y', 'end_x', 'end_y'};
    sac_table_types = {'char', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    blink_table_vars = {'eye', 'onset', 'offset', 'duration'};
    blink_table_types = {'char', 'double', 'double', 'double'};
    
    tline = fgetl(fid);
    
    while ischar(tline)
        
        % find the start of a trial
        if startsWith(tline, 'MSG') && contains(tline, 'TRIALID')
            
            % clean up the *previous* trial's tables
            if current_trial_num > 0
                eye_data(current_trial_num).fixations = eye_data(current_trial_num).fixations(1:fix_count, :);
                eye_data(current_trial_num).saccades = eye_data(current_trial_num).saccades(1:sac_count, :);
                eye_data(current_trial_num).blinks = eye_data(current_trial_num).blinks(1:blink_count, :);
            end
            
            % now start the new trial
            current_trial_num = sscanf(tline, 'MSG %*d TRIALID %d');
            if current_trial_num == 0, continue; end % skip trial 0
            
            eye_data(current_trial_num).trial_id = current_trial_num;
            eye_data(current_trial_num).task = task_label;
            
            % reset per-trial counters
            fix_count = 0;
            sac_count = 0;
            blink_count = 0;
            
            % pre-allocate event tables with correct types
            eye_data(current_trial_num).fixations = table('Size', [50, 7], ...
                'VariableTypes', fix_table_types, 'VariableNames', fix_table_vars);
            eye_data(current_trial_num).saccades = table('Size', [50, 8], ...
                'VariableTypes', sac_table_types, 'VariableNames', sac_table_vars);
            eye_data(current_trial_num).blinks = table('Size', [20, 4], ...
                'VariableTypes', blink_table_types, 'VariableNames', blink_table_vars);
        
        % find the stim onset time for this trial
        elseif current_trial_num > 0 && startsWith(tline, 'MSG') && contains(tline, 'STIM_ONSET')
            parts = split(tline);
            eye_data(current_trial_num).stim_onset_time = str2double(parts{2});
            eye_data(current_trial_num).stim_id = parts{4};
            
        % if we are in a trial, grab all events
        elseif current_trial_num > 0
            
            % --- parse fixations ---
            if startsWith(tline, 'EFIX')
                fix_data = sscanf(tline, 'EFIX %c %d %d %d %f %f %d');
                fix_count = fix_count + 1;
                eye_data(current_trial_num).fixations(fix_count,:) = {
                    char(fix_data(1)), fix_data(2), fix_data(3), fix_data(4), ...
                    fix_data(5), fix_data(6), fix_data(7)
                };
                
            % --- parse saccades ---
            elseif startsWith(tline, 'ESACC')
                sac_data = sscanf(tline, 'ESACC %c %d %d %d %f %f %f %f %*f %*d');
                sac_count = sac_count + 1;
                eye_data(current_trial_num).saccades(sac_count,:) = {
                    char(sac_data(1)), sac_data(2), sac_data(3), sac_data(4), ...
                    sac_data(5), sac_data(6), sac_data(7), sac_data(8)
                };
            
            % parse blinks
            elseif startsWith(tline, 'EBLINK')
                blink_data = sscanf(tline, 'EBLINK %c %d %d %d');
                blink_count = blink_count + 1;
                eye_data(current_trial_num).blinks(blink_count,:) = {
                    char(blink_data(1)), blink_data(2), blink_data(3), blink_data(4)
                };
            end
        end
        
        tline = fgetl(fid); 
    end
    fclose(fid);
    
    % clean up the struct
    if current_trial_num > 0
        eye_data(current_trial_num).fixations = eye_data(current_trial_num).fixations(1:fix_count, :);
        eye_data(current_trial_num).saccades = eye_data(current_trial_num).saccades(1:sac_count, :);
        eye_data(current_trial_num).blinks = eye_data(current_trial_num).blinks(1:blink_count, :);
    end
    eye_data = eye_data(1:current_trial_num);
    
    % append to combined data
    all_eye_data = [all_eye_data, eye_data];
    
    fprintf('Task %d (%s) complete: %d trials loaded\n', t, task_label, current_trial_num);
end

% rename for consistency
eye_data = all_eye_data;

fprintf('\nTotal trials loaded: %d\n', length(eye_data));


% create trial-level summary table
trial_summary = table();
for i = 1:length(eye_data)
    trial_summary = [trial_summary; {eye_data(i).trial_id, eye_data(i).task, ...
        eye_data(i).stim_id, eye_data(i).stim_onset_time, ...
        height(eye_data(i).fixations), height(eye_data(i).saccades), height(eye_data(i).blinks)}];
end
trial_summary.Properties.VariableNames = {'trial_id', 'task', 'stim_id', 'stim_onset_time', ...
    'n_fixations', 'n_saccades', 'n_blinks'};

writetable(trial_summary, fullfile(results_dir, sprintf('sub%03d_trial_summary.csv', subj_id)));

% save all fixations
all_fixations = table();
for i = 1:length(eye_data)
    if height(eye_data(i).fixations) > 0
        trial_fix = eye_data(i).fixations;
        trial_fix.trial_id = repmat(eye_data(i).trial_id, height(trial_fix), 1);
        trial_fix.task = repmat({eye_data(i).task}, height(trial_fix), 1);
        trial_fix.stim_id = repmat({eye_data(i).stim_id}, height(trial_fix), 1);
        all_fixations = [all_fixations; trial_fix];
    end
end
writetable(all_fixations, fullfile(results_dir, sprintf('sub%03d_fixations.csv', subj_id)));

% save all saccades
all_saccades = table();
for i = 1:length(eye_data)
    if height(eye_data(i).saccades) > 0
        trial_sac = eye_data(i).saccades;
        trial_sac.trial_id = repmat(eye_data(i).trial_id, height(trial_sac), 1);
        trial_sac.task = repmat({eye_data(i).task}, height(trial_sac), 1);
        trial_sac.stim_id = repmat({eye_data(i).stim_id}, height(trial_sac), 1);
        all_saccades = [all_saccades; trial_sac];
    end
end
writetable(all_saccades, fullfile(results_dir, sprintf('sub%03d_saccades.csv', subj_id)));

% save all blinks
all_blinks = table();
for i = 1:length(eye_data)
    if height(eye_data(i).blinks) > 0
        trial_blink = eye_data(i).blinks;
        trial_blink.trial_id = repmat(eye_data(i).trial_id, height(trial_blink), 1);
        trial_blink.task = repmat({eye_data(i).task}, height(trial_blink), 1);
        trial_blink.stim_id = repmat({eye_data(i).stim_id}, height(trial_blink), 1);
        all_blinks = [all_blinks; trial_blink];
    end
end
writetable(all_blinks, fullfile(results_dir, sprintf('sub%03d_blinks.csv', subj_id)));