
fprintf('attempting to open asc file: %s\n', eye_file);
fid = fopen(eye_file);
if fid == -1, error('file not found.'); end
fprintf('file opened successfully. parsing all events...\n');

% pre-allocate a struct
max_trials = 200;
eye_data(max_trials).trial_id = [];
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

fix_table_vars = {'eye', 'start_time', 'end_time', 'duration', 'x_pos', 'y_pos', 'pupil_size'};
fix_table_types = {'char', 'double', 'double', 'double', 'double', 'double', 'double'};
sac_table_vars = {'eye', 'start_time', 'end_time', 'duration', 'start_x', 'start_y', 'end_x', 'end_y'};
sac_table_types = {'char', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
blink_table_vars = {'eye', 'start_time', 'end_time', 'duration'};
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
eye_data = eye_data(1:current_trial_num); % remove extra pre-allocated rows

fprintf('parsing complete. created "eye_data" struct.\n');

% final check
fprintf('\n--- check data for trial 10 ---\n');
disp(eye_data(10));
fprintf('--- fixations for trial 10 ---\n');
disp(eye_data(10).fixations);