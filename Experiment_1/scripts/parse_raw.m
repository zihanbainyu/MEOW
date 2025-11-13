clear;
clc;
close all;
dbstop if error;

%==========================================================================
% SCRIPT 1: RAW ASC SAMPLE PARSER
%==========================================================================
fprintf('--- starting script 1: raw asc parser ---\n');
in = [];

subj_id = 501;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id);
results_dir = fullfile(base_dir, 'data', subj_folder);



for run_num = 1:4

    asc_file = fullfile(results_dir, sprintf('%03d_2_%01d.asc', subj_id, run_num));
    output_file = fullfile(results_dir, sprintf('%03d_2_%01d_raw_v2.mat', subj_id, run_num));

    % --- 1. open the file ---
    fid = fopen(asc_file);
    if fid == -1, error('asc file not found: %s', asc_file); end
    fprintf('file opened. parsing... \n');

    % --- 2. pre-allocate (this is the key to not crashing) ---
    % a 10-minute 1000hz file has 600,000 samples.
    % we'll set a 5-million-sample buffer (overkill is safe)
    n_buffer = 5000000;
    raw_timestamps = zeros(n_buffer, 1);
    raw_xpos = zeros(n_buffer, 1);
    raw_ypos = zeros(n_buffer, 1);
    raw_pupil = zeros(n_buffer, 1);

    % for events
    n_event_buffer = 5000;
    event_timestamps = zeros(n_event_buffer, 1);
    event_messages = cell(n_event_buffer, 1);

    sample_count = 0;
    event_count = 0;
    tline = fgetl(fid);

    % --- 3. read the file line by line ---
    while ischar(tline)

        if isempty(tline)
            tline = fgetl(fid);
            continue;
        end

        % if the first character is a number, it's a sample
        if tline(1) >= '0' && tline(1) <= '9'
            sample_count = sample_count + 1;
            % use sscanf for speed. it's built for this.
            % format: timestamp, x, y, pupil, (skip rest)
            try
                sample_data = sscanf(tline, '%d %f %f %d %*s');
                raw_timestamps(sample_count) = sample_data(1);
                raw_xpos(sample_count) = sample_data(2);
                raw_ypos(sample_count) = sample_data(3);
                raw_pupil(sample_count) = sample_data(4);
            catch
                % a sample line was malformed, or a blink ('.')
                % we set to nan
                sample_count = sample_count + 1;
                raw_timestamps(sample_count) = sscanf(tline, '%d %*s');
                raw_xpos(sample_count) = nan;
                raw_ypos(sample_count) = nan;
                raw_pupil(sample_count) = nan;
            end

            % if it starts with 'MSG', it's an event
        elseif startsWith(tline, 'MSG')
            event_count = event_count + 1;
            parts = split(tline); % split is fine for slow events
            event_timestamps(event_count) = str2double(parts{2});
            event_messages{event_count} = tline;
        end

        tline = fgetl(fid);
    end
    fclose(fid);

    % --- 4. truncate and save ---
    fprintf('parsing complete. read %d samples and %d events.\n', sample_count, event_count);
    raw_samples.timestamps = raw_timestamps(1:sample_count);
    raw_samples.xpos_px = raw_xpos(1:sample_count);
    raw_samples.ypos_px = raw_ypos(1:sample_count);
    raw_samples.pupil_raw = raw_pupil(1:sample_count);

    raw_samples.events.timestamps = event_timestamps(1:event_count);
    raw_samples.events.messages = event_messages(1:event_count);

    fprintf('saving raw .mat file to: %s\n', output_file);
    save(output_file, 'raw_samples');
    fprintf('--- script 1 complete ---\n');

    % --- 1. define constants (FROM YOU) ---
    MONITOR_WIDTH_CM = 60;
    MONITOR_WIDTH_PX = 1920;
    MONITOR_HEIGHT_PX = 1080;
    VIEWING_DIST_CM = 100;
    SAMPLE_RATE_HZ = 1000;
    TRIAL_TIMEOUT_MS = 1500; % 1.5s (from your 2-back script)

    % raw_file  = fullfile(results_dir, sprintf('%03d_2_1_raw_v2.mat', subj_id));
    % load(raw_file, 'raw_samples');
    % b) load the behavioral analysis file (subj_stats)
    analysis_file = fullfile(results_dir, sprintf('sub%03d_2_back_b%01d.mat', subj_id, run_num));
    load(analysis_file, 'results_2_back');
    results_2_back.corr_resp = cellstr(results_2_back.corr_resp);
    results_2_back.resp_key = cellstr(results_2_back.resp_key);
    na_idx = strcmp(results_2_back.resp_key, 'NA');
    results_2_back.resp_key(na_idx) = {'none'};
    results_2_back.correct = strcmp(results_2_back.corr_resp, results_2_back.resp_key);

    % --- 3. prepare the 'in' struct ---
    fprintf('building "in" struct...\n');


    % --- 3a. convert pixels to dva ---
    fprintf('converting pixels to dva...\n');
    [x_dva, y_dva] = pixels_to_dva(raw_samples.xpos_px, raw_samples.ypos_px, ...
        MONITOR_WIDTH_CM, MONITOR_WIDTH_PX, MONITOR_HEIGHT_PX, VIEWING_DIST_CM);

    % --- 3b. populate the simple fields ---
    in.xPos{run_num} = x_dva;
    in.yPos{run_num} = y_dva;
    in.pupilArea{run_num} = raw_samples.pupil_raw;
    in.sampleRate{run_num} = SAMPLE_RATE_HZ;

    % --- 3c. build 'startInds' (the trial boundaries) ---
    fprintf('building "startInds" matrix...\n');
    num_trials = height(results_2_back);
    in.startInds{run_num} = zeros(num_trials, 2);

    % get event timestamps
    event_msgs = raw_samples.events.messages;
    event_times = raw_samples.events.timestamps;

    for i = 1:num_trials
        % find the 'TRIALID' message for this trial
        trial_msg = sprintf('TRIALID %d', i);
        trial_event_idx = find(contains(event_msgs, trial_msg), 1, 'first');
        if isempty(trial_event_idx)
            warning('could not find TRIALID for trial %d. skipping.', i);
            in.startInds{run_num}(i, :) = [nan, nan];
            continue;
        end

        % find the 'STIM_ONSET' message for this trial
        stim_event_idx = find(contains(event_msgs(trial_event_idx:end), 'STIM_ONSET'), 1, 'first');
        if isempty(stim_event_idx)
            warning('could not find STIM_ONSET for trial %d. skipping.', i);
            in.startInds{run_num}(i, :) = [nan, nan];
            continue;
        end
        stim_event_idx = trial_event_idx + stim_event_idx - 1;
        stim_onset_time_elink = event_times(stim_event_idx);

        % get the reaction time
        rt_sec = results_2_back.rt(i);
        if isnan(rt_sec)
            rt_ms = TRIAL_TIMEOUT_MS;
        else
            rt_ms = rt_sec * 1000;
        end
        stim_end_time_elink = stim_onset_time_elink + rt_ms;

        % now find the *sample indices* for these times
        start_index = find(raw_samples.timestamps >= stim_onset_time_elink, 1, 'first');
        end_index = find(raw_samples.timestamps >= stim_end_time_elink, 1, 'first');

        if isempty(start_index) || isempty(end_index)
            warning('could not find sample indices for trial %d. skipping.', i);
            in.startInds{run_num}(i, :) = [nan, nan];
            continue;
        end

        in.startInds{run_num}(i, :) = [start_index, end_index];
    end

    % --- 3d. build 'trialTypes' (the 9-level code) ---
    fprintf('building "trialTypes" vector...\n');
    in.trialTypes{run_num} = zeros(1, num_trials);
    for i = 1:num_trials
        correct = results_2_back.correct(i);
        cond = results_2_back.condition(i);
        goal = results_2_back.goal(i);
        corr_resp = results_2_back.corr_resp(i);
    
        code = 0; % default
    
        if correct == 1
            switch cond
                case "compared"
                    if goal == "A-A" && corr_resp == "j"
                        code = 1;
                    elseif goal == "A-B" && corr_resp == "k"
                        code = 2;
                    end
    
                case "isolated"
                    if goal == "A-A" && corr_resp == "j"
                        code = 3;
                    elseif goal == "A-B" && corr_resp == "k"
                        code = 4;
                    end
    
                case "novel"
                    if goal == "A-A" && corr_resp == "j"
                        code = 5;
                    elseif goal == "A-B" && corr_resp == "k"
                        code = 6;
                    end
            end
        end
    
        in.trialTypes{run_num}(i) = code;
    end

    % --- 3e. clean up nan trials ---
    nan_trials = any(isnan(in.startInds{run_num}), 2);
    in.startInds{run_num}(nan_trials,:) = [];
    in.trialTypes{run_num}(nan_trials) = [];
    fprintf('removed %d nan/junk trials. final trial count: %d\n', sum(nan_trials), length(in.trialTypes{run_num}));

    % --- 4. save the final 'in' struct ---
    output_file = fullfile(results_dir, sprintf('sub%03d_2_%01d_input_v3.mat', subj_id, run_num));
    fprintf('saving final "in" struct to: %s\n', output_file);
    save(output_file, 'in');
    fprintf('--- script 2 complete. you are ready to run dataAnalysis.m ---\n');

end


%
% [d, f] = fitModel(in);
% save(('model_pred_v2.mat'), "d", "f");








function [x_dva, y_dva] = pixels_to_dva(x_px, y_px, screen_width_cm, screen_width_px, screen_height_px, dist_cm)
fprintf('  helper: converting %d samples to dva (bulletproof)...\n', length(x_px));

% 1. calculate pixels per cm
px_per_cm = screen_width_px / screen_width_cm;

% 2. find center of screen (in pixels)
center_x_px = screen_width_px / 2;
center_y_px = screen_height_px / 2;

% --- [THIS IS THE FUCKING FIX] ---
% Your raw data is full of NaNs and 0s (blinks/errors).
% We must *only* convert valid, on-screen pixels.

% First, create empty vectors
x_cm = nan(size(x_px));
y_cm = nan(size(y_px));

% Find *only* the valid, on-screen samples
% (Eyelink often uses < 0 for off-screen, so we check all bounds)
valid_idx = (x_px > 0 & x_px < screen_width_px & y_px > 0 & y_px < screen_height_px);

% 3. convert *only valid* pixels to cm *relative to center*
x_cm(valid_idx) = (x_px(valid_idx) - center_x_px) / px_per_cm;
y_cm(valid_idx) = (y_px(valid_idx) - center_y_px) / px_per_cm;

% All other samples (blinks, errors) remain NaN.
% --- [END FIX] ---

% 4. use trigonometry (atan) to get the angle
% atan(nan) is nan, which is what we want.
x_rad = atan(x_cm / dist_cm);
y_rad = atan(y_cm / dist_cm);

% 5. convert radians to degrees
x_dva = x_rad * (180 / pi);
y_dva = y_rad * (180 / pi);

fprintf('  helper: conversion complete. (nans are preserved)\n');
end