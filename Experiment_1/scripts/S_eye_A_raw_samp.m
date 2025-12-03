clear;
clc;
close all;


subj_id = 501;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id);
results_dir = fullfile(base_dir, 'data', subj_folder);
MONITOR_WIDTH_CM = 60;
MONITOR_WIDTH_PX = 1920;
MONITOR_HEIGHT_PX = 1080;
VIEWING_DIST_CM = 100;
SAMPLE_RATE_HZ = 1000;
TRIAL_TIMEOUT_MS = 1500;
task_name = 'results_2_back';
in = [];

nRuns = 4;
for b = 1:nRuns

    fprintf('block %01d..\n', b);

    asc_file = fullfile(results_dir, sprintf('%03d_2_%01d.asc', subj_id, b));
    output_raw_file = fullfile(results_dir, sprintf('%03d_2_%01d_raw.mat', subj_id, b));
    model_input = fullfile(results_dir, sprintf('%03d_2_pcdm_input.mat', subj_id));

    fid = fopen(asc_file);
    if fid == -1
        warning('Could not open ASC file: %s. Skipping block.', asc_file);
    end
    behav_file = fullfile(results_dir, sprintf('sub%03d_2_back_b%01d.mat', subj_id, b));
    load(behav_file, task_name);

    n_buffer = 500000;
    raw_timestamps = zeros(n_buffer, 1);
    raw_xpos = zeros(n_buffer, 1);
    raw_ypos = zeros(n_buffer, 1);
    raw_pupil = zeros(n_buffer, 1);
    n_event_buffer = 5000;
    event_timestamps = zeros(n_event_buffer, 1);
    event_messages = cell(n_event_buffer, 1);

    sample_count = 0;
    event_count = 0;
    tline = fgetl(fid);


    while ischar(tline)
        if isempty(tline), tline = fgetl(fid); continue; end

        if tline(1) >= '0' && tline(1) <= '9'
            sample_count = sample_count + 1;
            try
                sample_data = sscanf(tline, '%d %f %f %d %*s');
                raw_timestamps(sample_count) = sample_data(1);
                raw_xpos(sample_count) = sample_data(2);
                raw_ypos(sample_count) = sample_data(3);
                raw_pupil(sample_count) = sample_data(4);
            catch
                raw_timestamps(sample_count) = sscanf(tline, '%d %*s');
                raw_xpos(sample_count) = nan;
                raw_ypos(sample_count) = nan;
                raw_pupil(sample_count) = nan;
            end
        elseif startsWith(tline, 'MSG')
            event_count = event_count + 1;
            parts = split(tline);
            event_timestamps(event_count) = str2double(parts{2});
            event_messages{event_count} = tline;
        end
        tline = fgetl(fid);
    end
    fclose(fid);

    raw_samples = [];
    raw_samples.timestamps = raw_timestamps(1:sample_count);
    raw_samples.xpos_px = raw_xpos(1:sample_count);
    raw_samples.ypos_px = raw_ypos(1:sample_count);
    raw_samples.pupil_raw = raw_pupil(1:sample_count);
    raw_samples.events.timestamps = event_timestamps(1:event_count);
    raw_samples.events.messages = event_messages(1:event_count);

    fprintf('converting gaze positions...\n');

    [x_dva, y_dva] = pixels_to_dva(raw_samples.xpos_px, raw_samples.ypos_px, ...
        MONITOR_WIDTH_CM, MONITOR_WIDTH_PX, MONITOR_HEIGHT_PX, VIEWING_DIST_CM);

    in.xPos{b} = x_dva;
    in.yPos{b} = y_dva;
    in.pupilArea{b} = raw_samples.pupil_raw;
    in.sampleRate{b} = SAMPLE_RATE_HZ;
    num_trials = height(results_2_back);
    in.startInds{b} = zeros(num_trials, 2);

    % get event timestamps
    event_msgs = raw_samples.events.messages;
    event_times = raw_samples.events.timestamps;

    fprintf('getting event messages and timestamps...\n');

    for i = 1:num_trials
        trial_msg = sprintf('TRIALID %d', i);
        trial_event_idx = find(contains(event_msgs, trial_msg), 1, 'first');
        if isempty(trial_event_idx)
            warning('could not find TRIALID for trial %d. skipping.', i);
            in.startInds{b}(i, :) = [nan, nan];
            continue;
        end

        % find the 'STIM_ONSET' message for this trial
        stim_event_idx = find(contains(event_msgs(trial_event_idx:end), 'STIM_ONSET'), 1, 'first');
        if isempty(stim_event_idx)
            warning('could not find STIM_ONSET for trial %d. skipping.', i);
            in.startInds{b}(i, :) = [nan, nan];
            continue;
        end
        stim_event_idx = trial_event_idx + stim_event_idx - 1;
        stim_onset_time_elink = event_times(stim_event_idx);

        % get the reaction time
        stim_end_time_elink = stim_onset_time_elink + TRIAL_TIMEOUT_MS;

        % now find the *sample indices* for these times
        start_index = find(raw_samples.timestamps >= stim_onset_time_elink, 1, 'first');
        end_index = find(raw_samples.timestamps >= stim_end_time_elink, 1, 'first');

        if isempty(start_index) || isempty(end_index)
            warning('could not find sample indices for trial %d. skipping.', i);
            in.startInds{b}(i, :) = [nan, nan];
            continue;
        end

        in.startInds{b}(i, :) = [start_index, end_index];
    end

    fprintf('building trial type ...\n');
    in.trialTypes{b} = zeros(1, num_trials);
    for i = 1:num_trials
        cond = results_2_back.condition(i);
        goal = results_2_back.goal(i);

        if cond == "compared"
            if goal == "A-A", code = 1;
            elseif goal == "A-B", code = 2;
            elseif goal == "A-N", code = 3;
            else, code = 0; % junk
            end
        elseif cond == "isolated"
            if goal == "A-A", code = 4;
            elseif goal == "A-B", code = 5;
            elseif goal == "A-N", code = 6;
            else, code = 0; % junk
            end
        elseif cond == "novel"
            if goal == "A-A", code = 7;
            elseif goal == "A-B", code = 8;
            elseif goal == "A-N", code = 9;
            else, code = 0; % junk
            end
        else
            code = 0; % junk
        end
        in.trialTypes{b}(i) = code;
    end
    save(output_raw_file, 'in{b}');
    fprintf('raw file saved\n');

end

save(model_input, 'in');
fprintf('model input saved\n');

function [x_dva, y_dva] = pixels_to_dva(x_px, y_px, screen_width_cm, screen_width_px, screen_height_px, dist_cm)
% 1. calculate pixels per cm
px_per_cm = screen_width_px / screen_width_cm;

% 2. find center of screen (in pixels)
center_x_px = screen_width_px / 2;
center_y_px = screen_height_px / 2;

% First, create empty vectors
x_cm = nan(size(x_px));
y_cm = nan(size(y_px));

% Find *only* the valid, on-screen samples
% (Eyelink often uses < 0 for off-screen, so we check all bounds)
valid_idx = (x_px > 0 & x_px < screen_width_px & y_px > 0 & y_px < screen_height_px);

% 3. convert *only valid* pixels to cm *relative to center*
x_cm(valid_idx) = (x_px(valid_idx) - center_x_px) / px_per_cm;
y_cm(valid_idx) = (y_px(valid_idx) - center_y_px) / px_per_cm;


% 4. use trigonometry (atan) to get the angle
% atan(nan) is nan, which is what we want.
x_rad = atan(x_cm / dist_cm);
y_rad = atan(y_cm / dist_cm);

% 5. convert radians to degrees
x_dva = x_rad * (180 / pi);
y_dva = y_rad * (180 / pi);

fprintf('  helper: conversion complete. (nans are preserved)\n');
end