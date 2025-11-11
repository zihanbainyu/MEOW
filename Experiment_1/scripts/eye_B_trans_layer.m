clear;
clc;
close all;
dbstop if error;

%==========================================================================
% SCRIPT 2: TRANSLATION LAYER
% (Builds the 'in' struct for dataAnalysis.m)
%==========================================================================
fprintf('--- starting script 2: translation layer ---\n');

subj_id = 501;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id); 
results_dir = fullfile(base_dir, 'data', subj_folder);

% --- 1. define constants (FROM YOU) ---
MONITOR_WIDTH_CM = 60;
MONITOR_WIDTH_PX = 1920;
MONITOR_HEIGHT_PX = 1080;
VIEWING_DIST_CM = 100;
SAMPLE_RATE_HZ = 1000;
TRIAL_TIMEOUT_MS = 1500; % 1.5s (from your 2-back script)

% --- 2. load the two data files ---
% a) load the raw samples we just parsed
raw_sample_file = fullfile(results_dir, sprintf('%03d_2_1_raw_samples.mat', subj_id));
if ~exist(raw_sample_file, 'file'), error('run script 1 (parser) first.'); end
fprintf('loading raw samples: %s\n', raw_sample_file);
load(raw_sample_file, 'raw_samples');

% b) load the behavioral analysis file (subj_stats)
analysis_file = fullfile(results_dir, sprintf('sub%03d_concat.mat', subj_id));
fprintf('loading analysis file: %s\n', analysis_file);
load(analysis_file, 'final_data_output');

% --- 3. prepare the 'in' struct ---
fprintf('building "in" struct...\n');
in = [];
run_num = 1; % this is block 1

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
results_2_back = final_data_output.results_2_back;
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
in.trialTypes{run_num} = zeros(num_trials, 1);
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
    in.trialTypes{run_num}(i) = code;
end

% --- 3e. clean up nan trials ---
nan_trials = any(isnan(in.startInds{run_num}), 2);
in.startInds{run_num}(nan_trials,:) = [];
in.trialTypes{run_num}(nan_trials) = [];
fprintf('removed %d nan/junk trials. final trial count: %d\n', sum(nan_trials), length(in.trialTypes{run_num}));

% --- 4. save the final 'in' struct ---
output_file = fullfile(results_dir, sprintf('sub%03d_2_1_for_tool.mat', subj_id));
fprintf('saving final "in" struct to: %s\n', output_file);
save(output_file, 'in');
fprintf('--- script 2 complete. you are ready to run dataAnalysis.m ---\n');

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