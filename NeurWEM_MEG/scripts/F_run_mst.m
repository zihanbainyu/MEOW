%==========================================================================
%              Part 2: Post-task MST (old / similar / new)
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%
% MEG version: keeps eyetracking and sends MEG triggers, like the n-back runs.
% Each trial is a jittered fixation (the ITI, 1.0-1.5 s, same as the n-back)
% + a fixed 1.5 s image, during which responses are collected. The image is
% replaced directly by the next trial's fixation (no post-image blank). A
% lead-in fixation is held for p.timing.block_lead_in (2 s), matching the n-back.
%==========================================================================
function [results_table] = F_run_mst(p, el, sequence_mst)
%% ========================================================================
% SECTION 1: SET UP
% ========================================================================
is_eyetracking = p.eyetracking == 1;
results_table = sequence_mst;
num_trials = height(results_table);
results_table.resp_key = strings(num_trials, 1);
results_table.resp_key(:) = "NA";
results_table.rt = nan(num_trials, 1);
% Per-trial fixation and stimulus onsets, in seconds since the run-start
% trigger (= lead-in onset, so the 2 s lead-in is included).
results_table.fix_onset  = nan(num_trials, 1);
results_table.stim_onset = nan(num_trials, 1);

% define key names (Assuming KbName('UnifyKeyNames') was called in main)
old_key     = KbName({'1!','1'});      % OLD: top-row '1' and numpad '1' (button box)
similar_key = KbName({'2@','2'});      % SIMILAR: top-row '2' and numpad '2' (button box)
escape_key  = KbName(p.keys.quit);         % NEW = withhold (no press)
start_key   = KbName('f');

Screen('TextSize', p.window, p.text_size);
Screen('TextFont', p.window, 'Helvetica');

% Initialize KbQueue for fast, buffered responses during the trial loop
respKeys = zeros(1, 256);
respKeys([old_key, similar_key, escape_key]) = 1;
KbQueueCreate(p.keys.device, respKeys);
KbCheck(p.keys.device);

%% ========================================================================
% SECTION 2: RUN EXECUTION
% ========================================================================
%------------------------------------------------------------------
% 2A: Start of Run Screen -- show the MST instruction text (like the
%     one-back / two-back). The experimenter advances it with 'f'.
%------------------------------------------------------------------
instructions(p, 'mst', start_key);   % experimenter advances with 'f'
% Eyelink: start recording eye movements
if is_eyetracking
    if Eyelink('IsConnected') ~= 1
        error('EYELINK_FATAL: Connection lost before MST recording started.');
    end
    Eyelink('Command', 'set_offline_mode');
    WaitSecs(0.05);
    status = Eyelink('StartRecording');
    if status ~= 0
        error('EYELINK_FATAL: StartRecording failed with status %d.', status);
    end
    WaitSecs(0.1);
    Eyelink('Message', 'TRIAL_RESULT 0');
end
KbReleaseWait(p.keys.device);
% Suppress keystrokes from reaching the MATLAB command window for the rest of
% the run. KbQueue still captures responses (it reads the device directly).
% Restored at the end of the run.
ListenChar(2);
fprintf('MST starting.\n');
% Lead-in fixation before the first trial, held for p.timing.block_lead_in.
% Its onset is the run's t = 0 (run-start trigger).
draw_fixation_target(p);
draw_photodiode(p, false);
lead_in_onset = Screen('Flip', p.window);
meg_trigger(p, p.trig_codes.run_start_mst);
if is_eyetracking, Eyelink('Message', 'LEAD_IN_ONSET'); end
run_start = lead_in_onset;
WaitSecs('UntilTime', lead_in_onset + p.timing.block_lead_in);

KbQueueStart(p.keys.device);

for i = 1:num_trials
    trial_info = results_table(i,:);
    %------------------------------------------------------------------
    % 2B: Trial loop
    %------------------------------------------------------------------
    if is_eyetracking
        Eyelink('Message', 'TRIALID %d', i);
        Eyelink('command', 'record_status_message "MST, Trial %d"', i);
    end
    % --------- fixation ------------
    draw_fixation_target(p);
    draw_photodiode(p, false);
    fix_onset_time = Screen('Flip', p.window);
    meg_trigger(p, p.trig_codes.fixation);
    if is_eyetracking, Eyelink('Message', 'FIXATION_ONSET'); end
    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, results_table.stim_id(i));
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
    draw_photodiode(p, true);
    % Clear any events that occurred during fixation before presenting stimulus
    KbQueueFlush(p.keys.device);
    % --- present image for exactly mst_image_dur (responses collected during it) ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);
    meg_trigger(p, trial_info.trig);
    if is_eyetracking
        Eyelink('Message', 'SYNCTIME');
        Eyelink('Message', 'STIM_ONSET %s', char(trial_info.stim_id));
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', char(img_path), p.xCenter, p.yCenter);
    end

    key_pressed = "NA";
    response_time = NaN;
    responded = false;
    % Response window = the image presentation (mst_image_dur). The image stays
    % on screen throughout and is replaced by the next trial's fixation (or, for
    % the last trial, the tail fixation). No post-image blank.
    while GetSecs < stim_onset_time + p.timing.mst_image_dur
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        if pressed && ~responded
            responded = true;
            response_key_code = find(firstPress, 1);
            response_time = firstPress(response_key_code) - stim_onset_time;
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif any(response_key_code == old_key)
                key_pressed = string(p.keys.mst_old);
                meg_trigger(p, p.trig_codes.resp_same);
            elseif any(response_key_code == similar_key)
                key_pressed = string(p.keys.mst_similar);
                meg_trigger(p, p.trig_codes.resp_similar);
            else
                key_pressed = "invalid";
            end
            if is_eyetracking
                rt_ms = response_time * 1000;
                if ~isfinite(rt_ms) || rt_ms < 0
                    rt_log_value = -999;
                else
                    rt_log_value = round(rt_ms);
                end
                Eyelink('Message', 'RESPONSE KEY %s RT %d', char(key_pressed), rt_log_value);
            end
        end
    end
    Screen('Close', img_texture);
    % Eyelink: log trial variables
    if is_eyetracking
        Eyelink('Message', '!V TRIAL_VAR stimulus %s', char(trial_info.stim_id));
        Eyelink('Message', '!V TRIAL_VAR trial_type %s', char(trial_info.trial_type));
        Eyelink('Message', '!V TRIAL_VAR corr_resp %s', char(trial_info.corr_resp));
        Eyelink('Message', '!V TRIAL_VAR response %s', char(key_pressed));
        Eyelink('Message', 'TRIAL_RESULT 0');
    end
    %------------------------------------------------------------------
    % 2C: Record trial data
    %------------------------------------------------------------------
    results_table.resp_key(i) = key_pressed;
    results_table.rt(i) = response_time;
    results_table.fix_onset(i)  = fix_onset_time  - run_start;
    results_table.stim_onset(i) = stim_onset_time - run_start;
end % end of the trial loop

% --- Clear screen before stopping recording ---
Screen('FillRect', p.window, p.colors.bgcolor);
draw_photodiode(p, false);
Screen('Flip', p.window);
meg_trigger(p, p.trig_codes.run_end);
WaitSecs(0.05);

if is_eyetracking
    WaitSecs(0.1);
    Eyelink('StopRecording');
end
% Release the KbQueue resources after the trial loop
KbQueueRelease(p.keys.device);
ListenChar(0);   % restore keystrokes to MATLAB between runs

end

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function draw_fixation_target(p)
% Thaler, Schutz, Goodale & Gegenfurtner (2013, Vision Research):
% combined bullseye-and-crosshair target ("ABC"), the most stable for
% steady fixation. Outer disc (d1) split into quadrants by a background-
% coloured crosshair (width d2), with a central dot (d2) on top.
d1 = p.fix_dot_d1;   % outer disc diameter (px)
d2 = p.fix_dot_d2;   % central dot diameter / crosshair width (px)
cx = p.xCenter; cy = p.yCenter;
% outer disc
Screen('FillOval', p.window, p.fix_dot_color, [cx-d1/2, cy-d1/2, cx+d1/2, cy+d1/2]);
% crosshair in background colour, cutting the disc into four quadrants
Screen('FillRect', p.window, p.colors.bgcolor, [cx-d1/2, cy-d2/2, cx+d1/2, cy+d2/2]); % horizontal bar
Screen('FillRect', p.window, p.colors.bgcolor, [cx-d2/2, cy-d1/2, cx+d2/2, cy+d1/2]); % vertical bar
% central dot
Screen('FillOval', p.window, p.fix_dot_color, [cx-d2/2, cy-d2/2, cx+d2/2, cy+d2/2]);
end
