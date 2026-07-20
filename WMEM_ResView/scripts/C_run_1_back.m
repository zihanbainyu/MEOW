%==========================================================================
% Phase 1: 1-Back Encoding Task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================
function [results_table] = C_run_1_back(p, el, sequence_1_back_block, b)
%% ========================================================================
% SECTION 1: SET UP
% ========================================================================
% define parameters
is_eyetracking = p.eyetracking == 1;
results_table = sequence_1_back_block;
num_trials_in_block = height(results_table);
results_table.resp_key = strings(num_trials_in_block, 1);
results_table.resp_key(:) = "NA";
results_table.rt = nan(num_trials_in_block, 1);
results_table.fix_broken = false(num_trials_in_block, 1);  % gaze left window during viewing
results_table.gate_timeout = false(num_trials_in_block, 1); % onset gate timed out
% define key names (Assuming KbName('UnifyKeyNames') was called in main)
escape_key = KbName(p.keys.quit);
same_key = KbName(p.keys.same);
similar_key = KbName(p.keys.diff);
start_key = KbName('f');
current_block = results_table.block(1);
total_blocks = p.nBlocks;
%% ------------------------------------------------------------------------
% Apply Text and Keyboard Settings
% -------------------------------------------------------------------------
Screen('TextSize', p.window, p.text_size);
Screen('TextFont', p.window, 'Helvetica');
% Initialize KbQueue for fast, buffered responses during the trial loop
respKeys = zeros(1, 256);
respKeys([same_key, similar_key, escape_key]) = 1;
KbQueueCreate(p.keys.device, respKeys); % Create queue for the target device
% Clear any initial key presses before the experiment begins
KbCheck(p.keys.device);

%% ========================================================================
% SECTION 2: BLOCK & TRIAL EXECUTION
% ========================================================================
%------------------------------------------------------------------
% 2A: Start of Block Screen (Fixed DrawFormattedText syntax)
%------------------------------------------------------------------
DrawFormattedText(p.window, ...
    [sprintf('Block %01d of 4\n\n\n\n', b) ...
    '1-Back Task', ...
    '\n\n\n', ...
    'Reminder: Compare each image against the one from ONE TRIAL AGO.', ...
    '\n\n\n', ...
    'Please remain fixated at the dot at the center of the screen while viewing each image.', ...
    '\n\n\n', ...
    'Press ' p.keys.same ' = SAME.\n\n' ...
    'Press ' p.keys.diff ' = SIMILAR.\n\n' ...
    'No press = NEW.' ...
    '\n\n\n\n' ...
    'Please place your index and middle fingers on j and k, and press f to begin'], 'center', 'center', p.colors.black);
Screen('Flip', p.window);
KbReleaseWait(p.keys.device);
% Wait for 'g' key press on the specified device to start
while true
    [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
    if keyIsDown
        if keyCode(start_key), break;
        elseif keyCode(escape_key), error('USER_ABORT'); end
    end
    WaitSecs(0.001); % Yield CPU for a millisecond
end
% Eyelink: start recording eye movements
if is_eyetracking
    if Eyelink('IsConnected') ~= 1
        error('EYELINK_FATAL: Connection lost before block recording started.');
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
% Initial fixation before the first trial of the block
draw_fixation_target(p);
Screen('Flip', p.window);
WaitSecs(2);

KbQueueStart(p.keys.device);

% Determine which eye is tracked, for live gaze reads (sample index = eye+1)
eye_idx = 2; % default to right eye
if is_eyetracking
    tracked = Eyelink('EyeAvailable');
    if tracked == el.BINOCULAR, tracked = el.RIGHT_EYE; end
    if tracked == el.LEFT_EYE, eye_idx = 1; else, eye_idx = 2; end
end

for i = 1:height(results_table)
    trial_info = results_table(i,:);
    %------------------------------------------------------------------
    % 2B: Trial loop for the current block
    %------------------------------------------------------------------

    % Eyelink: mark trial start and perform drift check
    if is_eyetracking
        Eyelink('Message', 'TRIALID %d', i);
        Eyelink('command', 'record_status_message "Block %d, Trial %d"', current_block, i);
        % EyelinkDoDriftCorrection(el);
    end
    % --------- fixation ------------
    draw_fixation_target(p);
    fix_onset_time = Screen('Flip', p.window);
    if is_eyetracking
        Eyelink('Message', 'FIXATION_ONSET');
        % --- onset gate: wait for stable central fixation before stimulus ---
        gate_start = GetSecs;
        in_since = NaN;
        while true
            [gx, gy, valid] = get_gaze(eye_idx);
            now_t = GetSecs;
            if valid && hypot(gx - p.xCenter, gy - p.yCenter) <= p.fix_gate_tol_px
                if isnan(in_since), in_since = now_t; end
                if (now_t - in_since) >= p.fix_gate_ms / 1000, break; end
            else
                in_since = NaN; % excursion or blink resets the hold
            end
            if (now_t - gate_start) >= p.fix_gate_timeout
                results_table.gate_timeout(i) = true;
                Eyelink('Message', 'FIX_GATE_TIMEOUT');
                break;
            end
            WaitSecs(0.001);
        end
        Eyelink('Message', 'FIX_GATE_PASSED');
        fix_onset_time = GetSecs; % re-anchor so fix_duration runs from stable fixation
    end
    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, results_table.stim_id(i));
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
    draw_fixation_target(p);   % restricted-viewing target overlaid on the image
    % Clear any events that occurred during fixation before presenting stimulus
    KbQueueFlush(p.keys.device);
    % --- present image and collect response ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);
    % Eyelink: send messages exactly at stimulus onset
    if is_eyetracking
        Eyelink('Message', 'SYNCTIME');
        Eyelink('Message', 'STIM_ONSET %s', char(trial_info.stim_id));
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', char(img_path), p.xCenter, p.yCenter);
    end

    key_pressed = "NA";
    response_time = NaN;
    responded = false;
    fix_broken = false;
    out_since = NaN;
    % Use KbQueueCheck for efficient and accurate response collection
    while GetSecs < stim_onset_time + p.timing.image_dur
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        % --- gaze monitoring: soft break detection (blink-aware) ---
        if is_eyetracking && ~fix_broken
            [gx, gy, valid] = get_gaze(eye_idx);
            if valid
                if hypot(gx - p.xCenter, gy - p.yCenter) > p.fix_tol_px
                    if isnan(out_since), out_since = GetSecs; end
                    if (GetSecs - out_since) >= p.fix_break_ms / 1000
                        fix_broken = true;
                        Eyelink('Message', 'FIX_BREAK');
                    end
                else
                    out_since = NaN; % back inside the window
                end
            else
                out_since = NaN; % blink / missing sample never counts as a break
            end
        end
        if pressed && ~responded
            responded = true;
            % Find the key that was pressed (finds the first non-zero entry)
            response_key_code = find(firstPress, 1);
            response_time = firstPress(response_key_code) - stim_onset_time;
            % Determine which key was pressed
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif response_key_code == same_key
                key_pressed = string(p.keys.same);
            elseif response_key_code == similar_key
                key_pressed = string(p.keys.diff);
            else
                key_pressed = "invalid";
            end
            
            if is_eyetracking
                % Calculate RT in milliseconds
                rt_ms = response_time * 1000;
                
                % Check if RT is a valid finite number. If not, set a safe placeholder (e.g., -1).
                if ~isfinite(rt_ms) || rt_ms < 0 
                    rt_log_value = -999; % Log a clearly invalid number
                else
                    rt_log_value = round(rt_ms); % Use the rounded value
                end
            
                % Log the message using the safe integer value
                Eyelink('Message', 'RESPONSE KEY %s RT %d', char(key_pressed), rt_log_value); 
            end
        end
    end
    Screen('Close', img_texture);
    % Eyelink: log trial variables
    if is_eyetracking
        % send variables to Data Viewer for analysis
        Eyelink('Message', '!V TRIAL_VAR stimulus %s', char(trial_info.stim_id));
        Eyelink('Message', '!V TRIAL_VAR condition %s', char(trial_info.condition));
        Eyelink('Message', '!V TRIAL_VAR identity %s', char(trial_info.identity));
        Eyelink('Message', '!V TRIAL_VAR corr_resp %s', char(trial_info.corr_resp));
        Eyelink('Message', '!V TRIAL_VAR response %s', char(key_pressed));
        Eyelink('Message', '!V TRIAL_VAR block %d', trial_info.block);
        Eyelink('Message', '!V TRIAL_VAR fix_broken %d', fix_broken);
        Eyelink('Message', 'TRIAL_RESULT 0');
    end
    %------------------------------------------------------------------
    % 2C: Record trial data
    %------------------------------------------------------------------
    results_table.resp_key(i) = key_pressed;
    results_table.rt(i) = response_time;
    results_table.fix_broken(i) = fix_broken;
end % end of the trial loop

% --- Clear screen before stopping recording ---
Screen('FillRect', p.window, p.colors.bgcolor);
Screen('Flip', p.window);
WaitSecs(0.05);

if is_eyetracking
    WaitSecs(0.1);
    Eyelink('StopRecording');
end
% Release the KbQueue resources after the trial loop
KbQueueRelease(p.keys.device);

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

function [gx, gy, valid] = get_gaze(eye_idx)
% Read the newest linked gaze sample. Returns valid = false when no fresh
% sample is available or the eye is lost (blink / missing = -32768).
gx = NaN; gy = NaN; valid = false;
if Eyelink('NewFloatSampleAvailable') > 0
    s = Eyelink('NewestFloatSample');
    gx = s.gx(eye_idx); gy = s.gy(eye_idx);
    if ~isempty(gx) && gx ~= -32768 && gy ~= -32768
        valid = true;
    end
end
end