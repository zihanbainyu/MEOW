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
% Per-trial fixation and stimulus onsets, in seconds since the run-start
% trigger (= lead-in onset, so the 2 s lead-in is included).
results_table.fix_onset  = nan(num_trials_in_block, 1);
results_table.stim_onset = nan(num_trials_in_block, 1);
% define key names (Assuming KbName('UnifyKeyNames') was called in main)
escape_key = KbName(p.keys.quit);
same_key = KbName({'1!','1'});      % top-row '1' and numpad '1' (button box)
similar_key = KbName({'2@','2'});   % top-row '2' and numpad '2' (button box)
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
    [sprintf('Block %01d:  One-Back\n\n\n\n', b) ...
    'Compare each object to the one RIGHT BEFORE IT.', ...
    '\n\n\n\n' ...
    'Left button = same.     Right button = similar.     No press = different.' ...
    '\n\n\n\n' ...
    'Pay attention to every object, because it may reappear in later tasks.\n\n\n\n' ...
    'Please let the experimenter know when you are ready to begin.'], 'center', 'center', p.colors.black);
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
KbReleaseWait(p.keys.device);
% Suppress keystrokes from reaching the MATLAB command window for the rest of
% the run. KbQueue still captures responses (it reads the device directly).
% Restored at the end of the run.
ListenChar(2);
fprintf('1-back block %d starting.\n', current_block);
% Lead-in fixation before the first trial of the block, held for
% p.timing.block_lead_in. Its onset is the run's t = 0 (run-start trigger).
draw_fixation_target(p);
draw_photodiode(p, false);
lead_in_onset = Screen('Flip', p.window);
meg_trigger(p, p.trig_codes.run_start_1back);
if is_eyetracking, Eyelink('Message', 'LEAD_IN_ONSET'); end
run_start = lead_in_onset;
WaitSecs('UntilTime', lead_in_onset + p.timing.block_lead_in);

KbQueueStart(p.keys.device);

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
    draw_photodiode(p, false);
    fix_onset_time = Screen('Flip', p.window);
    meg_trigger(p, p.trig_codes.fixation);
    if is_eyetracking
        Eyelink('Message', 'FIXATION_ONSET');
    end
    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, results_table.stim_id(i));
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
    draw_photodiode(p, true);
    % Clear any events that occurred during fixation before presenting stimulus
    KbQueueFlush(p.keys.device);
    % --- present image and collect response ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);
    meg_trigger(p, trial_info.trig);
    % Eyelink: send messages exactly at stimulus onset
    if is_eyetracking
        Eyelink('Message', 'SYNCTIME');
        Eyelink('Message', 'STIM_ONSET %s', char(trial_info.stim_id));
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', char(img_path), p.xCenter, p.yCenter);
    end

    key_pressed = "NA";
    response_time = NaN;
    responded = false;
    % Use KbQueueCheck for efficient and accurate response collection
    while GetSecs < stim_onset_time + p.timing.image_dur
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        if pressed && ~responded
            responded = true;
            % Find the key that was pressed (finds the first non-zero entry)
            response_key_code = find(firstPress, 1);
            response_time = firstPress(response_key_code) - stim_onset_time;
            % Determine which key was pressed
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif any(response_key_code == same_key)
                key_pressed = string(p.keys.same);
                meg_trigger(p, p.trig_codes.resp_same);
            elseif any(response_key_code == similar_key)
                key_pressed = string(p.keys.diff);
                meg_trigger(p, p.trig_codes.resp_similar);
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