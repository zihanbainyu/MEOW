%==========================================================================
% Phase 1: 1-Back Encoding Task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================
function [results_table] = C_run_1_back(p, el, sequence_1_back_block)
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
    [sprintf('1-Back Task'), ...
    '\n\n\n', ...
    'Reminder: Compare the item to the one from ONE TRIAL AGO.', ...
    '\n\n\n', ...
    'Press ' p.keys.same ' = SAME.\n\n' ...
    'Press ' p.keys.diff ' = SIMILAR.\n\n' ...
    'No press = NEW.' ...
    '\n\n\n\n' ...
    'When you are ready, press f to begin'], 'center', 'center', p.colors.black);
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
xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
allCoords = [xCoords; yCoords];
Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
Screen('Flip', p.window);
WaitSecs(2);

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
    Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
    fix_onset_time = Screen('Flip', p.window);
    if is_eyetracking
        Eyelink('Message', 'FIXATION_ONSET');
    end
    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, results_table.stim_id(i));
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
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
        Eyelink('Message', 'TRIAL_RESULT 0');
    end
    %------------------------------------------------------------------
    % 2C: Record trial data
    %------------------------------------------------------------------
    results_table.resp_key(i) = key_pressed;
    results_table.rt(i) = response_time;
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