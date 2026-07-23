%==========================================================================
% Phase 2: Old/new recognition task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================
function [results_table] = E_run_recognition(p, el, sequence_recognition)
%% ========================================================================
% SECTION 1: SET UP
% ========================================================================
% define parameters
is_eyetracking = p.eyetracking == 1;
results_table = sequence_recognition;
num_trials = height(results_table);
% --- Pre-allocate results columns ---
results_table.resp_key = strings(num_trials, 1);
results_table.resp_key(:) = "NA";
results_table.rt = nan(num_trials, 1);
% --- Define key names ---
escape_key = KbName(p.keys.quit);
old_key = KbName(p.keys.same); % 'j' = OLD
new_key = KbName(p.keys.diff); % 'k' = NEW
start_key = KbName('f'); % 'f' to start
%% ------------------------------------------------------------------------
% Apply Text and Keyboard Settings
% -------------------------------------------------------------------------
Screen('TextSize', p.window, p.text_size);
Screen('TextFont', p.window, 'Helvetica');
% Initialize KbQueue for fast, buffered responses
respKeys = zeros(1, 256);
respKeys([old_key, new_key, escape_key]) = 1; % ONLY listen for these keys
KbQueueCreate(p.keys.device, respKeys); 
KbCheck(p.keys.device);
%% ========================================================================
% SECTION 2: BLOCK & TRIAL EXECUTION
% ========================================================================
%------------------------------------------------------------------
% 2A: Start of Block Screen 
%------------------------------------------------------------------
DrawFormattedText(p.window, ...
    ['Part 2: Memory Test\n\n\n' ...
    'You will see a series of images, one at a time.\n\n' ...
    'For each, please decide whether you have seen that exact image today.\n\n\n' ...
    'j = OLD (previously seen)     k = NEW (not seen)\n\n\n' ...
    'Please respond as quickly and accurately as you can.\n\n\n' ...
    'When you are ready, press f to begin.'], ...
    'center', 'center', p.colors.black, [], [], [], 1.2);
Screen('Flip', p.window);
KbReleaseWait(p.keys.device);
% Wait for 'f' key press to start
while true
    [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
    if keyIsDown
        if keyCode(start_key), break;
        elseif keyCode(escape_key), error('USER_ABORT'); end
    end
    WaitSecs(0.001); 
end
% Eyelink: start recording eye movements
if is_eyetracking
    if Eyelink('IsConnected') ~= 1
        error('EYELINK_FATAL: Connection lost.');
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
% Initial fixation before the first trial
draw_fixation_target(p);   % same target as the 1-back / 2-back tasks
Screen('Flip', p.window);
WaitSecs(2);
KbQueueStart(p.keys.device);
for i = 1:num_trials
    trial_info = results_table(i,:);
    
    %------------------------------------------------------------------
    % 2B: Trial loop
    %------------------------------------------------------------------
    if is_eyetracking
        Eyelink('Message', 'TRIALID %d', i);
        Eyelink('command', 'record_status_message "Recognition, Trial %d"', i);
    end
    
    % --------- fixation ------------
    draw_fixation_target(p);   % same target as the 1-back / 2-back tasks
    fix_onset_time = Screen('Flip', p.window);
    if is_eyetracking
        Eyelink('Message', 'FIXATION_ONSET');
    end
    
    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, trial_info.stim_id);
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
    
    KbQueueFlush(p.keys.device);
    % --- present image and collect response
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);
    
    if is_eyetracking
        Eyelink('Message', 'SYNCTIME');
        Eyelink('Message', 'STIM_ONSET %s', char(trial_info.stim_id));
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', char(img_path), p.xCenter, p.yCenter);
    end
    
    key_pressed = "NA";
    response_time = NaN;
    responded = false;

    % self-paced, but has a 4s limit
    deadline = stim_onset_time + 4.0;
    
    % Loop until 4s passes OR a valid response is made
    while GetSecs < deadline && ~responded
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        
        if pressed
            response_key_code = find(firstPress > 0, 1);
            if isempty(response_key_code)
                continue; % A non-monitored key was pressed
            end

            % A monitored key was pressed. Log the time.
            response_time = firstPress(response_key_code) - stim_onset_time;
            
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif response_key_code == old_key
                key_pressed = string(p.keys.same);
                responded = true; % Valid response, stop looping
            elseif response_key_code == new_key
                key_pressed = string(p.keys.diff);
                responded = true; % Valid response, stop looping
            else
                % An invalid (but monitored) key was pressed.
                % Do nothing, let the loop continue.
            end
            
            if is_eyetracking && responded
                rt_ms = response_time * 1000;
                if ~isfinite(rt_ms) || rt_ms < 0, rt_log_value = -999;
                else, rt_log_value = round(rt_ms); end
                Eyelink('Message', 'RESPONSE KEY %s RT %d', char(key_pressed), rt_log_value); 
            end
        end
        WaitSecs(0.001); % Yield
    end
    
    % Clear screen after response
    Screen('Flip', p.window); 
    Screen('Close', img_texture);
    
    % Eyelink: log trial variables
    if is_eyetracking
        Eyelink('Message', '!V TRIAL_VAR stimulus %s', char(trial_info.stim_id));
        Eyelink('Message', '!V TRIAL_VAR condition %s', char(trial_info.condition));
        Eyelink('Message', '!V TRIAL_VAR identity %s', char(trial_info.identity));
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
    
end % end of the trial loop

Screen('Flip', p.window);
WaitSecs(0.05);
if is_eyetracking
    WaitSecs(0.1);
    Eyelink('StopRecording');
end
KbQueueRelease(p.keys.device);
end

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