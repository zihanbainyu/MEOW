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
    ['Phase 2: Final Task\n\n' ...
    'You will now see a long series of everyday items, one at a time.\n\n' ...
    'Your task is simple:\n\n' ...
    'Have you seen this EXACT item at *any point* in the experiment today?\n\n\n' ...
    'Press  ''' p.keys.same '''  (OLD) if you HAVE seen it.\n\n' ...
    'Press  ''' p.keys.diff '''  (NEW) if you have NOT seen it.\n\n\n' ...
    '\n\n' ...
    'You will have 1 second for each item.\n\n' ...
    'Press f to begin the final task.'], ...
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
    edf_filename = sprintf('%d_r.edf', p.subj_id);
    Eyelink('OpenFile', edf_filename);
    fprintf('EYELINK: opened edf file: %s\n', edf_filename);
    Eyelink('command', 'add_file_preamble_text ''Recognition Task''');
    
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
xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
allCoords = [xCoords; yCoords];
Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
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
    Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
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
    % --- present image and collect response (FIXED DURATION) ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);
    
    if is_eyetracking
        Eyelink('Message', 'SYNCTIME');
        Eyelink('Message', 'STIM_ONSET %s', char(trial_info.stim_id));
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', char(img_path), p.xCenter, p.yCenter);
    end
    
    key_pressed = "NA";
    response_time = NaN;
    responded = false;
    stim_duration = 1.0; 
    
    while GetSecs < stim_onset_time + stim_duration
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        if pressed && ~responded
            responded = true; 
            response_key_code = find(firstPress, 1);
            response_time = firstPress(response_key_code) - stim_onset_time;
            
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif response_key_code == old_key
                key_pressed = string(p.keys.same);
            elseif response_key_code == new_key
                key_pressed = string(p.keys.diff);
            else
                key_pressed = "invalid";
            end
            
            if is_eyetracking
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