%==========================================================================
%                  Phase 2: 2-Back Retrieval Task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================
function [results_table] = D_run_2_back(p, el, sequence_2_back_block, b)
%% ========================================================================
%  SECTION 1: SET UP
%  ========================================================================

is_eyetracking = p.eyetracking == 1;
results_table = sequence_2_back_block;
num_trials_in_block = height(results_table);
results_table.resp_key = strings(num_trials_in_block, 1);
results_table.resp_key(:) = "NA";
results_table.rt = nan(num_trials_in_block, 1);

% Define key codes
escape_key = KbName(p.keys.quit);
same_key = KbName({'1!','1'});      % top-row '1' and numpad '1' (button box)
similar_key = KbName({'2@','2'});   % top-row '2' and numpad '2' (button box)
start_key = KbName('f');
% Scanner trigger key ("5" marks scan onset). Configurable via p.keys.trigger;
% defaults to the top-row 5 that most trigger boxes emulate.
if isfield(p.keys, 'trigger'), trig_name = p.keys.trigger; else, trig_name = '5%'; end
trigger_key = KbName(trig_name);

current_block = results_table.block(1);
total_blocks = p.nBlocks;

%% ------------------------------------------------------------------------
%  Apply Text Settings & Initialize KbQueue
% -------------------------------------------------------------------------
Screen('TextSize', p.window, p.text_size); % Set font size
Screen('TextFont', p.window, 'Helvetica'); % Set font style

% Initialize KbQueue for fast, buffered responses during the trial loop
respKeys = zeros(1, 256);
respKeys([same_key, similar_key, escape_key]) = 1;
KbQueueCreate(p.keys.device, respKeys); % Create queue for the target device

% Clear any initial key presses before the experiment begins
KbCheck(p.keys.device);
%% ========================================================================
%  SECTION 2: BLOCK & TRIAL EXECUTION
%  ========================================================================
%------------------------------------------------------------------
% 2A: Start of Block Screen
%------------------------------------------------------------------
DrawFormattedText(p.window, ...
    [sprintf('Block %01d of %01d:  Two-Back\n\n\n\n', b, total_blocks) ...
    'Please compare each object to the one TWO IMAGES BACK.', ...
    '\n\n\n\n' ...
    'Index finger = same.     Middle finger = similar.     No press = different.' ...
    '\n\n\n\n' ...
    'Please let the experimenter know when you are ready to begin.'], 'center', 'center', p.colors.black);
Screen('Flip', p.window);
KbReleaseWait(p.keys.device);

% Wait for 'f' key press on the specified device to start
while true
    [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
    if keyIsDown
        if keyCode(start_key), break;
        elseif keyCode(escape_key), error('USER_ABORT'); end
    end
    WaitSecs(0.001); % Yield CPU for a millisecond
end
KbReleaseWait(p.keys.device);

% Eyelink: start recording eye movements for the block
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

% --- Wait for the scanner trigger ("5") before starting the run ---
wait_for_scanner(p, trigger_key, escape_key);
% Suppress keystrokes from reaching the MATLAB command window for the rest of
% the run, so the scanner's per-TR 5s stop flooding it. KbQueue still captures
% responses (it reads the device directly). Restored at the end of the run.
ListenChar(2);
fprintf('Scanner trigger received; 2-back block %d starting.\n', current_block);

% Lead-in fixation before the first trial of the block. The fixation target is
% held for p.timing.block_lead_in, giving the BOLD signal time to settle before
% the first real trial.
draw_fixation_target(p);
lead_in_onset = Screen('Flip', p.window);
if is_eyetracking, Eyelink('Message', 'LEAD_IN_ONSET'); end
WaitSecs('UntilTime', lead_in_onset + p.timing.block_lead_in);

% Start the KbQueue to collect responses during the trial loop
KbQueueStart(p.keys.device);

%------------------------------------------------------------------
% 2B: Trial loop for the current block
%------------------------------------------------------------------

for i = 1:height(results_table)

    trial_info = results_table(i,:);

    % Eyelink: mark trial start and perform drift check
    if is_eyetracking
        Eyelink('Message', 'TRIALID %d', i);
        Eyelink('command', 'record_status_message "Block %d, Trial %d"', current_block, i);
        % EyelinkDoDriftCorrection(el);
    end

    % ---------fixation------------
    draw_fixation_target(p);
    fix_onset_time = Screen('Flip', p.window);
    if is_eyetracking, Eyelink('Message', 'FIXATION_ONSET'); end

    % ---------stimulus------------
    img_path = fullfile(p.stim_dir, results_table.stim_id(i));
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);

    % Clear any events that occurred during fixation before presenting stimulus
    KbQueueFlush(p.keys.device);

    % --- present image and collect response ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);

    %   Eyelink: send messages exactly at stimulus onset (Applied char() fix)
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

            % Find the key that was pressed first
            response_key_code = find(firstPress, 1);
            response_time = firstPress(response_key_code) - stim_onset_time;

            % Determine which key was pressed
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif any(response_key_code == same_key)
                key_pressed = string(p.keys.same);
            elseif any(response_key_code == similar_key)
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

    %   Eyelink: log trial variables and mark trial end
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
end % end of trial loop

% --- Tail fixation: hold the target so the final trial's response is fully
% captured before the run ends (mirrors the lead-in). ---
if isfield(p.timing, 'block_tail'), tail_dur = p.timing.block_tail; else, tail_dur = 6; end
draw_fixation_target(p);
tail_onset = Screen('Flip', p.window);
if is_eyetracking, Eyelink('Message', 'TAIL_ONSET'); end
WaitSecs('UntilTime', tail_onset + tail_dur);

% --- Clear screen before stopping recording ---
Screen('FillRect', p.window, p.colors.bgcolor);
Screen('Flip', p.window);
WaitSecs(0.05);

% Eyelink BLOCK RECORDING END
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
function wait_for_scanner(p, trigger_key, escape_key)
% Hold after the experimenter starts the run until the scanner sends its first
% trigger ("5"), then return so the run begins. Escape still aborts. Runs
% before KbQueueStart, so KbCheck is safe here.
Screen('TextSize', p.window, p.text_size);
Screen('TextFont', p.window, 'Helvetica');
DrawFormattedText(p.window, 'Waiting for the scanner to start...', 'center', 'center', p.colors.black);
Screen('Flip', p.window);
KbReleaseWait(p.keys.device); % clear the experimenter's start key
while true
    [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
    if keyIsDown
        if keyCode(trigger_key)
            break;
        elseif keyCode(escape_key)
            error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
        end
    end
    WaitSecs(0.001); % yield CPU
end
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