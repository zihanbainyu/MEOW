%==========================================================================
%              Part 2: Post-task MST (old / similar / new)
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%
% Scanner version: keeps eyetracking and waits for the scanner trigger, like
% the n-back runs. Each trial is a jittered fixation (0.75 +/- 0.25 s, same as
% the n-back) + a fixed 1.5 s image + a 0.5 s blank during which responses are
% still accepted (2 s response window in total). Lead-in and tail fixations are
% held for p.timing.block_lead_in / block_tail (6 s), matching the n-back.
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

% define key names (Assuming KbName('UnifyKeyNames') was called in main)
old_key     = KbName({'1!','1'});      % OLD: top-row '1' and numpad '1' (button box)
similar_key = KbName({'2@','2'});      % SIMILAR: top-row '2' and numpad '2' (button box)
escape_key  = KbName(p.keys.quit);         % NEW = withhold (no press)
start_key   = KbName('f');
% Scanner trigger key ("5" marks scan onset). Configurable via p.keys.trigger;
% defaults to the top-row 5 that most trigger boxes emulate.
if isfield(p.keys, 'trigger'), trig_name = p.keys.trigger; else, trig_name = '5%'; end
trigger_key = KbName(trig_name);

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
% 2A: Start of Run Screen -- show the MST instruction figure (like the
%     one-back / two-back). Displays ins_mst and advances on 'f'.
%------------------------------------------------------------------
run_instructions(p, {'ins_mst'});
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
% --- Wait for the scanner trigger ("5") before starting the run ---
wait_for_scanner(p, trigger_key, escape_key);
% Suppress keystrokes from reaching the MATLAB command window for the rest of
% the run, so the scanner's per-TR 5s stop flooding it. KbQueue still captures
% responses (it reads the device directly). Restored at the end of the run.
ListenChar(2);
fprintf('Scanner trigger received; MST starting.\n');
% Lead-in fixation before the first trial. The fixation target is held for
% p.timing.block_lead_in, giving the BOLD signal time to settle.
draw_fixation_target(p);
lead_in_onset = Screen('Flip', p.window);
if is_eyetracking, Eyelink('Message', 'LEAD_IN_ONSET'); end
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
    fix_onset_time = Screen('Flip', p.window);
    if is_eyetracking, Eyelink('Message', 'FIXATION_ONSET'); end
    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, results_table.stim_id(i));
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
    % Clear any events that occurred during fixation before presenting stimulus
    KbQueueFlush(p.keys.device);
    % --- present image for exactly mst_image_dur, then a responsive blank ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);
    if is_eyetracking
        Eyelink('Message', 'SYNCTIME');
        Eyelink('Message', 'STIM_ONSET %s', char(trial_info.stim_id));
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', char(img_path), p.xCenter, p.yCenter);
    end

    key_pressed = "NA";
    response_time = NaN;
    responded = false;
    blank_flipped = false;
    blank_time = stim_onset_time + p.timing.mst_image_dur - 0.5 * p.ifi;
    resp_window_end = stim_onset_time + p.timing.mst_image_dur + p.timing.mst_blank_dur;
    % Single response window spanning image + blank (2 s). The image is
    % replaced by a blank exactly at mst_image_dur; responses are collected
    % throughout both phases.
    while GetSecs < resp_window_end
        % Flip to blank once the image duration has elapsed (frame-accurate).
        if ~blank_flipped && GetSecs >= blank_time
            Screen('FillRect', p.window, p.colors.bgcolor);
            Screen('Flip', p.window, blank_time);
            Screen('Close', img_texture);
            if is_eyetracking, Eyelink('Message', 'BLANK_ONSET'); end
            blank_flipped = true;
        end
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        if pressed && ~responded
            responded = true;
            response_key_code = find(firstPress, 1);
            response_time = firstPress(response_key_code) - stim_onset_time;
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif any(response_key_code == old_key)
                key_pressed = string(p.keys.mst_old);
            elseif any(response_key_code == similar_key)
                key_pressed = string(p.keys.mst_similar);
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
    % Safety: guarantee the texture is closed / blank shown even if the image
    % duration was shorter than one refresh (should not happen at 1.5 s).
    if ~blank_flipped
        Screen('FillRect', p.window, p.colors.bgcolor);
        Screen('Flip', p.window);
        Screen('Close', img_texture);
        if is_eyetracking, Eyelink('Message', 'BLANK_ONSET'); end
    end
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
end % end of the trial loop

% --- Tail fixation before the run ends (mirrors the lead-in) ---
if isfield(p.timing, 'block_tail'), tail_dur = p.timing.block_tail; else, tail_dur = 6; end
draw_fixation_target(p);
tail_onset = Screen('Flip', p.window);
if is_eyetracking, Eyelink('Message', 'TAIL_ONSET'); end
WaitSecs('UntilTime', tail_onset + tail_dur);

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
