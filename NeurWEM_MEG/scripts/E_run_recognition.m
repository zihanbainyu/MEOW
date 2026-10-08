%==========================================================================
%              Part 3: Final recognition test (old / new)
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%
% Out-of-scanner version, run on the laptop: no scanner trigger, no
% eyetracking, no button box. The participant judges each object OLD or NEW
% with the laptop keyboard. Each trial is a short fixation (fix_duration,
% 0.5 s) followed by the image, which stays on screen until the participant
% responds (self-paced, with a safety timeout). Both OLD and NEW require a
% keypress. A brief blank inter-trial interval follows each response.
%
% Response mapping (matches the n-back keys stored in the setup):
%   p.keys.same ('1') = OLD   (corr_resp "1")
%   p.keys.diff ('2') = NEW   (corr_resp "2")
%==========================================================================
function [results_table] = E_run_recognition(p, sequence_recognition)
%% ========================================================================
% SECTION 1: SET UP
% ========================================================================
results_table = sequence_recognition;
results_table.corr_resp = string(results_table.corr_resp);  % setup stores this as char; match the MST table's string type
num_trials = height(results_table);
results_table.resp_key = strings(num_trials, 1);
results_table.resp_key(:) = "NA";
results_table.rt = nan(num_trials, 1);

% define key names (assumes KbName('UnifyKeyNames') was called by the launcher)
old_key    = KbName({'1!','1'});   % OLD: top-row '1' and numpad '1'
new_key    = KbName({'2@','2'});   % NEW: top-row '2' and numpad '2'
escape_key = KbName(p.keys.quit);
start_key  = KbName('space');

% response window: image stays up until a response, or up to max_resp (4 s),
% after which the trial is scored no-response ("NA") and advances. Then a blank ITI.
if isfield(p.timing, 'rec_max_resp'), max_resp = p.timing.rec_max_resp; else, max_resp = 4; end
if isfield(p.timing, 'rec_isi'),      isi      = p.timing.rec_isi;      else, isi      = 0.3; end

Screen('TextSize', p.window, p.text_size);
Screen('TextFont', p.window, 'Helvetica');

% Initialize KbQueue for fast, buffered responses during the trial loop
respKeys = zeros(1, 256);
respKeys([old_key, new_key, escape_key]) = 1;
KbQueueCreate(p.keys.device, respKeys);
KbCheck(p.keys.device);

%% ========================================================================
% SECTION 2: RUN EXECUTION
% ========================================================================
%------------------------------------------------------------------
% 2A: Instruction / start screen -- participant advances with SPACE
%------------------------------------------------------------------
instr = ['Final memory test\n\n\n' ...
    'You will see one object at a time.\n\n' ...
    'For each object, decide whether you saw it earlier in the experiment.\n\n\n' ...
    'Press  1  if the object is OLD  (you saw it before).\n\n' ...
    'Press  2  if the object is NEW  (you did not see it before).\n\n\n' ...
    'You have up to ' num2str(max_resp) ' seconds to respond on each object.\n\n' ...
    'Answer as quickly and accurately as you can.\n\n\n' ...
    'Press the SPACE bar to begin.'];
DrawFormattedText(p.window, instr, 'center', 'center', p.colors.black, [], [], [], 1.4);
Screen('Flip', p.window);
KbReleaseWait(p.keys.device);
while true
    [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
    if keyIsDown
        if keyCode(start_key)
            break;
        elseif keyCode(escape_key)
            error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
        end
    end
    WaitSecs(0.001);
end

KbQueueStart(p.keys.device);

for i = 1:num_trials
    trial_info = results_table(i,:);
    %------------------------------------------------------------------
    % 2B: Trial loop
    %------------------------------------------------------------------
    % --------- fixation ------------
    draw_fixation_target(p);
    fix_onset_time = Screen('Flip', p.window);
    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, results_table.stim_id(i));
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
    % Clear any events that occurred during fixation before presenting stimulus
    KbQueueFlush(p.keys.device);
    % --- present image; it stays on screen until a response (or timeout) ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);

    key_pressed = "NA";
    response_time = NaN;
    while GetSecs < stim_onset_time + max_resp
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        if pressed
            response_key_code = find(firstPress, 1);
            response_time = firstPress(response_key_code) - stim_onset_time;
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif any(response_key_code == old_key)
                key_pressed = string(p.keys.same);   % "1" = OLD
                break;
            elseif any(response_key_code == new_key)
                key_pressed = string(p.keys.diff);   % "2" = NEW
                break;
            end
        end
        WaitSecs(0.001);
    end
    Screen('Close', img_texture);

    % --------- blank inter-trial interval ------------
    Screen('FillRect', p.window, p.colors.bgcolor);
    Screen('Flip', p.window);
    WaitSecs(isi);

    %------------------------------------------------------------------
    % 2C: Record trial data
    %------------------------------------------------------------------
    results_table.resp_key(i) = key_pressed;
    results_table.rt(i) = response_time;
end % end of the trial loop

% --- End screen ---
DrawFormattedText(p.window, 'You have completed the memory test.\n\nThank you!', ...
    'center', 'center', p.colors.black, [], [], [], 1.4);
Screen('Flip', p.window);
WaitSecs(2);

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
