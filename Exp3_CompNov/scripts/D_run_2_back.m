%==========================================================================
%                  Phase 2: 2-Back Retrieval Task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%
% Exp3_CompNov: no eye tracking. All Eyelink calls, the fixation onset
% gate, live gaze monitoring and the restricted-viewing fixation overlay
% have been removed. The fixation target itself is kept as the inter-trial
% fixation point.
%==========================================================================
function [results_table] = D_run_2_back(p, sequence_2_back_block, b)
%% ========================================================================
%  SECTION 1: SET UP
%  ========================================================================

results_table = sequence_2_back_block;
num_trials_in_block = height(results_table);
results_table.resp_key = strings(num_trials_in_block, 1);
results_table.resp_key(:) = "NA";
results_table.rt = nan(num_trials_in_block, 1);

% Define key codes
escape_key = KbName(p.keys.quit);
same_key = KbName(p.keys.same);
similar_key = KbName(p.keys.diff);
start_key = KbName('f');

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
    [sprintf('Block %01d of %01d:  2-Back\n\n\n\n', b, total_blocks) ...
    'Please compare each image to the one shown TWO TRIALS AGO.', ...
    '\n\n\n\n' ...
    'Press ' p.keys.same ' = SAME.     Press ' p.keys.diff ' = SIMILAR.     No press = NEW.' ...
    '\n\n\n\n' ...
    'Please rest your fingers on j and k, then press f to begin.'], 'center', 'center', p.colors.black);
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

% Lead-in fixation before the first trial of the block. This replaces the
% junk trials that used to open each block (see A_subject_setup.m), giving
% the signal time to settle before the first real trial.
draw_fixation_target(p);
lead_in_onset = Screen('Flip', p.window);
WaitSecs('UntilTime', lead_in_onset + p.timing.block_lead_in);

% Start the KbQueue to collect responses during the trial loop
KbQueueStart(p.keys.device);

%------------------------------------------------------------------
% 2B: Trial loop for the current block
%------------------------------------------------------------------

for i = 1:height(results_table)

    trial_info = results_table(i,:);

    % ---------fixation------------
    draw_fixation_target(p);
    fix_onset_time = Screen('Flip', p.window);

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
            elseif response_key_code == same_key
                key_pressed = string(p.keys.same);
            elseif response_key_code == similar_key
                key_pressed = string(p.keys.diff);
            else
                key_pressed = "invalid";
            end
        end
    end
    Screen('Close', img_texture);

    %------------------------------------------------------------------
    % 2C: Record trial data
    %------------------------------------------------------------------
    results_table.resp_key(i) = key_pressed;
    results_table.rt(i) = response_time;
end % end of trial loop

% --- Clear screen at end of block ---
Screen('FillRect', p.window, p.colors.bgcolor);
Screen('Flip', p.window);
WaitSecs(0.05);

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
