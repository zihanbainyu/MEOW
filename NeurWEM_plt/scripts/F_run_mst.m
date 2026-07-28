%==========================================================================
%              Post-task MST: old / similar / new recognition
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%
% Standard MST on the 1-back repeat foils (built in A_subject_setup.m, P6B):
%   old  = the exact repeat item seen in the 1-back            -> 'j'
%   lure = its similar pairmate, never shown                   -> 'k'
%   new  = a fresh foil, never shown                           -> 'l'
%
% Fixed 1.5 s image with the response collected DURING the image (fixed
% window, not self-paced), jittered fixation ITI, and a 10 s lead-in / tail
% like the n-back runs.
%==========================================================================
function [results_table] = F_run_mst(p, sequence_mst)

%% SET UP
results_table = sequence_mst;
num_trials = height(results_table);
results_table.resp_key = strings(num_trials, 1);
results_table.resp_key(:) = "NA";
results_table.rt = nan(num_trials, 1);

% keys
old_key     = KbName(p.keys.mst_old);       % 'j' = OLD
similar_key = KbName(p.keys.mst_similar);   % 'k' = SIMILAR
new_key     = KbName(p.keys.mst_new);       % 'l' = NEW
escape_key  = KbName(p.keys.quit);
start_key   = KbName('f');

Screen('TextSize', p.window, p.text_size);
Screen('TextFont', p.window, 'Helvetica');

respKeys = zeros(1, 256);
respKeys([old_key, similar_key, new_key, escape_key]) = 1;
KbQueueCreate(p.keys.device, respKeys);
KbCheck(p.keys.device);

%% INSTRUCTION
% Text placeholder. Once an instruction image exists, replace this whole
% block with:  run_instructions(p, {'ins_mst'});   (as done for recognition)
DrawFormattedText(p.window, ...
    ['Final Memory Test\n\n\n' ...
    'You will see objects one at a time.\n\n' ...
    'For each object, decide:\n\n' ...
    'j = OLD (the exact object you saw earlier)\n' ...
    'k = SIMILAR (the same object, but changed)\n' ...
    'l = NEW (an object you did not see)\n\n\n' ...
    'Respond while the image is on the screen.\n\n' ...
    'When you are ready, press f to begin.'], ...
    'center', 'center', p.colors.black, [], [], [], 1.2);
Screen('Flip', p.window);
KbReleaseWait(p.keys.device);
while true
    [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
    if keyIsDown
        if keyCode(start_key), break;
        elseif keyCode(escape_key), error('USER_ABORT'); end
    end
    WaitSecs(0.001);
end

% Lead-in fixation
draw_fixation_target(p);
lead_in_onset = Screen('Flip', p.window);
WaitSecs('UntilTime', lead_in_onset + p.timing.block_lead_in);

KbQueueStart(p.keys.device);

%% TRIAL LOOP
for i = 1:num_trials
    trial_info = results_table(i,:);

    % --------- fixation ------------
    draw_fixation_target(p);
    fix_onset_time = Screen('Flip', p.window);

    % --------- stimulus preparation ------------
    img_path = fullfile(p.stim_dir, trial_info.stim_id);
    if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
    img_data = imread(img_path);
    img_texture = Screen('MakeTexture', p.window, img_data);
    Screen('DrawTexture', p.window, img_texture, [], [], 0);
    KbQueueFlush(p.keys.device);

    % --- present image (fixed 1.5 s) and collect response during it ---
    stim_onset_time = Screen('Flip', p.window, fix_onset_time + trial_info.fix_duration - 0.5 * p.ifi);

    key_pressed = "NA";
    response_time = NaN;
    responded = false;
    while GetSecs < stim_onset_time + p.timing.image_dur
        [pressed, firstPress] = KbQueueCheck(p.keys.device);
        if pressed && ~responded
            response_key_code = find(firstPress > 0, 1);
            if isempty(response_key_code), continue; end
            response_time = firstPress(response_key_code) - stim_onset_time;
            if response_key_code == escape_key
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            elseif response_key_code == old_key
                key_pressed = string(p.keys.mst_old);     responded = true;
            elseif response_key_code == similar_key
                key_pressed = string(p.keys.mst_similar); responded = true;
            elseif response_key_code == new_key
                key_pressed = string(p.keys.mst_new);     responded = true;
            end
        end
    end
    Screen('Close', img_texture);

    % --------- record ------------
    results_table.resp_key(i) = key_pressed;
    results_table.rt(i) = response_time;
end

%% TAIL
if isfield(p.timing, 'block_tail'), tail_dur = p.timing.block_tail; else, tail_dur = 10; end
draw_fixation_target(p);
tail_onset = Screen('Flip', p.window);
WaitSecs('UntilTime', tail_onset + tail_dur);

Screen('FillRect', p.window, p.colors.bgcolor);
Screen('Flip', p.window);
WaitSecs(0.05);
KbQueueRelease(p.keys.device);
end

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function draw_fixation_target(p)
% Thaler, Schutz, Goodale & Gegenfurtner (2013, Vision Research): combined
% bullseye-and-crosshair target. Same target used across the run scripts.
d1 = p.fix_dot_d1;   % outer disc diameter (px)
d2 = p.fix_dot_d2;   % central dot diameter / crosshair width (px)
cx = p.xCenter; cy = p.yCenter;
Screen('FillOval', p.window, p.fix_dot_color, [cx-d1/2, cy-d1/2, cx+d1/2, cy+d1/2]);
Screen('FillRect', p.window, p.colors.bgcolor, [cx-d1/2, cy-d2/2, cx+d1/2, cy+d2/2]);
Screen('FillRect', p.window, p.colors.bgcolor, [cx-d2/2, cy-d1/2, cx+d2/2, cy+d1/2]);
Screen('FillOval', p.window, p.fix_dot_color, [cx-d2/2, cy-d2/2, cx+d2/2, cy+d2/2]);
end
