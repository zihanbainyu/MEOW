function [practice_results] = C_run_encoding_practice(p, practice_schedule)
%==========================================================================
%                  Practice: 1-Back Encoding Task
%==========================================================================
% Author: Zihan Bai
    %% Initialize results table
    practice_results = practice_schedule;
    nTrials = height(practice_results);
    
    practice_results.response_key = strings(nTrials,1);
    practice_results.correct = nan(nTrials,1);
    practice_results.rt = nan(nTrials,1);
    
    % Keyboard keys (start_key replaces space_key for consistency)
    start_key  = KbName('g');
    escape_key = KbName(p.keys.quit);
    same_key   = KbName(p.keys.same);
    diff_key   = KbName(p.keys.diff);
    
    % Apply current window settings (assuming window is open in main)
    Screen('TextSize', p.window, p.text_size);
    Screen('TextFont', p.window, 'Helvetica');
    
    % Setup KbQueue for responses
    responseKeys = zeros(1, 256);
    responseKeys([same_key, diff_key, escape_key]) = 1;
    KbQueueCreate(p.keys.device, responseKeys); % Create queue for the target device

    %% Instructions
    instruction_text = ['Practice Round: Part A (1-Back)\n\n' ...
        '\n\n' ...
        'When you are ready, press g to begin.\n\n'];
        
    DrawFormattedText(p.window, instruction_text, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    
    % Wait for G to start the practice block
    KbReleaseWait(p.keys.device);
    while true
        [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
        if keyIsDown
            if keyCode(start_key), break; % Use start_key
            elseif keyCode(escape_key), error('USER_ABORT'); end
        end
        WaitSecs(0.001); 
    end
    
    %% Initial fixation
    xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
    yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
    allCoords = [xCoords; yCoords];
    Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
    Screen('Flip', p.window);
    WaitSecs(2);
    
    %% Trial loop
    KbQueueStart(p.keys.device); % Start listening before the trial loop
    
    for i = 1:nTrials
        
        % Clear queue from any previous instruction key presses
        KbQueueFlush(p.keys.device); 
        
        % Load path
        img_path = fullfile(p.stim_dir, 'practice', practice_results.stimulus_id{i});
        if ~exist(img_path, 'file')
            error('Cannot find image: %s', img_path);
        end
    
        % Fixation
        Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
        fix_onset = Screen('Flip', p.window);
        
        % Load & draw image
        img_data = imread(img_path);
        img_tex = Screen('MakeTexture', p.window, img_data);
        Screen('DrawTexture', p.window, img_tex, [], [], 0);
        
        % Stimulus onset (wait for fixation duration)
        stim_onset = Screen('Flip', p.window, fix_onset + practice_results.fix_duration(i) - 0.5 * p.ifi);
    
        % Collect response (using KbQueueCheck for accuracy)
        key_pressed = "NA";
        rt = NaN;
        responded = false;
        
        while GetSecs < stim_onset + p.timing.image_dur
            [pressed, firstPress] = KbQueueCheck(p.keys.device);
            
            if pressed && ~responded
                responded = true;
                
                response_key_code = find(firstPress, 1);
                rt = firstPress(response_key_code) - stim_onset;
                
                if response_key_code == escape_key
                    error('USER_ABORT:ExperimentAborted','Experiment aborted by user.');
                elseif response_key_code == same_key
                    key_pressed = string(p.keys.same);
                elseif response_key_code == diff_key
                    key_pressed = string(p.keys.diff);
                else
                    key_pressed = "invalid";
                end
            end
        end
        Screen('Close', img_tex);
    
        % Record
        practice_results.response_key(i) = key_pressed;
        practice_results.rt(i) = rt;
    
        correct_resp = string(practice_results.correct_response(i));
        if (key_pressed == "NA" && correct_resp == "none") || (key_pressed == correct_resp)
            practice_results.correct(i) = 1;
        else
            practice_results.correct(i) = 0;
        end
    end
    
    % Release resources after the loop finishes
    KbQueueRelease(p.keys.device);
    
    %% Save practice results
    try
        filename = sprintf('sub%03d_enc_practice.mat', p.subj_id);
        filepath = fullfile(p.results_dir, filename);
        save(filepath, 'practice_results');
        fprintf('Practice data saved: %d trials.\n', nTrials);
    catch ME
        warning('Could not save practice data: %s', ME.message);
    end
    
    %% End of practice screen
    DrawFormattedText(p.window, ...
        ['Great! You have completed the practice round.\n\n' ...
         'If you have any questions, please find the experimenter.'], ...
        'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    
    % Wait for G to exit practice screen
    KbReleaseWait(p.keys.device);
    while true
        [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
        if keyIsDown
            if keyCode(start_key)
                break; % exit loop
            elseif keyCode(escape_key)
                error('USER_ABORT');
            end
        end
        WaitSecs(0.001);
    end
    return
end