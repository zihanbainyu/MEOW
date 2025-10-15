%==========================================================================
%                  Phase 2: 2-Back Test Task
%==========================================================================
% author: Zihan Bai

function [results_table] = D_run_retrieval(p, test_schedule)

    %% PREPARE FOR THE BLOCK
    %----------------------------------------------------------------------
   
    % make a copy of the schedule to fill in with results
    results_table = test_schedule;
      
    % add new columns to this table to store participant responses
    num_trials = height(results_table);
    results_table.response_key = strings(num_trials, 1);
    results_table.response_key(:) = "NA";
    results_table.correct = nan(num_trials, 1);
    results_table.rt = nan(num_trials, 1);

    % define key names using KbName for robustness
    KbName('UnifyKeyNames');
    escape_key = KbName(p.keys.quit);
    same_key = KbName(p.keys.same);
    diff_key = KbName(p.keys.diff);
    space_key = KbName('space');
    
    % define the number of trials per block before a break
    trials_per_block = 30; % More frequent breaks for a longer, harder task
    
    %% DISPLAY "GET READY" SCREEN
    %----------------------------------------------------------------------
    DrawFormattedText(p.window, 'Phase 2: 2-Back Task\n\nRemember, you are now comparing each picture to the one from TWO trials ago.\n\nPlease press SPACE to start.', 'center', 'center', p.colors.white);
    Screen('Flip', p.window);
    
    % wait for space
    while true
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(space_key)
                break;
            elseif keyCode(escape_key)
                sca;
                return;
            end
        end
    end
    
    % initial fixation before the first trial
    Screen('TextSize', p.window, 100); % make the fixation cross larger
    DrawFormattedText(p.window, '+', 'center', 'center', p.colors.white);
    block_start_time = Screen('Flip', p.window);
    WaitSecs(3); % initial fixation for 3 second
    
    %% MAIN TRIAL LOOP
    %----------------------------------------------------------------------
    for i = 1:30
        %num_trials
        % --- check for a self-paced break ---
        if mod(i, trials_per_block) == 1 && i > 1
            DrawFormattedText(p.window, 'You may now take a short break.\n\nPress SPACE when you are ready to continue.', ...
                'center', 'center', p.colors.white);
            Screen('Flip', p.window);
            KbWait([], 2); 
            while true
                [~, ~, keyCode] = KbCheck;
                if keyCode(space_key), break; end
            end
            DrawFormattedText(p.window, 'Get Ready!', 'center', 'center', p.colors.white);
            Screen('Flip', p.window);
            WaitSecs(1.5);
        end
        
        % --- Step 1: Fixation ---
        DrawFormattedText(p.window, '+', 'center', 'center', p.colors.white);
        fix_onset_time = Screen('Flip', p.window);
        
        % --- Step 2: Prepare the image stimulus ---
        img_path = fullfile(p.stim_dir, results_table.stimulus_id(i));
        if ~exist(img_path, 'file')
            error('cannot find image file: %s', img_path);
        end
        
        img_data = imread(img_path);
        img_texture = Screen('MakeTexture', p.window, img_data);
        Screen('DrawTexture', p.window, img_texture, [], [], 0);
        
        % --- Step 3: Present image and collect response ---
        stim_onset_time = Screen('Flip', p.window, fix_onset_time + p.timing.fix_dur - 0.5 * p.ifi);
        
        key_pressed = "NA";
        response_time = NaN;
        responded = false;
        
        while GetSecs < stim_onset_time + p.timing.image_dur
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown && ~responded
                responded = true;
                response_time = secs - stim_onset_time;
                if keyCode(escape_key)
                    error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
                elseif keyCode(same_key)
                    key_pressed = string(p.keys.same);
                elseif keyCode(diff_key)
                    key_pressed = string(p.keys.diff);
                else
                    key_pressed = "invalid";
                end
            end
        end
        
        Screen('Close', img_texture);
        
        %% RECORD RESULTS FOR THIS TRIAL
        %------------------------------------------------------------------
        results_table.response_key(i) = key_pressed;
        results_table.rt(i) = response_time;
        
        % --- Score the response for accuracy ---
        is_correct = 0;
        
        % get the correct response for this trial from the schedule.
        % we use trial_type_final because the shuffling can create
        % incidental repeats that must be scored correctly.
        correct_resp_type = string(results_table.trial_type_final(i));
        
        % determine the correct key based on the final trial type
        correct_key_for_trial = "NA"; % Default to no response
        if correct_resp_type == "repeat"
            correct_key_for_trial = string(p.keys.same);
        elseif correct_resp_type == "lure"
            correct_key_for_trial = string(p.keys.diff);
        end
        
        % compare participant's response to what was expected
        if key_pressed == correct_key_for_trial
            is_correct = 1;
        % Special case: a correct 'none' response where a 'none' response was not designed
        % e.g. for a trial_type_designed == 'new', the correct response is 'none'.
        elseif correct_resp_type == "new" && key_pressed == "NA"
             is_correct = 1;
        end
        
        results_table.correct(i) = is_correct;
        

    end % end of the main trial loop
    
    % reset screen text size for instructions
    Screen('TextSize', p.window, 36);
    
    fprintf('Phase 2 complete.\n');

end % end of the function