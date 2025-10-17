%==========================================================================
%                  Phase 1: 1-Back Encoding Task
%==========================================================================
% author: Zihan Bai
%==========================================================================

function [results_table] = run_encoding(p, encoding_schedule)  
    %% Set up for the block
    % make a copy of the schedule to fill in with results
    results_table = encoding_schedule;
      
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
    trials_per_block = 30; % CHANGE THIS AFTER TESTING
    
    %% Reminder Screen
    DrawFormattedText(p.window, ['Part 1: 1-back Encoding\n\n' ...
        'You will now start the actual task.\n\n' ...
        'Again, it is very important that you respond as accurately and fast as you can.\n\n' ...
        'Press place your fingers on appropriate keys and press SPACE when you are ready.'], 'center', 'center', p.colors.white);
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
    block_start_time = Screen('Flip', p.window); % use this to time the first trial
    WaitSecs(2); % initial fixation for 3 second
    
    %% Main trial loop
    for i = 1:20 % CHANGE THIS AFTER TESTING TO num_trials
        % num_trials
        
        %------------------------------------------------------------------
        % Eyelink message: start of trial
        %------------------------------------------------------------------
        if p.eyetracking == 1
            Eyelink('Message', 'TRIALID %d COND %s ROLE %s', ...
                results_table.trial_phase(i), ...
                results_table.condition{i}, ...
                results_table.role{i});
        end

        %------------------------------------------------------------------
        % Check for a self-paced break
        %------------------------------------------------------------------
        % we also make sure not to break before the very first trial
        if mod(i, trials_per_block) == 1 && i > 1
            if p.eyetracking == 1
                Eyelink('Message', 'DISPLAY_BREAK_SCREEN');
            end
            DrawFormattedText(p.window, 'You may now take a short break.\n\nPlease press SPACE when you are ready to continue.', 'center', 'center', p.colors.white);
            Screen('Flip', p.window);
            while true
                [~, ~, keyCode] = KbCheck;
                if keyCode(space_key), break; end
            end
            if p.eyetracking == 1
                Eyelink('Message', 'END_BREAK_SCREEN'); 
            end
            DrawFormattedText(p.window, 'Get Ready!', 'center', 'center', p.colors.white);
            Screen('Flip', p.window);
            WaitSecs(1.5);
        end
        
        %------------------------------------------------------------------
        % Fixation
        %------------------------------------------------------------------
        DrawFormattedText(p.window, '+', 'center', 'center', p.colors.white);
        fix_onset_time = Screen('Flip', p.window);
        if p.eyetracking == 1
            Eyelink('Message', 'FIXATION_ONSET'); 
        end
        
        %------------------------------------------------------------------
        % Pepare the image stimulus
        %------------------------------------------------------------------
        img_path = fullfile(p.stim_dir, results_table.stimulus_id(i));

        if ~exist(img_path, 'file') % check if image exists to prevent a crash
            error('cannot find image file: %s', img_path);
        end
        
        img_data = imread(img_path);
        img_texture = Screen('MakeTexture', p.window, img_data);
        
        % draw the texture to the back buffer (it won't be shown yet)
        Screen('DrawTexture', p.window, img_texture, [], [], 0);
        
        % present image and collect response ---
        % flip the screen at the precise moment the fixation duration is over.
        % the return value 'stim_onset_time' is our critical timestamp (t=0) for this trial.
        stim_onset_time = Screen('Flip', p.window, fix_onset_time + p.timing.fix_dur - 0.5 * p.ifi);
        
        %------------------------------------------------------------------
        % Eyelink message: stimulus onset
        %------------------------------------------------------------------
        if p.eyetracking == 1
            Eyelink('Message', 'SYNCTIME'); % high-precision sync event
            Eyelink('Message', 'STIM_ONSET %s', results_table.stimulus_id{i});
        end

        % prepare for response collection
        key_pressed = "NA";
        response_time = NaN;
        responded = false;
        
        % loop for the entire duration the image is on screen
        while GetSecs < stim_onset_time + p.timing.image_dur
            
            % check the keyboard
            [keyIsDown, secs, keyCode] = KbCheck;
            
            if keyIsDown && ~responded % only log the first keypress
                responded = true; % prevent logging multiple presses
                
                % calculate reaction time
                response_time = secs - stim_onset_time;
                
                % figure out which key was pressed
                if keyCode(escape_key)
                    % if escape is pressed, save data and gracefully exit
                    error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
                elseif keyCode(same_key)
                    key_pressed = string(p.keys.same);
                elseif keyCode(diff_key)
                    key_pressed = string(p.keys.diff);
                else
                    key_pressed = "invalid"; % pressed a wrong key
                end

                % --- EYELINK MESSAGE: RESPONSE ---
                if p.eyetracking
                    Eyelink('Message', 'RESPONSE KEY %s RT %.0f', key_pressed, response_time * 1000);
                end
            end
        end
        
        % we are done with this trial's image texture, so we close it
        % this frees up video memory
        Screen('Close', img_texture);
        
        %% Record results for each trial
        % store the collected response and RT in our results table
        results_table.response_key(i) = key_pressed;
        results_table.rt(i) = response_time;
        
        %------------------------------------------------------------------
        % Accuracy
        %------------------------------------------------------------------
        is_correct = 0; % default to incorrect
        
        % get the correct response from our schedule (e.g., "j", "k", or "none")
        correct_resp = string(results_table.correct_response(i));
        
        % compare participant's response to the correct one
        if key_pressed == "NA" && correct_resp == "none"
            is_correct = 1;
        elseif key_pressed == correct_resp
            is_correct = 1;
        end
        
        results_table.correct(i) = is_correct;

        %------------------------------------------------------------------
        % Eyelink message: end of trial
        %------------------------------------------------------------------
        % (send this just before the next trial's fixation appears)
        if p.eyetracking == 1
            Eyelink('Message', 'TRIAL_END %d', results_table.trial_phase(i));
        end
        
    end % end of the main trial loop
    
    % reset screen text size for instructions
    Screen('TextSize', p.window, 26);
    Screen('TextFont', p.window, 'Open Sans');
    
    DrawFormattedText(p.window, ['Excellent. You have now completed Part 1.\n\n'...
        'Please press SPACE to close this screen'], 'center', 'center', p.colors.white);
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

end % end of the function