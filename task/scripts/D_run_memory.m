%==========================================================================
%                  Phase 2: 2-Back Retrieval Task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================

function [results_table] = D_run_memory(p, test_schedule_block)
    %% ========================================================================
    %  SECTION 1: SET UP
    %  ========================================================================
    
    % receives a schedule for only one block.
    results_table = test_schedule_block;
      
    % add new columns to this table to store participant responses
    num_trials_in_block = height(results_table);
    results_table.response_key = strings(num_trials_in_block, 1);
    results_table.response_key(:) = "NA";
    results_table.correct = nan(num_trials_in_block, 1);
    results_table.rt = nan(num_trials_in_block, 1);

    % define key names
    KbName('UnifyKeyNames');
    escape_key = KbName(p.keys.quit);
    same_key = KbName(p.keys.same);
    diff_key = KbName(p.keys.diff);
    space_key = KbName('space');
    
    % get current and total block number for instructions
    current_block = results_table.block(1);
    total_blocks = p.nBlocks;
    
    %% ========================================================================
    %  SECTION 2: BLOCK & TRIAL EXECUTION
    %  ========================================================================

    %------------------------------------------------------------------
    % 2A: Start of Block Screen
    %------------------------------------------------------------------
    block_instruction = sprintf('Phase 2: Block %d of %d.\n\nPress SPACE to begin.', current_block, total_blocks);
    DrawFormattedText(p.window, block_instruction, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    
    while true
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(space_key), break;
            elseif keyCode(escape_key), error('USER_ABORT'); end
        end
    end

    % initial fixation before the first trial of the block
    DrawFormattedText(p.window, '+', 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    WaitSecs(2);

    %------------------------------------------------------------------
    % 2B: Trial loop for the current block
    %------------------------------------------------------------------
    for i = 1:height(results_table)
        
        if p.eyetracking == 1, Eyelink('Message', 'TRIALID %d', i); end

        % ---------fixation------------
        DrawFormattedText(p.window, '+', 'center', 'center', p.colors.black);
        fix_onset_time = Screen('Flip', p.window);
        if p.eyetracking == 1, Eyelink('Message', 'FIXATION_ONSET'); end

        % ---------stimulus------------
        img_path = fullfile(p.stim_dir, results_table.stimulus_id(i));
        if ~exist(img_path, 'file'), error('cannot find image file: %s', img_path); end
        img_data = imread(img_path);
        img_texture = Screen('MakeTexture', p.window, img_data);
        Screen('DrawTexture', p.window, img_texture, [], [], 0);

        % --- present image and collect response ---
        % Use the pre-generated JITTERED fixation duration for this trial
        stim_onset_time = Screen('Flip', p.window, fix_onset_time + results_table.fix_duration(i) - 0.5 * p.ifi);
        if p.eyetracking == 1, Eyelink('Message', 'STIM_ONSET %s', results_table.stimulus_id{i}); end

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
                elseif keyCode(same_key), key_pressed = string(p.keys.same);
                elseif keyCode(diff_key), key_pressed = string(p.keys.diff);
                else, key_pressed = "invalid";
                end
                
                if p.eyetracking == 1, Eyelink('Message', 'RESPONSE KEY %s RT %.0f', key_pressed, response_time * 1000); end
            end
        end
        Screen('Close', img_texture);

        %------------------------------------------------------------------
        % 2C: Record trial data
        %------------------------------------------------------------------
        results_table.response_key(i) = key_pressed;
        results_table.rt(i) = response_time;
        
        correct_resp = string(results_table.correct_response(i));
        if (key_pressed == "NA" && correct_resp == "none") || (key_pressed == correct_resp)
            results_table.correct(i) = 1;
        else
            results_table.correct(i) = 0;
        end

        if p.eyetracking == 1, Eyelink('Message', 'TRIAL_END %d', i); end

    end % end of trial loop

    %% ========================================================================
    %  SECTION 3: SAVE BLOCK DATA
    %  ========================================================================
    try
        block_filename = sprintf('sub%03d_mem_b%d.mat', p.subj_id, current_block);
        block_filepath = fullfile(p.results_dir, block_filename);
        
        save(block_filepath, 'results_table');
        fprintf('Retrieval block %d data saved.\n', current_block);
    catch ME
        warning('SAVE_FAILED: Could not save retrieval data for block %d. Reason: %s', current_block, ME.message);
    end

end % end of function