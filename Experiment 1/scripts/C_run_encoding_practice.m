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
    
    % Keyboard keys
    KbName('UnifyKeyNames');
    escape_key = KbName(p.keys.quit);
    same_key   = KbName(p.keys.same);
    diff_key   = KbName(p.keys.diff);
    space_key  = KbName('g');
    
    %% Instructions
    DrawFormattedText(p.window, ...
        ['Practice Round: 1-Back Task\n\n' ...
        '\n\n'...
        'You will view a series of items, one after another\n\n' ...
        'Your task is to compare each item to the one that came right before it\n\n' ...
        '\n\n'...
        'Press j (repeat) if the item is exactly the same as the one just before it\n\n' ...
        'Press k (similar) if the item is similar but not identical as the one just before it\n\n' ...
        'Do not press any key if the item is completely new\n\n' ...
        '\n\n' ...
        'You are always comparing the current image to the one just before it'
        'When you are ready, press g to begin'], ...
        'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    
    while true
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(space_key), break;
            elseif keyCode(escape_key), error('USER_ABORT'); end
        end
    end
    
    %% Initial fixation
    xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
    yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
    allCoords = [xCoords; yCoords];
    Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
    Screen('Flip', p.window);
    WaitSecs(2);
    
    %% Trial loop
    for i = 1:nTrials
        img_path = fullfile(p.stim_dir, 'practice', practice_results.stimulus_id{i});
        if ~exist(img_path, 'file')
            error('Cannot find image: %s', img_path);
        end
    
        % Fixation
        Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
        fix_onset = Screen('Flip', p.window);
        WaitSecs(practice_results.fix_duration(i));
    
        % Load & draw image
        img_data = imread(img_path);
        img_tex = Screen('MakeTexture', p.window, img_data);
        Screen('DrawTexture', p.window, img_tex, [], [], 0);
        stim_onset = Screen('Flip', p.window);
    
        % Collect response
        key_pressed = "NA";
        rt = NaN;
        responded = false;
        while GetSecs < stim_onset + p.timing.image_dur
            [keyIsDown, t, keyCode] = KbCheck;
            if keyIsDown && ~responded
                responded = true;
                rt = t - stim_onset;
    
                if keyCode(escape_key)
                    error('USER_ABORT:ExperimentAborted','Experiment aborted by user.');
                elseif keyCode(same_key)
                    key_pressed = string(p.keys.same);
                elseif keyCode(diff_key)
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
    
    %% Save practice results
    try
        filename = sprintf('sub%03d_enc_practice.mat', p.subj_id);
        filepath = fullfile(p.results_dir, filename);
        save(filepath, 'practice_results');
        fprintf('Practice data saved: %d trials.\n', nTrials);
    catch ME
        warning('Could not save practice data: %s', ME.message);
    end
    
    Screen('TextSize', p.window, 32);
    DrawFormattedText(p.window, ...
        ['Great! You have completed the practice round\n\n' ...
         'If you have any questions, please find the experimenter'], ...
        'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    
    wait_for_end = true;
    while wait_for_end
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(space_key)
                wait_for_end = false; % exit loop
            elseif keyCode(escape_key)
                error('USER_ABORT');
            end
        end
    end
    return

end
