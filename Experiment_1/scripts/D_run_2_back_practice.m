function D_run_2_back_practice(p)
%==========================================================================
%                  Practice: 2-Back (Minimal - No Save)
%==========================================================================
% Author: Zihan Bai
    
    % --- 1. Generate the Practice Schedule In-Line ---
    sequence_2_back_practice = gen_2_back_practice(p); % 
    
    nTrials = height(sequence_2_back_practice);
    
    % Keyboard keys
    start_key  = KbName('g');
    escape_key = KbName(p.keys.quit);
    same_key   = KbName(p.keys.same);
    similar_key   = KbName(p.keys.diff);
    
    % Apply current window settings
    Screen('TextSize', p.window, p.text_size);
    Screen('TextFont', p.window, 'Helvetica');
    
    % Setup KbQueue for responses (Still needed for the ESC key check)
    responseKeys = zeros(1, 256);
    responseKeys([same_key, similar_key, escape_key]) = 1;
    KbQueueCreate(p.keys.device, responseKeys);
    
    %% 2. Instructions & Start
    DrawFormattedText(p.window, ['Practice Round: 2-Back' ...
        '\n\n\n' ...
        'When you are ready, press g to begin.\n\n'], 'center', 'center', p.colors.black);
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
    
    %% 3. Initial fixation
    xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
    yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
    allCoords = [xCoords; yCoords];
    Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
    Screen('Flip', p.window);
    WaitSecs(2);
    
    %% 4. Trial loop
    KbQueueStart(p.keys.device);
    
    for i = 1:nTrials
        
        KbQueueFlush(p.keys.device); 
        
        img_path = fullfile(p.stim_dir, sequence_2_back_practice.stimulus_id{i});
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
        stim_onset = Screen('Flip', p.window, fix_onset + sequence_2_back_practice.fix_duration(i) - 0.5 * p.ifi);
    
        % Collect response (only checking for escape key, otherwise ignoring response)
        while GetSecs < stim_onset + p.timing.image_dur
            [pressed, firstPress] = KbQueueCheck(p.keys.device);
            
            if pressed
                response_key_code = find(firstPress, 1);
                
                if response_key_code == escape_key
                    error('USER_ABORT:ExperimentAborted','Experiment aborted by user.');
                end
                % Responses (same_key/diff_key) are ignored as results are not saved
            end
        end
        Screen('Close', img_tex);
    
    end
    
    % Release resources
    KbQueueRelease(p.keys.device);
    
    %% 5. End of practice screen
    DrawFormattedText(p.window, ...
        ['Great! You have completed the practice round.\n\n' ...
         'If you have any questions, please find the experimenter.\n\n' ...
         'Otherwise, press f to start the actual task!'], ...
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

    % Clear screen before returning to main
    Screen('FillRect', p.window, p.colors.bgcolor);
    Screen('Flip', p.window);
    
    % Function must not return practice_results, as it is now a procedure-only function
    return
end

%==========================================================================
%                  Sub-Function: Generate Practice Schedule (2-Back)
%==========================================================================
function [sequence_2_back_practice] = gen_2_back_practice(p)
    
    A_files = dir(fullfile(p.stim_dir,'prac_*_A.png'));
    B_files = dir(fullfile(p.stim_dir,'prac_*_B.png'));

    % Convert names to strings
    A_names = string({A_files.name});
    B_names = string({B_files.name});
    
    % Convert names to strings
    A_names = string({A_files.name});
    B_names = string({B_files.name});
    seq = { ...
        A_names(5), "foil", "A", "none"; ... 
        A_names(7), "foil", "A", "none"; ...
        A_names(6), "same", "A", "none"; ... 
        A_names(8), "similar", "A", "none"; ...
        A_names(6), "same", "A", "j";...
        B_names(8), "similar", "B", "k";   ... 
        A_names(17), "similar", "A", "none"; ...
        A_names(11), "foil", "A", "none";...
        B_names(17), "similar", "B", "k"; ...
        B_names(23), "same", "A", "none";...
        A_names(9), "foil", "A", "none";...
        B_names(23), "same", "A", "j";...
        A_names(18), "foil", "A", "none";...

    };
    
    % Only the columns required to run the task (stim_id and timing) are kept
    col_names = {'stim_id','condition','identity','corr_resp','fix_duration'};
    variable_types = {'string','string','string','string','double'}; 
    
    sequence_2_back_practice = table('Size',[0,5], 'VariableTypes', variable_types, ...
                              'VariableNames', col_names);
    
    for i = 1:size(seq,1)
        stim_id = seq{i,1};
        cond = seq{i,2};
        role = seq{i,3};
        correct_resp = seq{i,4};
        
        % Jittered fixation 0.5-1.0s
        fix_dur = 0.75 + (rand*0.5 - 0.25);
        
        % Append row
        new_row = table(stim_id, cond, role, correct_resp, fix_dur, ...
                    'VariableNames', col_names);
        sequence_2_back_practice = [sequence_2_back_practice; new_row];
    end
end