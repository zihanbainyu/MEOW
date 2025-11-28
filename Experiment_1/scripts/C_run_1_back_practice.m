function C_run_1_back_practice(p)
    
    sequence_1_back_practice = gen_1_back_practice(p);
    nTrials = height(sequence_1_back_practice);
    
    % Keys
    start_key = KbName('f');
    escape_key = KbName(p.keys.quit);
    same_key = KbName(p.keys.same);
    similar_key = KbName(p.keys.diff);
    
    Screen('TextSize', p.window, p.text_size);
    Screen('TextFont', p.window, 'Helvetica');
    
    % Setup KbQueue
    respKeys = zeros(1, 256);
    respKeys([same_key, similar_key, escape_key]) = 1;
    KbQueueCreate(p.keys.device, respKeys);
    
    %% Instructions
    DrawFormattedText(p.window, ...
        ['Practice: 1-Back Task\n\n\n' ...
        'You will receive feedback after each response.\n\n' ...
        'When you are ready, press f to begin.'], ...
        'center', 'center', p.colors.black);
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
    
    %% Initial fixation
    xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
    yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
    allCoords = [xCoords; yCoords];
    Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
    Screen('Flip', p.window);
    WaitSecs(2);
    
    %% Trial loop with feedback
    KbQueueStart(p.keys.device);
    n_correct = 0;
    
    for i = 1:nTrials
        
        KbQueueFlush(p.keys.device);
        
        % Load image
        img_path = fullfile(p.stim_dir, sequence_1_back_practice.stim_id{i});
        if ~exist(img_path, 'file')
            error('Cannot find: %s', img_path);
        end
        img_data = imread(img_path);
        img_tex = Screen('MakeTexture', p.window, img_data);
        
        % Fixation
        Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
        fix_onset = Screen('Flip', p.window);
        
        % Stimulus
        Screen('DrawTexture', p.window, img_tex, [], [], 0);
        stim_onset = Screen('Flip', p.window, fix_onset + sequence_1_back_practice.fix_duration(i) - 0.5*p.ifi);
        
        % Collect response
        responded = false;
        response_key = '';
        
        while GetSecs < stim_onset + p.timing.image_dur
            [pressed, firstPress] = KbQueueCheck(p.keys.device);
            
            if pressed && ~responded
                key_code = find(firstPress, 1);
                
                if key_code == escape_key
                    error('USER_ABORT');
                elseif key_code == same_key
                    response_key = 'j';
                    responded = true;
                elseif key_code == similar_key
                    response_key = 'k';
                    responded = true;
                end
            end
        end
        
        % If no response, set to 'none'
        if ~responded
            response_key = 'none';
        end
        
        % Evaluate correctness
        correct_resp = sequence_1_back_practice.corr_resp{i};
        is_correct = strcmp(response_key, correct_resp);
        
        if is_correct
            n_correct = n_correct + 1;
        end
        
        % FEEDBACK SCREEN (1.5s)
        Screen('FillRect', p.window, p.colors.bgcolor);
        
        if is_correct
            % Green checkmark
            feedback_text = 'Correct';
            feedback_color = [0 150 0];
        else
            % Red X + correct answer
            if strcmp(correct_resp, 'j')
                correct_text = 'j (SAME)';
            elseif strcmp(correct_resp, 'k')
                correct_text = 'k (SIMILAR)';
            else
                correct_text = 'No response';
            end
            
            feedback_text = sprintf('Incorrect!\n\nCorrect response is: %s', correct_text);
            feedback_color = [200 0 0];
        end
        
        DrawFormattedText(p.window, feedback_text, 'center', 'center', feedback_color);
        Screen('Flip', p.window);
        if is_correct
            WaitSecs(0.75);  % Quick confirmation
        else
            WaitSecs(2);    % Time to read correct answer
        end
        
        Screen('Close', img_tex);
    end
    
    KbQueueRelease(p.keys.device);
    
    %% Performance summary
    accuracy = (n_correct / nTrials) * 100;
    
    if accuracy >= 100
        summary_text = sprintf(...
            ['Great! Practice complete.\n\n' ...
            'Accuracy: %.0f%%\n\n\n\n' ...
            'During the actual task, you will not receive feedback.\n\n' ...
            'When you are ready, press f to begin the actual task.'], accuracy);
        summary_color = p.colors.black;
    else
        summary_text = sprintf(...
            ['Practice complete.\n\n' ...
            'Accuracy: %.0f%%\n\n\n' ...
            'Please find the experimenter for any questions.\n\n' ...
            'Press f to try practice again.'], accuracy);
        summary_color = [200 0 0];
    end
    
    DrawFormattedText(p.window, summary_text, 'center', 'center', summary_color);
    Screen('Flip', p.window);
    
    KbReleaseWait(p.keys.device);
    while true
        [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
        if keyIsDown
            if keyCode(start_key)
                if accuracy < 80
                    % Repeat practice
                    C_run_1_back_practice(p);
                end
                break;
            elseif keyCode(escape_key)
                error('USER_ABORT');
            end
        end
        WaitSecs(0.001);
    end
    
    Screen('FillRect', p.window, p.colors.bgcolor);
    Screen('Flip', p.window);
end

%==========================================================================
%                  Generate Practice Schedule (unchanged)
%==========================================================================
function [sequence_1_back_practice] = gen_1_back_practice(p)

    A_files = dir(fullfile(p.stim_dir,'prac_*_A.png'));
    B_files = dir(fullfile(p.stim_dir,'prac_*_B.png'));
    
    A_names = string({A_files.name});
    B_names = string({B_files.name});
    
    seq = { ...
        A_names(1), "foil", "A", "none"; ...
        A_names(3), "same", "A", "none"; ...
        A_names(3), "same", "B", "j"; ...
        A_names(4), "foil", "A", "none"; ...
        A_names(2), "similar", "A", "none"; ...
        B_names(2), "similar", "B", "k" ...
    };
    
    col_names = {'stim_id','condition','identity','corr_resp','fix_duration'};
    variable_types = {'string','string','string','string','double'};
    
    sequence_1_back_practice = table('Size',[0,5], 'VariableTypes', variable_types, ...
        'VariableNames', col_names);
    
    for i = 1:size(seq,1)
        stim_id = seq{i,1};
        condition = seq{i,2};
        identity = seq{i,3};
        corr_resp = seq{i,4};
        fix_dur = 0.75 + (rand*0.5 - 0.25);
        
        new_row = table(stim_id, condition, identity, corr_resp, fix_dur, ...
            'VariableNames', col_names);
        sequence_1_back_practice = [sequence_1_back_practice; new_row];
    end
end
