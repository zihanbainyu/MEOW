
function D_run_2_back_practice(p)
    
    sequence_2_back_practice = gen_2_back_practice(p);
    nTrials = height(sequence_2_back_practice);
    
    % Keys
    start_key = KbName('f');
    escape_key = KbName(p.keys.quit);
    same_key = KbName({'1!','1'});      % top-row '1' and numpad '1'
    similar_key = KbName({'2@','2'});   % top-row '2' and numpad '2'
    
    Screen('TextSize', p.window, p.text_size);
    Screen('TextFont', p.window, 'Helvetica');
    
    % Setup KbQueue
    responseKeys = zeros(1, 256);
    responseKeys([same_key, similar_key, escape_key]) = 1;
    KbQueueCreate(p.keys.device, responseKeys);

    %% Initial fixation
    draw_fixation_target(p);
    Screen('Flip', p.window);
    WaitSecs(2);
    
    %% Trial loop with feedback
    KbQueueStart(p.keys.device);
    n_correct = 0;
    
    for i = 1:nTrials
        
        KbQueueFlush(p.keys.device);
        
        % Load image
        img_path = fullfile(p.stim_dir, sequence_2_back_practice.stim_id{i});
        if ~exist(img_path, 'file')
            error('Cannot find: %s', img_path);
        end
        img_data = imread(img_path);
        img_tex = Screen('MakeTexture', p.window, img_data);
        
        % Fixation
        draw_fixation_target(p);
        fix_onset = Screen('Flip', p.window);

        % Stimulus (dontclear = 1 keeps the image in the buffer so the
        % response highlight is drawn on top of it, without a redraw)
        Screen('DrawTexture', p.window, img_tex, [], [], 0);
        stim_onset = Screen('Flip', p.window, fix_onset + sequence_2_back_practice.fix_duration(i) - 0.5*p.ifi, 1);
        
        % Collect response
        responded = false;
        response_key = '';

        while GetSecs < stim_onset + p.timing.image_dur && ~responded
            [pressed, firstPress] = KbQueueCheck(p.keys.device);
            
            if pressed && ~responded
                key_code = find(firstPress, 1);
                
                if key_code == escape_key
                    error('USER_ABORT');
                elseif any(key_code == same_key)
                    response_key = '1';
                    responded = true;
                elseif any(key_code == similar_key)
                    response_key = '2';
                    responded = true;
                end
            end
        end
        
        % If no response, set to 'none'
        if ~responded
            response_key = 'none';
        end
        
        % Evaluate correctness
        correct_resp = sequence_2_back_practice.corr_resp{i};
        is_correct = strcmp(response_key, correct_resp);
        
        if is_correct
            n_correct = n_correct + 1;
        end
        
        % Highlight the image already on screen (no redraw): frame it green
        % if correct, red if incorrect. Incorrect also shows the response and
        % the correct answer below the image.
        [ih, iw, ~] = size(img_data);
        frame_rect = CenterRectOnPointd([0 0 iw ih], p.centerX, p.centerY) + [-10 -10 10 10];

        if is_correct
            Screen('FrameRect', p.window, [0 1 0], frame_rect, 8);   % green
            Screen('Flip', p.window);
            WaitSecs(0.75);   % quick confirmation, minimal disturbance
        else
            Screen('FrameRect', p.window, [1 0 0], frame_rect, 8);   % red
            fb = sprintf('You pressed: %s\n\nCorrect answer: %s', ...
                         resp_label(response_key), resp_label(correct_resp));
            DrawFormattedText(p.window, fb, 'center', frame_rect(4) + 50, p.colors.black);
            Screen('Flip', p.window);
            WaitSecs(2);      % time to read the correct answer
        end
        
        Screen('Close', img_tex);
    end
    
    KbQueueRelease(p.keys.device);
    
    %% Performance summary
    accuracy = (n_correct / nTrials) * 100;
    pass_thresh = 100;   % percent correct required to move on to the real task

    if accuracy >= pass_thresh
        summary_text = sprintf(...
            ['Accuracy: %.0f%%\n\n\n\n' ...
            'The real task provides no feedback.\n\n' ...
            'press any button to end.'], accuracy);
        summary_color = p.colors.black;
    else
        summary_text = sprintf(...
            ['Accuracy: %.0f%%\n\n\n\n' ...
            'The practice will be repeated once more.\n\n' ...
            'Press any button to continue. If anything is unclear, please ask the experimenter.'], accuracy);
        summary_color = p.colors.black;
    end

    DrawFormattedText(p.window, summary_text, 'center', 'center', summary_color);
    Screen('Flip', p.window);

    KbReleaseWait(p.keys.device);
    while true
        [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
        if keyIsDown
            if any(keyCode([same_key, similar_key]))   % any button advances
                if accuracy < pass_thresh
                    D_run_2_back_practice(p);   % repeat practice
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
function [sequence_2_back_practice] = gen_2_back_practice(p)
    
    A_files = dir(fullfile(p.stim_dir,'prac_*_A.png'));
    B_files = dir(fullfile(p.stim_dir,'prac_*_B.png'));

    A_names = string({A_files.name});
    B_names = string({B_files.name});
    
    seq = { ...
        A_names(5), "foil", "A", "none"; ... 
        A_names(6), "same", "A", "none"; ... 
        A_names(11), "similar", "A", "none"; ...
        A_names(6), "same", "A", "1";...
        B_names(11), "similar", "B", "2";   ... 
        A_names(12), "similar", "A", "none"; ...
        A_names(13), "same", "A", "none";...
        B_names(12), "similar", "B", "2"; ...
        A_names(13), "same", "A", "1";...
        A_names(14), "similar", "A", "none";...
        A_names(15), "foil", "A", "none";...
        B_names(14), "similar", "B", "2";...
    };
    
    col_names = {'stim_id','condition','identity','corr_resp','fix_duration'};
    variable_types = {'string','string','string','string','double'}; 
    
    sequence_2_back_practice = table('Size',[0,5], 'VariableTypes', variable_types, ...
                              'VariableNames', col_names);
    
    for i = 1:size(seq,1)
        stim_id = seq{i,1};
        cond = seq{i,2};
        role = seq{i,3};
        correct_resp = seq{i,4};
        fix_dur = 0.75 + (rand*0.5 - 0.25);
        
        new_row = table(stim_id, cond, role, correct_resp, fix_dur, ...
                    'VariableNames', col_names);
        sequence_2_back_practice = [sequence_2_back_practice; new_row];
    end
end

function label = resp_label(key)
% Human-readable label for a response code, for feedback text.
switch key
    case '1',  label = 'SAME';
    case '2',  label = 'SIMILAR';
    otherwise, label = 'nothing (NEW)';
end
end

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