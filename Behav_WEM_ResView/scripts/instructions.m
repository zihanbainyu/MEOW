%==========================================================================
%                      Instruction Function
%==========================================================================
% author: Zihan Bai
%==========================================================================
function instructions(p, instruction_type)
    line_spacing = 1; % Set line spacing to 1.0

    % Set font settings for this function, ensuring consistency (assumes window is open)
    Screen('TextSize', p.window, p.text_size);
    Screen('TextFont', p.window, 'Helvetica');
    
    % Define key names (start_key replaces space_key for consistency)
    start_key = KbName('f');
    escape_key = KbName(p.keys.quit);

    % Use a 'switch' statement to select the correct set of instructions
    switch instruction_type
        
        case 'welcome'
            DrawFormattedText(p.window, ['Welcome, and thank you for participating.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            DrawFormattedText(p.window, [
                'This study has two parts.\n\n\n\n' ...
                'Part 1: You will view and respond to a series of objects in 4 blocks.\n\n' ...
                'Because your eye movements are recorded, please keep your gaze\n\n' ...
                'on the central dot throughout.\n\n\n\n' ...
                'Part 2: At the end of the experiment, you will indicate which objects you remember.\n\n' ...
                'Your eyes are not tracked in this part, so you may look freely.\n\n\n\n' ...
                'Breaks are provided between blocks.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);


        case 'calibration'
            DrawFormattedText(p.window, ['First, the eye tracker will be calibrated.\n\n\n\n' ...
                'Please rest your head on the chinrest and sit comfortably.\n\n' ...
                'A dot will appear and move across the screen.\n\n' ...
                'Please follow it with your eyes, looking directly at its center.\n\n\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

        case '1_back'
            % pre-load example images
            try
                img_repeat = imread(fullfile(p.stim_dir, 'instr_A.png'));
                img_similar = imread(fullfile(p.stim_dir, 'instr_B.png'));
                img_new = imread(fullfile(p.stim_dir, 'instr_B_foil.png'));
            catch
                error('could not find example images!');
            end
            
            tex_repeat = Screen('MakeTexture', p.window, img_repeat);
            tex_similar = Screen('MakeTexture', p.window, img_similar);
            tex_new = Screen('MakeTexture', p.window, img_new);

            % --- Layout constants ---
            IMAGE_SIZE_DISPLAY  = 250;  % Size in pixels
            IMAGE_SPACING_HALF  = 350;  % Distance from center to image
            TEXT_OFFSET_Y       = 225;  % Label height above center
            TEXT_OFFSET_CAPTION = 350;  % Caption height below center
            X_LEFT  = p.xCenter - IMAGE_SPACING_HALF; % "One trial ago" position
            X_RIGHT = p.xCenter + IMAGE_SPACING_HALF; % "Current image" position
            Y_CENTER  = p.yCenter;
            Y_LABEL   = Y_CENTER - TEXT_OFFSET_Y;
            Y_CAPTION = Y_CENTER - TEXT_OFFSET_CAPTION;
            Y_CAPTION_2 = Y_CENTER + TEXT_OFFSET_Y;

            % Setup screen: task + fixation in one
            DrawFormattedText(p.window, ['1-Back Task\n\n\n' ...
                'Please compare each image to the one shown ONE TRIAL AGO,\n\n' ...
                'keeping your gaze on the central dot.\n\n\n' ...
                'Press f to view three examples.'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % ---- Example 1: SAME ----
            DrawFormattedText(p.window, '1 trial ago', X_LEFT - 70, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current image', X_RIGHT - 70, Y_LABEL, p.colors.black);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_CENTER));
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_CENTER));
            DrawFormattedText(p.window, 'SAME', 'center', Y_CAPTION, p.colors.black);
            DrawFormattedText(p.window, ['Identical to the previous image. Press ' p.keys.same '.'], 'center', Y_CAPTION_2, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % ---- Example 2: SIMILAR ----
            DrawFormattedText(p.window, '1 trial ago', X_LEFT - 70, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current image', X_RIGHT - 70, Y_LABEL, p.colors.black);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_CENTER));
            Screen('DrawTexture', p.window, tex_similar, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_CENTER));
            DrawFormattedText(p.window, 'SIMILAR', 'center', Y_CAPTION, p.colors.black);
            DrawFormattedText(p.window, ['Similar but not identical. Press ' p.keys.diff '.'], 'center', Y_CAPTION_2, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % ---- Example 3: NEW ----
            DrawFormattedText(p.window, '1 trial ago', X_LEFT - 70, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current image', X_RIGHT - 70, Y_LABEL, p.colors.black);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_CENTER));
            Screen('DrawTexture', p.window, tex_new, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_CENTER));
            DrawFormattedText(p.window, 'NEW', 'center', Y_CAPTION, p.colors.black);
            DrawFormattedText(p.window, 'Completely different. Press nothing.', 'center', Y_CAPTION_2, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % (practice round follows immediately and shows its own "press f to begin" screen)

            % Clean up textures
            Screen('Close', [tex_repeat, tex_similar, tex_new]);


        case '2_back'
            % pre-load example images
            try
                img_repeat = imread(fullfile(p.stim_dir, 'instr_A.png'));
                img_similar = imread(fullfile(p.stim_dir, 'instr_B.png'));
                img_foil = imread(fullfile(p.stim_dir, 'instr_foil_A.png'));
                img_new = imread(fullfile(p.stim_dir, 'instr_foil_1.png'));
            catch
                error('could not find example images!');
            end
            
            tex_repeat = Screen('MakeTexture', p.window, img_repeat);
            tex_similar = Screen('MakeTexture', p.window, img_similar);
            tex_new = Screen('MakeTexture', p.window, img_new);
            tex_foil = Screen('MakeTexture', p.window, img_foil);
            
            % --- Layout constants ---
            IMAGE_SIZE_DISPLAY  = 200;  % Slightly smaller to fit 3 in a row
            IMAGE_SPACING       = 300;  % Distance between centers
            TEXT_OFFSET_Y       = 225;  % Label height above center
            TEXT_OFFSET_CAPTION = 350;  % Caption height below center
            X_MID   = p.xCenter;                  % 1 trial ago (center)
            X_LEFT  = p.xCenter - IMAGE_SPACING;  % 2 trials ago (left)
            X_RIGHT = p.xCenter + IMAGE_SPACING;  % Current image (right)
            Y_CENTER  = p.yCenter;
            Y_LABEL   = Y_CENTER - TEXT_OFFSET_Y;
            Y_CAPTION = Y_CENTER - TEXT_OFFSET_CAPTION;
            Y_CAPTION_2 = Y_CENTER + TEXT_OFFSET_Y;

            % Setup screen: the one change + same keys/fixation
            DrawFormattedText(p.window, ['2-Back Task\n\n\n' ...
                'Please compare each image to the one shown TWO TRIALS AGO.\n\n' ...
                'The keys and fixation remain the same.\n\n\n' ...
                'Press f to view three examples.'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % ---- Example 1: SAME ----
            DrawFormattedText(p.window, '2 trials ago', X_LEFT - 60, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, '1 trial ago', X_MID - 60, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current image', X_RIGHT - 40, Y_LABEL, [0 130 0]);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_CENTER));
            Screen('DrawTexture', p.window, tex_foil, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_MID, Y_CENTER));
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_CENTER));
            DrawFormattedText(p.window, 'SAME', 'center', Y_CAPTION, p.colors.black);
            DrawFormattedText(p.window, ['Identical to the image two trials ago. Press ' p.keys.same '.'], 'center', Y_CAPTION_2, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % ---- Example 2: SIMILAR ----
            DrawFormattedText(p.window, '2 trials ago', X_LEFT - 60, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, '1 trial ago', X_MID - 60, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current image', X_RIGHT - 40, Y_LABEL, [0 130 0]);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_CENTER));
            Screen('DrawTexture', p.window, tex_foil, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_MID, Y_CENTER));
            Screen('DrawTexture', p.window, tex_similar, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_CENTER));
            DrawFormattedText(p.window, 'SIMILAR', 'center', Y_CAPTION, p.colors.black);
            DrawFormattedText(p.window, ['Similar but not identical. Press ' p.keys.diff '.'], 'center', Y_CAPTION_2, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % ---- Example 3: NEW ----
            DrawFormattedText(p.window, '2 trials ago', X_LEFT - 60, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, '1 trial ago', X_MID - 60, Y_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current image', X_RIGHT - 40, Y_LABEL, [0 130 0]);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_CENTER));
            Screen('DrawTexture', p.window, tex_foil, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_MID, Y_CENTER));
            Screen('DrawTexture', p.window, tex_new, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_CENTER));
            DrawFormattedText(p.window, 'NEW', 'center', Y_CAPTION, p.colors.black);
            DrawFormattedText(p.window, 'Completely different. Press nothing.', 'center', Y_CAPTION_2, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % (practice round follows immediately and shows its own "press f to begin" screen)

            % Clean up textures
            Screen('Close', [tex_repeat, tex_similar, tex_new, tex_foil]);

        case 'goodbye'
        %==================================================================
        % Case 3: the final screen of the experiment
        %==================================================================
            end_text = 'This concludes the study. Thank you very much for your participation.\n\nPlease notify the experimenter.';
            DrawFormattedText(p.window, end_text, 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            
            % Wait for 2 seconds before allowing an exit, to prevent accidental closure
            WaitSecs(2);
            KbWait(p.keys.device); % Wait for *any* key press on the device to signal final exit
            
        otherwise
            error('Invalid instruction_type provided to instructions function: %s', instruction_type);
    
    end % switch statement
end % function

%==========================================================================
%                  HELPER FUNCTION
%==========================================================================
function waitForKeyPress(p, continue_key, escape_key)
   
    KbReleaseWait(p.keys.device); % Wait until all keys on the specified keyboard are released
    
    while true
        [keyIsDown, ~, keyCode] = KbCheck(p.keys.device); % Check ONLY the specified keyboard
        
        if keyIsDown
            if keyCode(continue_key)
                break; % Exit the loop if the continue key is pressed
            elseif keyCode(escape_key)
                sca; % Clean up screen
                error('Experiment terminated by us=er.');
            end
        end
        WaitSecs(0.001); % Yield CPU
    end
    
    KbReleaseWait(p.keys.device); % Wait for the key to be released before proceeding
end