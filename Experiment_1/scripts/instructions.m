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
            DrawFormattedText(p.window, ['Welcome to our study!\n\n' ...
                'Thank you for participating.\n\n',], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            DrawFormattedText(p.window, [
                'This experiment has two main phases.\n\n' ...
                '         Phase 1: Four Blocks of Memory Tasks\n\n' ...
                'Phase 2: One Final Task\n\n\n\n' ...
                'We will begin with Phase 1.'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            DrawFormattedText(p.window, [
                'Phase 1\n\n' ...
                'Each block has two parts:\n\n' ...
                '1-Back Task and 2-Back Task\n\n\n\n' ...
                'You may take a rest after each task.\n\n' ...
                'We will record your eye movements during these tasks.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);


        case 'calibration'
            DrawFormattedText(p.window, ['Before any further instructions, we need to calibrate the eye-tracker.\n\n' ...
                '\n\n', ...
                'Please sit comfortably and rest your head on the chinrest.\n\n' ...
                'There will be dots in the next few screens.\n\n' ...
                'Follow the dot naturally with your eyes. Look directly at its center until it moves.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

        case '1_back'
            % pre-load example images
            try
                img_repeat = imread(fullfile(p.stim_dir, 'instr_A.png'));
                img_similar = imread(fullfile(p.stim_dir, 'instr_B.png'));
            catch
                error('could not find example images!');
            end
            
            tex_repeat = Screen('MakeTexture', p.window, img_repeat);
            tex_similar = Screen('MakeTexture', p.window, img_similar);

            % Task Overview Screen 1
            DrawFormattedText(p.window, ['1-Back Task\n\n' ...
                '\n\n' ...
                'In this task, you will view a series of images, one after another.\n\n' ...
                'Your task is to compare each image to the one from ONE TRIAL AGO.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % Response Rules Screen 2
            DrawFormattedText(p.window, ['You have three possible responses:\n\n' ...
                '\n\n' ...
                'Press  ' p.keys.same '  (SAME)  with your index finger if the image is exactly the same as the one from one trial ago.\n\n' ...
                'Press  ' p.keys.diff '  (SIMILAR) with your middle finger if the image is similar but not identical to the one from one trial ago.\n\n' ...
                'Do not press any key if the image is completely new (NEW).\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % Demonstration Screen 3 (Setup)
            DrawFormattedText(p.window, 'Before you begin, here is a demonstration.', 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % --- Define Layout Constants (Using values from 2-Back for consistency) ---
            IMAGE_SIZE_DISPLAY = 300;  % 300x300 pixels
            IMAGE_SPACING_HALF = 325;  % Half the distance between image centers (350/2 is too small, use 350)
            TEXT_OFFSET_Y      = 175;  % Vertical offset for labels above images
            TEXT_OFFSET_CAPTION = 200; % Vertical offset for captions below images
            ROW_HEIGHT         = 450;  % Vertical spacing between the center of Row 1 and Row 2
            
            % --- Horizontal Offsets (centered on p.xCenter) ---
            X_LEFT_ENC   = p.xCenter - IMAGE_SPACING_HALF; % 350 pixels left of center
            X_RIGHT_ENC  = p.xCenter + IMAGE_SPACING_HALF; % 350 pixels right of center

            % --- Row 1 Vertical Positions (SAME/REPEAT Example) ---
            Y_ROW1_CENTER = p.yCenter - 275; 
            Y_ROW1_LABEL  = Y_ROW1_CENTER - TEXT_OFFSET_Y;
            Y_ROW1_CAPTION= Y_ROW1_CENTER + TEXT_OFFSET_CAPTION;
            
            % --- Row 2 Vertical Positions (SIMILAR/LURE Example) ---
            Y_ROW2_CENTER = Y_ROW1_CENTER + ROW_HEIGHT;
            Y_ROW2_LABEL  = Y_ROW2_CENTER - TEXT_OFFSET_Y;
            Y_ROW2_CAPTION= Y_ROW2_CENTER + TEXT_OFFSET_CAPTION;

            
            % ========================================================================
            % Example 1: SAME/REPEAT (Centered around X_LEFT_ENC and X_RIGHT_ENC)
            % ========================================================================
            % --- Labels ---
            % Position text labels above the image centers
            DrawFormattedText(p.window, 'Previous Image', X_LEFT_ENC - 50, Y_ROW1_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current Image', X_RIGHT_ENC - 50, Y_ROW1_LABEL, p.colors.black);
            
            % --- Images ---
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT_ENC, Y_ROW1_CENTER));
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT_ENC, Y_ROW1_CENTER));
            
            % --- Instruction text (Caption) ---
            DrawFormattedText(p.window, ['You would press  ' p.keys.same '  because it is exactly the SAME.'], 'center', Y_ROW1_CAPTION, p.colors.black);
        
            
            % ========================================================================
            % Example 2: SIMILAR/LURE (Centered around X_LEFT_ENC and X_RIGHT_ENC)
            % ========================================================================
            % --- Labels ---
            DrawFormattedText(p.window, 'Previous Image', X_LEFT_ENC - 50, Y_ROW2_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current Image', X_RIGHT_ENC - 50, Y_ROW2_LABEL, p.colors.black);
            
            % --- Images ---
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT_ENC, Y_ROW2_CENTER));
            Screen('DrawTexture', p.window, tex_similar, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT_ENC, Y_ROW2_CENTER));
            
            % --- Instruction text (Caption) ---
            DrawFormattedText(p.window, ['You would press  ' p.keys.diff '  because it is SIMILAR but not identical.'], 'center', Y_ROW2_CAPTION, p.colors.black);
            
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % Fixation Cross Screen 4
            DrawFormattedText(p.window, ['Before each image, you will see a fixation cross.\n\n' ...
                'Please keep your eyes on it.'], 'center', p.yCenter - 250, p.colors.black);
            xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
            yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
            allCoords = [xCoords; yCoords];
            Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);
                    
            % Final Instructions Screen 5
            DrawFormattedText(p.window, ['Please respond as quickly and accurately as you can.\n\n' ...
                'Please keep your eyes on the fixation cross (+) between images.\n\n' ...
                '\n\n' ...
                'If you have any questions, please ask the experimenter now.\n\n' ...
                'Otherwise, we will begin with a short practice round.\n\n' ...
                'When you are ready, press f to begin.'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);

            % Clean up textures
            Screen('Close', [tex_repeat, tex_similar]);


        case '2_back'
            % pre-load example images
            try
                img_repeat = imread(fullfile(p.stim_dir, 'instr_A.png'));
                img_similar = imread(fullfile(p.stim_dir, 'instr_B.png'));
                img_new = imread(fullfile(p.stim_dir, 'instr_A_foil.png'));
            catch
                error('could not find example images!');
            end
            
            tex_repeat = Screen('MakeTexture', p.window, img_repeat);
            tex_similar = Screen('MakeTexture', p.window, img_similar);
            tex_new = Screen('MakeTexture', p.window, img_new);
            
            % Task Overview Screen 1
            DrawFormattedText(p.window, ['2-Back Task\n\n' ...
                '\n\n' ...
                'This task is similar, but with one important change.\n\n' ...
                'Your task is now to compare each image to the one from TWO TRIALS AGO.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);
            
            % Response Rules Screen 2
            DrawFormattedText(p.window, ['Therefore, the response rules now are:\n\n' ...
                '\n\n' ...
                'Press  ' p.keys.same '  (SAME) with your index finger if the image is exactly the same as the one from two trials ago.\n\n' ...
                'Press  ' p.keys.diff '  (SIMILAR) with your middle finger if the image is similar but not identical to the one from two trials ago.\n\n' ...
                'Do not press any key if the image is completely new (NEW).\n\n' ...
                '\n\n' ...
                'Important: If an image repeats or resembles the one from only one trial ago, do nothing!\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);
            
            % Demonstration Screen 3 (Setup)
            DrawFormattedText(p.window, 'Before you begin, here is a brief demonstration.', 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);
            
            % --- Define Layout Constants ---
            IMAGE_SIZE_DISPLAY = 300; % Displayed size in pixels
            IMAGE_SPACING = 350;      % Horizontal spacing between image centers
            TEXT_OFFSET_Y = 175;      % Vertical offset for labels above images
            TEXT_OFFSET_CAPTION = 200; % Vertical offset for captions below images
            ROW_HEIGHT = 450;         % Vertical spacing between the center of Row 1 and Row 2
            
            % --- Horizontal Offsets (centered on p.xCenter) ---
            X_LEFT   = p.xCenter - IMAGE_SPACING;
            X_CENTER = p.xCenter;
            X_RIGHT  = p.xCenter + IMAGE_SPACING;
            
            % --- Row 1 Vertical Positions ---
            Y_ROW1_CENTER = p.yCenter - 275; % Adjusted to start higher up
            Y_ROW1_LABEL  = Y_ROW1_CENTER - TEXT_OFFSET_Y;
            Y_ROW1_CAPTION= Y_ROW1_CENTER + TEXT_OFFSET_CAPTION;
            
            % --- Row 2 Vertical Positions ---
            Y_ROW2_CENTER = Y_ROW1_CENTER + ROW_HEIGHT;
            Y_ROW2_LABEL  = Y_ROW2_CENTER - TEXT_OFFSET_Y;
            Y_ROW2_CAPTION= Y_ROW2_CENTER + TEXT_OFFSET_CAPTION;
            
            
            % ========================================================================
            %                   Example 1: SAME (2-back)
            % ========================================================================
            
            % --- Labels ---
            % X_LEFT - 100 to align the text leftwards so its center is over the image center
            DrawFormattedText(p.window, 'Two Trials Ago', X_LEFT - 100, Y_ROW1_LABEL, p.colors.black); 
            DrawFormattedText(p.window, 'Previous Image', X_CENTER - 50, Y_ROW1_LABEL, p.colors.black); 
            DrawFormattedText(p.window, 'Current Image', X_RIGHT - 50, Y_ROW1_LABEL, p.colors.black); 
            
            % --- Images ---
            % Note: CenterRectOnPointd only takes one X and one Y point
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_ROW1_CENTER));  % two trials ago
            Screen('DrawTexture', p.window, tex_new,  [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_CENTER, Y_ROW1_CENTER));  % previous
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_ROW1_CENTER));  % current
            
            % --- Instruction text (Caption) ---
            DrawFormattedText(p.window, ['You would press  ' p.keys.same '  because it is exactly the SAME as the image two trials ago.'], ...
                'center', Y_ROW1_CAPTION, p.colors.black);
            
            
            % ========================================================================
            %                   Example 2: SIMILAR (2-back)
            % ========================================================================
            
            % --- Labels ---
            % The error was here: too many X/Y arguments passed
            DrawFormattedText(p.window, 'Two Trials Ago', X_LEFT - 100, Y_ROW2_LABEL, p.colors.black); 
            DrawFormattedText(p.window, 'Previous Image', X_CENTER - 50, Y_ROW2_LABEL, p.colors.black);
            DrawFormattedText(p.window, 'Current Image', X_RIGHT - 50, Y_ROW2_LABEL, p.colors.black);
            
            % --- Images ---
            Screen('DrawTexture', p.window, tex_repeat,  [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_LEFT, Y_ROW2_CENTER)); % two trials ago
            Screen('DrawTexture', p.window, tex_new,   [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_CENTER, Y_ROW2_CENTER)); % previous
            Screen('DrawTexture', p.window, tex_similar, [], CenterRectOnPointd([0 0 IMAGE_SIZE_DISPLAY IMAGE_SIZE_DISPLAY], X_RIGHT, Y_ROW2_CENTER)); % current
            
            % --- Instruction text (Caption) ---
            DrawFormattedText(p.window, ['You would press  ' p.keys.diff '  because it is SIMILAR to the image two trials ago, but not identical.'], ...
                'center', Y_ROW2_CAPTION, p.colors.black);
            
            % --- Flip screen and wait ---
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);
        
            % Final Instructions Screen 4
            DrawFormattedText(p.window, ['Please respond as quickly and accurately as you can.\n\n' ...
                'Please keep your eyes on the fixation cross (+) between images.\n\n' ...
                '\n\n' ...
                'If you have any questions, please ask the experimenter now.\n\n' ...
                'Otherwise, we will begin with a short practice round.\n\n' ...
                'When you are ready, press f to begin.'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, start_key, escape_key);
            
            % Clean up textures
            Screen('Close', [tex_repeat, tex_similar, tex_new]);

        case 'goodbye'
        %==================================================================
        % Case 3: the final screen of the experiment
        %==================================================================
            end_text = 'You are all done.\n\nThank you very much for your time and effort!\n\nPlease find the experimenter.';
            DrawFormattedText(p.window, end_text, 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            
            % Wait for 5 seconds before allowing an exit, to prevent accidental closure
            WaitSecs(5);
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
                error('Experiment terminated by user.');
            end
        end
        WaitSecs(0.001); % Yield CPU
    end
    
    KbReleaseWait(p.keys.device); % Wait for the key to be released before proceeding
end