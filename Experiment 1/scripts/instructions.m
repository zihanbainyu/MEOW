%==========================================================================
%                      Instruction Function
%==========================================================================
% author: Zihan Bai
%==========================================================================

function instructions(p, instruction_type)
    line_spacing = 1;

    % Define key names for waiting for responses
    KbName('UnifyKeyNames');
    space_key = KbName('g');
    escape_key = KbName(p.keys.quit);

    % Use a 'switch' statement to select the correct set of instructions
    switch instruction_type

        case 'welcome'
        DrawFormattedText(p.window, ['Welcome to our experiment!\n\n' ...
            'Thank you for participating.\n\n' ...
            '\n\n' ...
            'This experiment has 4 blocks. Each block has three parts:\n\n' ...
            '1. A 1-Back memory task\n\n' ...
            '2. A short rest\n\n' ...
            '3. A 2-Back memory task\n\n' ...
            '\n\n' ...
            'You will receive specific instructions before each task begins.\n\n' ...
            '  \n\n' ...
            'When you are ready, please let the experimenter know.\n'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
        Screen('Flip', p.window);
        waitForKeyPress(p, space_key, escape_key);
        
        %==================================================================
        % case 1: encoding
        %==================================================================
        case 'encoding'
            % pre-load example images
            try
                img_repeat = imread(fullfile(p.stim_dir, 'example_targ.png'));
                img_similar = imread(fullfile(p.stim_dir, 'example_lure.png'));
                img_new = imread(fullfile(p.stim_dir, 'example_new.png'));
            catch
                error('could not find example images!');
            end
            
            tex_repeat = Screen('MakeTexture', p.window, img_repeat);
            tex_similar = Screen('MakeTexture', p.window, img_similar);
            tex_new = Screen('MakeTexture', p.window, img_new);

            % task overview
            DrawFormattedText(p.window, ['Part 1: 1-Back Memory Task\n\n' ...
                '\n\n' ...
                'In this task, you will view a series of images, one after another.\n\n' ...
                'Your task is to compare each image to the one that came right before it.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

            DrawFormattedText(p.window, ['You have three possible responses:\n\n' ...
                '\n\n' ...
                'Press  j  (repeat) if the image is exactly the same as the one just before it.\n\n' ...
                'Press  k  (similar) if the image is similar but not identical as the one just before it.\n\n' ...
                'Do not press any key if the image is completely new.\n\n' ...
                '\n\n' ...
                'You are always comparing the current image to the one just before it.'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

            % demonstration
            DrawFormattedText(p.window, 'Before you begin, we will show you a brief demonstration to give you a feel for the task.', 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

            % Example 1: REPEAT
            DrawFormattedText(p.window, 'Previous Image', p.xCenter - 450, p.yCenter - 425, p.colors.black);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 300 300], p.xCenter - 350, p.yCenter - 225));
            DrawFormattedText(p.window, 'Current Image', p.xCenter + 250, p.yCenter - 425, p.colors.black);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 300 300], p.xCenter + 350, p.yCenter - 225));
            DrawFormattedText(p.window, 'You would press  j  because it is exactly the same.', 'center', p.yCenter + 50, p.colors.black);
        
            % Example 2: LURE
            DrawFormattedText(p.window, 'Previous Image', p.xCenter - 450, p.yCenter + 75, p.colors.black);
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 300 300], p.xCenter - 350, p.yCenter + 275));
            DrawFormattedText(p.window, 'Current Image', p.xCenter + 250, p.yCenter + 75, p.colors.black);
            Screen('DrawTexture', p.window, tex_similar, [], CenterRectOnPointd([0 0 300 300], p.xCenter + 350, p.yCenter + 275));
            DrawFormattedText(p.window, 'You would press  k  because it is similar but not identical.', 'center', p.yCenter + 525, p.colors.black);
            
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

            % fixation crosss
            DrawFormattedText(p.window, ['Before each image, you will see a fixation cross.\n\n' ...
                'Please keep your eyes on it.'], 'center', p.yCenter - 250, p.colors.black);
            xCoords = [-p.fix_cross_size, p.fix_cross_size, 0, 0];
            yCoords = [0, 0, -p.fix_cross_size, p.fix_cross_size];
            allCoords = [xCoords; yCoords];
            Screen('DrawLines', p.window, allCoords, p.fix_cross_width, p.colors.black, [p.xCenter p.yCenter]);

            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);
                    
         
            % quesitons
            DrawFormattedText(p.window, ['Please respond as quickly and accurately as you can.\n\n' ...
                'It is also important that you keep your head and body still.\n\n' ...
                'You will have several chances to take breaks.\n\n' ...
                '\n\n' ...
                'If you have any questions, please ask the experimenter now.\n\n' ...
                'Otherwise, we will begin with a short practice round.\n\n' ...
                'When you are ready, press g to begin.'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

            if exist('tex_repeat', 'var'),  Screen('Close', tex_repeat);  end
            if exist('tex_similar', 'var'), Screen('Close', tex_similar); end
            if exist('tex_new', 'var'),     Screen('Close', tex_new);     end

        %==================================================================
        % case 2: memory
        %==================================================================
        case 'memory'
            % pre-load example images
            try
                img_repeat = imread(fullfile(p.stim_dir, 'example_targ.png'));
                img_similar = imread(fullfile(p.stim_dir, 'example_lure.png'));
                img_new = imread(fullfile(p.stim_dir, 'example_new.png'));
            catch
                error('could not find example images!');
            end
            
            tex_repeat = Screen('MakeTexture', p.window, img_repeat);
            tex_similar = Screen('MakeTexture', p.window, img_similar);
            tex_new = Screen('MakeTexture', p.window, img_new);

            
            % task overview
            DrawFormattedText(p.window, ['Part 2: 2-Back Memory Task\n\n' ...
                '\n\n' ...
                'This task is similar, but with one important change.\n\n' ...
                'Your task is now to compare each image to the one from two trials ago.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

            DrawFormattedText(p.window, ['Therefore, the response rules now are:\n\n' ...
                '\n\n' ...
                'Press  j  (repeat) if the image is exactly the same as the one from two trials ago.\n\n' ...
                'Press  k  (similar) if the image is similar but not identical as the one from two trials ago.\n\n' ...
                'Do not press any key if the image is completely new.\n\n' ...
                '\n\n' ...
                'Important: If an image is repeats or resembles from only one trial ago, do nothing!\n\n' ...
                'You are only comparing the current image to the one from two trials ago.\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

            % demonstration
            DrawFormattedText(p.window, 'Before you begin, we will show you a brief demonstration to give you a feel for the task.', 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);


            % ========================================================================
            %                   Example 1: REPEAT (2-back)
            % ========================================================================
            
            % --- Labels ---
            DrawFormattedText(p.window, 'Two Trials Ago', p.xCenter - 800, p.yCenter - 400, p.colors.black);
            DrawFormattedText(p.window, 'Previous Image', p.xCenter - 100, p.yCenter - 400, p.colors.black);
            DrawFormattedText(p.window, 'Current Image', p.xCenter + 600, p.yCenter - 400, p.colors.black);
            
            % --- Images ---
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 300 300], p.xCenter - 700, p.yCenter - 200));  % two trials ago
            Screen('DrawTexture', p.window, tex_new,  [], CenterRectOnPointd([0 0 300 300], p.xCenter,        p.yCenter - 200));  % previous
            Screen('DrawTexture', p.window, tex_repeat, [], CenterRectOnPointd([0 0 300 300], p.xCenter + 700, p.yCenter - 200));  % current
            
            % --- Instruction text ---
            DrawFormattedText(p.window, 'You would press  j  because it is exactly the same as the image two trials ago.', ...
                'center', p.yCenter + 25, p.colors.black);
            
            
            % ========================================================================
            %                   Example 2: LURE (2-back)
            % ========================================================================
            
            % --- Labels ---
            DrawFormattedText(p.window, 'Two Trials Ago', p.xCenter - 800, p.yCenter + 200, p.colors.black);
            DrawFormattedText(p.window, 'Previous Image', p.xCenter - 100, p.yCenter + 200, p.colors.black);
            DrawFormattedText(p.window, 'Current Image', p.xCenter + 600, p.yCenter + 200, p.colors.black);
            
            % --- Images ---
            Screen('DrawTexture', p.window, tex_repeat,  [], CenterRectOnPointd([0 0 300 300], p.xCenter - 700, p.yCenter + 400)); % two trials ago
            Screen('DrawTexture', p.window, tex_new,   [], CenterRectOnPointd([0 0 300 300], p.xCenter,        p.yCenter + 400)); % previous
            Screen('DrawTexture', p.window, tex_similar, [], CenterRectOnPointd([0 0 300 300], p.xCenter + 700, p.yCenter + 400)); % current
            
            % --- Instruction text ---
            DrawFormattedText(p.window, 'You would press  k  because it is similar to the image two trials ago, but not identical.', ...
                'center', p.yCenter + 625, p.colors.black);
            
            % --- Flip screen and wait ---
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);

        
            DrawFormattedText(p.window, ['Please respond as quickly and accurately as you can.\n\n' ...
                'It is also important that you keep your head and body still.\n\n' ...
                'You will have several chances to take breaks.\n\n' ...
                '\n\n' ...
                'If you have any questions, please ask the experimenter now.\n\n' ...
                'Otherwise, we will begin with a short practice round.\n\n' ...
                'When you are ready, press g to begin.'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(p, space_key, escape_key);
            
            if exist('tex_repeat', 'var'),  Screen('Close', tex_repeat);  end
            if exist('tex_similar', 'var'), Screen('Close', tex_similar); end
            if exist('tex_new', 'var'),     Screen('Close', tex_new);     end


        case 'goodbye'
        %==================================================================
        % Case 3: the final screen of the experiment
        %==================================================================
            end_text = 'You are all done.\n\nThank you very much for your time and effort!\n\nPlease find the experimenter';
            DrawFormattedText(p.window, end_text, 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            
            % Wait for 5 seconds before allowing an exit, to prevent accidental closure
            WaitSecs(5);
            KbWait(p.keys.device, 2);

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
    end
    
    KbReleaseWait(p.keys.device); % Wait for the key to be released before proceeding
end