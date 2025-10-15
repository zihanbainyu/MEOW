%==========================================================================
%                      Instruction Function
%==========================================================================
% author: Zihan Bai
%==========================================================================

function instructions(p, instruction_type)

    % Set common text parameters
    Screen('TextSize', p.window, 32);
    %Screen('TextFont', p.window, 'Open Sans');
    line_spacing = 1;

    % Define key names for waiting for responses
    KbName('UnifyKeyNames');
    space_key = KbName('g');
    escape_key = KbName(p.keys.quit);

    % Use a 'switch' statement to select the correct set of instructions
    switch instruction_type

        case 'welcome'
        DrawFormattedText(p.window, ['Welcome to our study!\n\n' ...
            '  \n\n' ...
            'Today, you will view items and make responses across two phases,\n\n' ...
            'We will explain each phase in detail before it begins.\n\n' ...
            '  \n\n' ...
            'When you are ready, please let the experimenter know.\n'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
        Screen('Flip', p.window);
        waitForKeyPress(space_key, escape_key);
        
        %==================================================================
        % case 1: encoding
        %==================================================================
        case 'encoding'

            DrawFormattedText(p.window, ['Phase I: 1-Back Repeat-Similar Task\n\n' ...
            'In this task, you will see a continuous stream of items\n\n' ...
            'For each item, you need to compare it against the one shown immediately before it\n\n'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);

            DrawFormattedText(p.window, ['You have three possible responses:\n\n' ...
                '  \n\n' ...
                'If the object is exactly the same as the previous one, press j (repeat)\n\n' ...
                'If the object is similar but not identical, press k (similar)\n\n' ...
                'If the object is new, do not press any key\n\n'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);
            
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
        
            DrawFormattedText(p.window, 'Before you begin, I am going to show you a brief demonstration to give you a feel for the task', 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);
            
            DrawFormattedText(p.window, 'Right before each item appears, you will see a fixation cross at the center of the screen:', ...
                'center', p.yCenter - 250, p.colors.black, [], [], [], line_spacing);
            Screen('TextSize', p.window, 100);
            DrawFormattedText(p.window, '+', 'center', 'center', p.colors.black); 
            Screen('TextSize', p.window, 32);
            DrawFormattedText(p.window, 'Please look at it, not anywhere else', ...
                'center', p.yCenter + 250, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);

            DrawFormattedText(p.window, 'Then, an item may show up:', 'center', p.yCenter - 250, p.colors.black, [], [], [], line_spacing);
            Screen('DrawTexture', p.window, tex_repeat, [], [], 0);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);
            
            DrawFormattedText(p.window, 'If the next item looks like this:', ...
                'center', p.yCenter - 250, p.colors.black, [], [], [], line_spacing);
            Screen('DrawTexture', p.window, tex_repeat, [], [], 0);
            DrawFormattedText(p.window, ['You would press j (repeat)'], ...
                'center', p.yCenter + 250, p.colors.black);
            DrawFormattedText(p.window, 'Because it looks exactly the same as the previous one', ...
                'center', p.yCenter + 300, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);


            DrawFormattedText(p.window, 'However, if the next item looks like this:', ...
                'center', p.yCenter - 250, p.colors.black, [], [], [], line_spacing);
            Screen('DrawTexture', p.window, tex_similar, [], [], 0);
            DrawFormattedText(p.window, ['You would press k (similar)'], ...
                'center', p.yCenter + 250, p.colors.black);
            DrawFormattedText(p.window, 'Because it looks somewhat similar, but not exactly the same', ...
                'center', p.yCenter + 300, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);

            DrawFormattedText(p.window, 'Or, if the next item looks like this:', ...
                'center', p.yCenter - 250, p.colors.black, [], [], [], line_spacing);
            Screen('DrawTexture', p.window, tex_new, [], [], 0);
            DrawFormattedText(p.window, 'You would not press any key.', ...
                'center', p.yCenter + 250, p.colors.black);
            DrawFormattedText(p.window, 'Because this item is completely new and unrelated to the previous one', ...
                'center', p.yCenter + 300, p.colors.black);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);
            
            DrawFormattedText(p.window, ['Please try your best to respond as quickly and accurately as possible\n\n' ...
            'Keep your head and body still during this part\n\n' ...
            'You will have several chances to take breaks\n\n'], ...
            'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);

            % Final confirmation
            DrawFormattedText(p.window, ['If you have any questions, please ask the experimenter now\n\n' ...
                'Otherwise, we will begin with a short practice round\n\n' ...
                'When you are ready, press g to begin'], ...
                'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            waitForKeyPress(space_key, escape_key);         
  
            
        case 'goodbye'
        %==================================================================
        % Case 3: the final screen of the experiment
        %==================================================================
            end_text = 'You are all done.\n\nThank you very much for your time and effort!\n\nPlease find the experimenter';
            DrawFormattedText(p.window, end_text, 'center', 'center', p.colors.black, [], [], [], line_spacing);
            Screen('Flip', p.window);
            
            % Wait for 5 seconds before allowing an exit, to prevent accidental closure
            WaitSecs(5);
            KbWait([], 2);

        otherwise
            error('Invalid instruction_type provided to instructions function: %s', instruction_type);
    
    end % switch statement

end % function


%==========================================================================
%                  HELPER FUNCTION
%==========================================================================
function waitForKeyPress(space_key, escape_key)
    % Flush any previous keypresses
    KbReleaseWait; % Wait until all keys are released
    KbQueueCreate();
    KbQueueStart();
    while true
        [pressed, firstPress] = KbQueueCheck();
        if pressed
            if firstPress(space_key)
                break;
            elseif firstPress(escape_key)
                sca; % exit screen
                error('Experiment terminated by user.');
            end
        end
    end
    KbQueueStop();
    KbQueueRelease();
    KbReleaseWait; % Wait again for key release before next screen
end
