%==========================================================================
%   run_instructions.m
%==========================================================================
% Instruction screens drawn on the SAME grey background as the task.
%
% The seven text screens (ins_start, ins_1, ins_2, ins_6, ins_10,
% ins_prac_1back, ins_prac_2back) are defined as text in instruction_content()
% below -- edit the wording there, no image editing needed. The six object-
% example screens (ins_3/4/5/7/8/9) still load their PNGs from p.instr_dir
% because they contain object photographs.
%==========================================================================
function run_instructions(p, names)

    start_key  = KbName('f');
    escape_key = KbName(p.keys.quit);

    if ischar(names) || isstring(names)
        names = cellstr(names);
    end

    for k = 1:numel(names)
        name = char(names{k});
        Screen('TextFont', p.window, 'Helvetica');
        Screen('FillRect', p.window, p.colors.bgcolor);

        [is_text, title, body, footer] = instruction_content(name);
        tex = [];
        if is_text
            draw_text_slide(p, title, body, footer);
        else
            % object-example slide: still a PNG (contains object photos)
            fname = fullfile(p.instr_dir, [name '.png']);
            if ~exist(fname, 'file')
                error('run_instructions: instruction image not found: %s', fname);
            end
            img   = imread(fname);
            tex   = Screen('MakeTexture', p.window, img);
            scale = min(p.screenX / size(img,2), p.screenY / size(img,1));
            dest  = CenterRectOnPointd([0 0 size(img,2)*scale size(img,1)*scale], ...
                                       p.centerX, p.centerY);
            Screen('DrawTexture', p.window, tex, [], dest);
        end

        Screen('Flip', p.window);

        % 'f' advances, escape aborts
        KbReleaseWait(p.keys.device);
        while true
            [down, ~, kc] = KbCheck(p.keys.device);
            if down
                if kc(start_key)
                    break;
                elseif kc(escape_key)
                    sca;
                    error('Experiment terminated by user.');
                end
            end
            WaitSecs(0.001);
        end
        KbReleaseWait(p.keys.device);
        if ~isempty(tex), Screen('Close', tex); end
    end
end

%% ------------------------------------------------------------------------
function draw_text_slide(p, title, body, footer)
% Title (large, centred, top), body (wrapped, left, middle), footer (centred,
% bottom) -- all in black on the grey background already filled by the caller.
    if ~isempty(title)
        Screen('TextSize', p.window, round(p.text_size * 1.35));
        DrawFormattedText(p.window, title, 'center', round(0.15 * p.screenY), p.colors.black);
    end
    Screen('TextSize', p.window, p.text_size);
    if ~isempty(body)
        % wrapat = 60 chars wraps the prose paragraphs; intentional short lines
        % below are kept under this width so they are not re-wrapped.
        DrawFormattedText(p.window, body, round(0.14 * p.screenX), round(0.30 * p.screenY), ...
            p.colors.black, 60, [], [], 1.4);
    end
    if ~isempty(footer)
        DrawFormattedText(p.window, footer, 'center', round(0.84 * p.screenY), p.colors.black);
    end
    Screen('TextSize', p.window, p.text_size);
end

%% ------------------------------------------------------------------------
function [is_text, title, body, footer] = instruction_content(name)
% Returns the text of each instruction screen. Edit wording here. Emphasis is
% carried by CAPITALS (DrawFormattedText has no inline bold). Screens not
% listed are object-example PNGs (is_text = false).
title = ''; body = ''; footer = ''; is_text = true;
switch name
    case 'ins_start'
        title  = 'Welcome to our study!';
        footer = 'Press f to start with instructions.';

    case 'ins_1'
        body = sprintf([ ...
            'Inside the scanner:\n\n' ...
            'Part 1 -- You will view everyday objects and respond to each ' ...
            'one, in 4 blocks. In every block, you will perform a one-back ' ...
            'task, followed by a two-back task.\n\n' ...
            'Part 2 -- After Part 1, you will perform a memory test.\n\n\n' ...
            'Outside the scanner:\n\n' ...
            'Part 3 -- At the end, we will bring you out of the scanner and ' ...
            'you will perform another memory test.\n\n\n' ...
            'Please stay focused during each task. Short breaks will be ' ...
            'provided every ~6 mins.']);
        footer = 'Press f to continue.';

    case 'ins_2'
        title = 'Part 1 -- One-Back';
        body  = sprintf([ ...
            'In this task, objects appear one at a time. Your job is to ' ...
            'compare each object to the one RIGHT BEFORE IT, and decide:\n\n\n' ...
            'If it looks exactly the SAME as the one before it, press with ' ...
            'your INDEX finger.\n\n' ...
            'If it looks SIMILAR (slightly altered details, colours, ' ...
            'orientations, etc.) but not identical, press with your MIDDLE ' ...
            'finger.\n\n' ...
            'If it looks completely DIFFERENT, do not press anything.']);
        footer = 'Press f to view examples.';

    case 'ins_6'
        title = 'Part 1 -- Two-Back';
        body  = sprintf([ ...
            'After the one-back task, your job is now to compare each image ' ...
            'to the one TWO IMAGES BACK, using the same keys:\n\n\n' ...
            'SAME: press 1 with your INDEX finger.\n\n' ...
            'SIMILAR: press 2 with your MIDDLE finger.\n\n' ...
            'DIFFERENT: do not press anything.']);
        footer = 'Press f to view examples.';

    case 'ins_10'
        title = 'Important note';
        body  = sprintf([ ...
            'The comparison moves forward one image at a time. For example:\n\n' ...
            'When image 3 appears, compare it to image 1.\n' ...
            'When image 4 appears, compare it to image 2.\n' ...
            'When image 5 appears, compare it to image 3.\n' ...
            '...\n\n' ...
            'Every image may become a target two steps later, so hold each ' ...
            'one in mind. Never drop the image in between -- it may be the ' ...
            'next thing you compare to.']);
        footer = 'Press f to continue.';

    case 'ins_prac_1back'
        title = 'Practice: One-Back';
        body  = sprintf([ ...
            'You will complete a few practice trials. After each image, you ' ...
            'will see whether your answer was correct.\n\n' ...
            'If it was correct, a green rectangle will be highlighted on the ' ...
            'screen. If it was incorrect, a red rectangle will be ' ...
            'highlighted, and you will be shown the correct answer.\n\n' ...
            'This feedback is provided ONLY during practice.']);
        footer = 'Press f to start.';

    case 'ins_prac_2back'
        title = 'Practice: Two-Back';
        body  = sprintf([ ...
            'You will complete a few practice trials. After each image, you ' ...
            'will see whether your answer was correct.\n\n' ...
            'If it was correct, a green rectangle will be highlighted on the ' ...
            'screen. If it was incorrect, a red rectangle will be ' ...
            'highlighted, and you will be shown the correct answer.\n\n' ...
            'This feedback is provided ONLY during practice.']);
        footer = 'Press f to start.';

    otherwise
        % ins_3, ins_4, ins_5, ins_7, ins_8, ins_9 -> object-example PNGs
        is_text = false;
end
end
