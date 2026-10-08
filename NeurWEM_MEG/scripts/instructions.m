%==========================================================================
%                      Instructions
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%
% instructions(p, which)               subject advances with button 1 or 2
% instructions(p, which, advance_keys) e.g. KbName('f') for experimenter-paced
%
% which: 'oneback', 'oneback_pp', 'twoback', 'twoback_pp', 'final'
%        (practice, outside the MEG), 'calibration', 'mst' (in the MEG).
%
% A page line that is exactly "<<IMG:key>>" is replaced by a centered row of
% example images (the *_pp post-practice recaps show the same/similar objects
% the subject just saw in practice). Text and images are kept inside a
% centered safe box; check the fractions below on the MEG projection screen.
% Escape aborts.
%==========================================================================
function instructions(p, which, advance_keys)

    if nargin < 3 || isempty(advance_keys)
        advance_keys = KbName({'1!','1','2@','2'});   % subject's buttons
    end
    escape_key = KbName(p.keys.quit);

    [pages, imgsets] = instruction_pages(p, which);

    % ---- Safe-area / text parameters (tune here) -------------------------
    safe_w_frac     = 0.60;   % usable fraction of window WIDTH
    safe_h_frac     = 0.80;   % usable fraction of window HEIGHT
    instr_text_size = 24;     % instruction font (px); smaller => more text fits
    vSpacing        = 1.5;    % line-spacing multiplier
    % ---------------------------------------------------------------------

    Screen('TextSize', p.window, instr_text_size);
    Screen('TextFont', p.window, 'Helvetica');

    safe_px_w = safe_w_frac * p.screenX;
    safe_px_h = safe_h_frac * p.screenY;
    wrapat    = max(10, floor(safe_px_w / (0.5 * instr_text_size)));

    for k = 1:numel(pages)
        txt = char(pages{k});

        if contains(txt, '<<IMG:')
            draw_rich_page(p, txt, imgsets, safe_px_w, safe_px_h, instr_text_size, vSpacing, wrapat);
        else
            Screen('FillRect', p.window, p.colors.bgcolor);
            [~, ~, bbox] = DrawFormattedText(p.window, txt, 'center', 'center', ...
                p.colors.black, wrapat, [], [], vSpacing);
            box_w = bbox(3) - bbox(1);
            box_h = bbox(4) - bbox(2);
            if box_w > safe_px_w || box_h > safe_px_h
                warning(['Instruction page %d renders %.0fx%.0f px but the safe ' ...
                         'box is %.0fx%.0f px; it may be clipped. Shorten it or ' ...
                         'lower instr_text_size.'], k, box_w, box_h, safe_px_w, safe_px_h);
            end
        end
        Screen('Flip', p.window);

        KbReleaseWait(p.keys.device);
        while true
            [down, ~, kc] = KbCheck(p.keys.device);
            if down
                if any(kc(advance_keys))
                    break;
                elseif kc(escape_key)
                    error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
                end
            end
            WaitSecs(0.001);
        end
        KbReleaseWait(p.keys.device);
    end

    Screen('TextSize', p.window, p.text_size);
end

%==========================================================================
%                  PAGE CONTENT
%==========================================================================
function [pages, imgsets] = instruction_pages(p, which)
    imgsets = struct();

    switch lower(which)

        case 'oneback'
            pages = {
                sprintf(['Part 1: One-Back\n\n' ...
                'You will view a stream of objects, one at a time.\n\n' ...
                'Your job is to compare each object to the one RIGHT BEFORE it (one-back).\n\n' ...
                'Same object:        press the LEFT button\n\n' ...
                'Similar object:     press the RIGHT button\n\n' ...
                'Totally different:  do not press\n\n' ...
                'Pay attention to every object, because it may reappear in later tasks.\n\n\n\n' ...
                'Press any button to start a short practice.'])
            };

        case 'oneback_pp'
            imgsets = prac_example_images(p, '1back');
            pages = {
                sprintf(['As you may have noticed in the practice:\n\n\n\n' ...
                'This object looked exactly the SAME as the one before it,\n\n' ...
                '<<IMG:same>>\n\n' ...
                'so you would press the LEFT button.\n\n\n\n' ...
                'This object looked SIMILAR, but not identical, to the one before it,\n\n' ...
                '<<IMG:similar>>\n\n' ...
                'so you would press the RIGHT button.\n\n' ...
                'Press any button to continue.'])
            };

        case 'twoback'
            pages = {
                sprintf(['Part 1: Two-Back\n\n' ...
                'After that, you will again view a stream of objects.\n\n' ...
                'This time, compare each object to the one TWO IMAGES BACK.\n\n' ...
                'Same object:        press the LEFT button\n\n' ...
                'Similar object:     press the RIGHT button\n\n' ...
                'Totally different:  do not press\n\n\n\n' ...
                'Some objects from the one-back may return here.\n\n\n\n' ...
                'Press any button to continue.'])
                sprintf(['Important note:\n\n' ...
                'The comparison moves forward one object at a time. For example:\n\n' ...
                'when image 3 appears, compare it to image 1,\n\n' ...
                'when image 4 appears, compare it to image 2,\n\n' ...
                'when image 5 appears, compare it to image 3,\n\n' ...
                'and so on.\n\n' ...
                'Press any button to start a short practice.'])
            };

        case 'twoback_pp'
            imgsets = prac_example_images(p, '2back');
            pages = {
                sprintf(['As you may have noticed in the practice:\n\n\n\n' ...
                'This object looked exactly the SAME as the one two images back,\n\n' ...
                '<<IMG:same>>\n\n' ...
                'so you would press the LEFT button.\n\n\n\n' ...
                'This object looked SIMILAR, but not identical, to the one two images back,\n\n' ...
                '<<IMG:similar>>\n\n' ...
                'so you would press the RIGHT button.\n\n' ...
                'Press any button to continue.'])
            };

        case 'final'
            pages = {
                sprintf(['You have completed the instructions and practice.\n\n' ...
                'In the actual task, you will do a one-back, then a two-back,\n\n' ...
                'repeated over four separate runs.\n\n\n\n' ...
                'Press any button to continue.'])
            };

        case 'calibration'
            pages = {
                sprintf(['Before any further instructions, we need to calibrate the eye-tracker.\n\n\n\n' ...
                'Please lay comfortably and stay still.\n\n' ...
                'There will be dots in the next few screens.\n\n' ...
                'Follow the dot with your eyes.\n\n' ...
                'Fixate directly at its center until it moves.'])
            };

        case 'mst'
            pages = {
                sprintf(['Part 2: Memory Test\n\n' ...
                'Finally, you will view a series of objects.\n\n' ...
                'For each one, decide whether it is old, similar, or new\n\n' ...
                'compared with the objects you saw in Part 1.\n\n\n\n' ...
                'Old:      press the LEFT button\n\n' ...
                'Similar:  press the RIGHT button\n\n' ...
                'New:      do not press\n\n\n\n' ...
                'Please let the experimenter know when you are ready.'])
            };

        otherwise
            error('instructions: unknown page set "%s".', which);
    end
end

%==========================================================================
%                  PRACTICE EXAMPLE IMAGES
%==========================================================================
% The exact same/similar objects used in practice, so the recap shows what the
% subject just saw. Indices mirror gen_1_back_practice / gen_2_back_practice:
%   1-back: same = A(1) shown twice; similar = A(2)->B(2)
%   2-back: same = A(10) shown twice; similar = A(13)->B(13)
% Keep these in sync if the practice sequences change.
function imgs = prac_example_images(p, task)

    A = dir(fullfile(p.stim_dir, 'prac_*_A.png'));
    B = dir(fullfile(p.stim_dir, 'prac_*_B.png'));
    A = string(sort({A.name}));   % zero-padded names -> alphabetical == numeric
    B = string(sort({B.name}));

    switch lower(task)
        case '1back'
            i_same = 1;   i_sim = 2;
        case '2back'
            i_same = 10;  i_sim = 13;
    end

    need = max(i_same, i_sim);
    if numel(A) < need || numel(B) < i_sim
        error(['prac_example_images: not enough prac images in %s ' ...
               '(need A>=%d, B>=%d; found A=%d, B=%d).'], ...
               p.stim_dir, need, i_sim, numel(A), numel(B));
    end

    imgs.same    = {char(A(i_same)), char(A(i_same))};
    imgs.similar = {char(A(i_sim)),  char(B(i_sim))};
end

%==========================================================================
%                  TEXT + IMAGE PAGE LAYOUT
%==========================================================================
% Alternating text blocks and image rows, vertically centered inside the safe
% box. A line that is exactly "<<IMG:key>>" becomes a row of imgsets.(key).
function draw_rich_page(p, txt, imgsets, safe_w, safe_h, tsize, vspace, wrapat)
    Screen('TextSize', p.window, tsize);
    Screen('TextFont', p.window, 'Helvetica');

    lines = strsplit(txt, newline, 'CollapseDelimiters', false);
    types = {};  datas = {};  buf = {};
    for i = 1:numel(lines)
        tok = regexp(strtrim(lines{i}), '^<<IMG:(\w+)>>$', 'tokens', 'once');
        if ~isempty(tok)
            if ~isempty(buf)
                types{end+1} = 'text';  datas{end+1} = strjoin(buf, newline);  buf = {};
            end
            types{end+1} = 'img';  datas{end+1} = tok{1};
        else
            buf{end+1} = lines{i};
        end
    end
    if ~isempty(buf)
        types{end+1} = 'text';  datas{end+1} = strjoin(buf, newline);
    end

    img_gap  = round(0.04 * p.screenX);
    img_side = round(min(0.16 * p.screenY, (safe_w - img_gap) / 2));
    gap_v    = round(0.015 * p.screenY);
    n = numel(types);

    % measurement pass (drawn to the back buffer, cleared before the real draw)
    h = zeros(1, n);
    for i = 1:n
        if strcmp(types{i}, 'text')
            [~, ny] = DrawFormattedText(p.window, datas{i}, 'center', 0, ...
                p.colors.black, wrapat, [], [], vspace);
            h(i) = ny;
        else
            h(i) = img_side;
        end
    end
    total_h = sum(h) + gap_v * max(0, n - 1);
    if total_h > safe_h
        warning(['Recap page content is %.0f px tall but the safe box is %.0f px; ' ...
                 'it may be clipped. Shorten the text or lower img_side/instr_text_size.'], ...
                 total_h, safe_h);
    end

    Screen('FillRect', p.window, p.colors.bgcolor);
    y = round(p.screenY / 2 - total_h / 2);
    for i = 1:n
        if strcmp(types{i}, 'text')
            [~, ny] = DrawFormattedText(p.window, datas{i}, 'center', y, ...
                p.colors.black, wrapat, [], [], vspace);
            y = ny + gap_v;
        else
            key = datas{i};
            if ~isfield(imgsets, key) || isempty(imgsets.(key))
                error('instructions: no images supplied for <<IMG:%s>>.', key);
            end
            files = imgsets.(key);
            row_w = numel(files) * img_side + (numel(files) - 1) * img_gap;
            x0 = round(p.screenX / 2 - row_w / 2);
            for j = 1:numel(files)
                imgpath = fullfile(p.stim_dir, files{j});
                if ~exist(imgpath, 'file')
                    error('instructions: image not found: %s', imgpath);
                end
                tex = Screen('MakeTexture', p.window, imread(imgpath));
                Screen('DrawTexture', p.window, tex, [], [x0, y, x0 + img_side, y + img_side]);
                Screen('Close', tex);
                x0 = x0 + img_side + img_gap;
            end
            y = y + img_side + gap_v;
        end
    end
end
