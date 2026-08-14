%==========================================================================
%   run_text_instructions.m
%==========================================================================
% Text-based replacement for run_instructions: shows one or more full-screen
% instruction pages. Pass a cell array of char strings (one element = one page;
% see instruction_text.m for the content).
%
% Advancing: each page advances on the subject's button box (1 or 2) by
% default. Pass advance_keys (e.g. KbName('f')) to make a page experimenter-
% advanced, as the MST start screen does. Escape aborts the experiment.
%
% IMAGES: pass imgsets (a struct) to embed pictures. A line in a page that is
% exactly "<<IMG:key>>" is replaced by a centered row of the images in
% imgsets.(key), stacked in reading order with the surrounding text. Used by the
% post-practice recap (imgsets from prac_example_images). Pages with no markers
% are drawn as plain centered text.
%
% SCANNER SAFE AREA
% -----------------
% On the projector only the CENTRAL region of the image is visible: the bore
% clips left/right and corners, the eye tracker blocks a bottom strip. From the
% projector-screen measurements (eye-to-screen 83.5 cm) the screen is wide
% (~60 x 36.2 cm) but the visible region is roughly a 36 cm SQUARE, so usable
% HEIGHT ~ full screen height but usable WIDTH ~ 0.6 of it. Text/images are kept
% inside a centered safe box (tighter horizontally) and long lines wrap to it.
% Fractions are conservative and assume the image fills the screen -- VERIFY in
% the scanner and shrink them if any edge is cut off.
%==========================================================================
function run_text_instructions(p, pages, advance_keys, imgsets)

    if nargin < 3 || isempty(advance_keys)
        advance_keys = KbName({'1!','1','2@','2'});   % subject's button box
    end
    if nargin < 4 || isempty(imgsets)
        imgsets = struct();
    end
    escape_key = KbName(p.keys.quit);

    if ischar(pages) || isstring(pages)
        pages = cellstr(pages);
    end

    % ---- Safe-area / text parameters (tune here) -------------------------
    safe_w_frac     = 0.60;   % usable fraction of window WIDTH  (visible region ~square on a wide screen)
    safe_h_frac     = 0.80;   % usable fraction of window HEIGHT (small bottom loss to eye tracker)
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
            % Mixed text + image page (post-practice recap).
            draw_rich_page(p, txt, imgsets, safe_px_w, safe_px_h, instr_text_size, vSpacing, wrapat);
        else
            % Plain centered text page.
            Screen('FillRect', p.window, p.colors.bgcolor);
            [~, ~, bbox] = DrawFormattedText(p.window, txt, 'center', 'center', ...
                p.colors.black, wrapat, [], [], vSpacing);
            box_w = bbox(3) - bbox(1);
            box_h = bbox(4) - bbox(2);
            if box_w > safe_px_w || box_h > safe_px_h
                warning(['Instruction page %d renders %.0fx%.0f px but the scanner ' ...
                         'safe box is %.0fx%.0f px; it may be clipped. Shorten it or ' ...
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
end

% -------------------------------------------------------------------------
function draw_rich_page(p, txt, imgsets, safe_w, safe_h, tsize, vspace, wrapat)
% Lay out a page of alternating text blocks and image rows, vertically centered
% inside the safe box. A line that is exactly "<<IMG:key>>" becomes a row of the
% images in imgsets.(key); all other lines are text.
    Screen('TextSize', p.window, tsize);
    Screen('TextFont', p.window, 'Helvetica');

    % --- parse into ordered text / image elements (empty lines preserved) ---
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

    % --- sizes ---
    img_gap  = round(0.04 * p.screenX);
    img_side = round(min(0.16 * p.screenY, (safe_w - img_gap) / 2));   % square images
    gap_v    = round(0.015 * p.screenY);
    n = numel(types);

    % --- measurement pass: get each element's height (drawn to the back buffer,
    %     then cleared before the real draw, so nothing is shown) ---
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

    % --- draw pass: centered vertically ---
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
                error('run_text_instructions: no images supplied for <<IMG:%s>>.', key);
            end
            files = imgsets.(key);
            row_w = numel(files) * img_side + (numel(files) - 1) * img_gap;
            x0 = round(p.screenX / 2 - row_w / 2);
            for j = 1:numel(files)
                imgpath = fullfile(p.stim_dir, files{j});
                if ~exist(imgpath, 'file')
                    error('run_text_instructions: image not found: %s', imgpath);
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
