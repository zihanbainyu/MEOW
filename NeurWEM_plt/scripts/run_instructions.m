%==========================================================================
%   run_instructions.m  -- show full-screen instruction images in the task
%==========================================================================
% Presents the given instruction images (in p.instr_dir), scaled to fit the
% screen, one per keypress. 'f' advances, escape quits.
%
%   names : cell array of image basenames WITHOUT extension, shown in order,
%           e.g. {'ins_start','ins_1','ins_prac_1back'}. Each <name>.png must
%           exist in p.instr_dir.
%
% Same fit-to-screen logic as preview_ins.m. p.instr_dir is set in main.m.
%==========================================================================
function run_instructions(p, names)

    start_key  = KbName('f');
    escape_key = KbName(p.keys.quit);

    if ischar(names) || isstring(names)
        names = cellstr(names);
    end

    for k = 1:numel(names)
        fname = fullfile(p.instr_dir, [char(names{k}) '.png']);
        if ~exist(fname, 'file')
            error('run_instructions: instruction image not found: %s', fname);
        end

        img = imread(fname);
        tex = Screen('MakeTexture', p.window, img);

        scale = min(p.screenX / size(img, 2), p.screenY / size(img, 1));
        dest  = CenterRectOnPointd([0 0 size(img,2)*scale size(img,1)*scale], ...
                                   p.centerX, p.centerY);

        Screen('FillRect', p.window, p.colors.bgcolor);
        Screen('DrawTexture', p.window, tex, [], dest);
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
        Screen('Close', tex);
    end
end
