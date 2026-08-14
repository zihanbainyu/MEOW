%==========================================================================
%   run_instructions.m
%==========================================================================

function run_instructions(p, names, advance_keys)

    % In-scanner: the subject advances each instruction with the button box
    % (button 1 or 2). Top-row and numpad names are both accepted. Pass the
    % optional advance_keys (e.g. KbName('f')) to override this for an
    % experimenter-advanced screen, such as the MST start page.
    if nargin < 3 || isempty(advance_keys)
        advance_keys = KbName({'1!','1','2@','2'});
    end
    escape_key   = KbName(p.keys.quit);

    if ischar(names) || isstring(names)
        names = cellstr(names);
    end

    exts = {'.jpg', '.jpeg', '.png'};   % prefer .jpg, fall back to .png
    for k = 1:numel(names)
        base  = char(names{k});
        fname = '';
        for e = 1:numel(exts)
            cand = fullfile(p.instr_dir, [base exts{e}]);
            if exist(cand, 'file'), fname = cand; break; end
        end
        if isempty(fname)
            error('run_instructions: instruction image not found: %s(.jpg/.jpeg/.png)', ...
                  fullfile(p.instr_dir, base));
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
                if any(kc(advance_keys))
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
