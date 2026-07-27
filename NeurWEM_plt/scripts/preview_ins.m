%==========================================================================
%   preview_ins.m  -- quick check that the instruction images fit the screen
%==========================================================================
% Opens a full-screen window on the task background grey and steps through
% every stimulus/ins_*.png (sorted by name), each scaled to fit (aspect
% ratio preserved, centred). Any key = next image, ESC = quit early.
% Nothing else in the experiment is touched.
%
% Run from anywhere:  >> preview_insfffff
%==========================================================================
function preview_ins()

    here = fileparts(mfilename('fullpath'));
    stim = fullfile(here, '..', 'stimulus');

    files = dir(fullfile(stim, 'ins_*.png'));
    if isempty(files)
        error('No ins_*.png found in %s', stim);
    end
    [~, ord] = sort({files.name});     % filename order = show order
    files = files(ord);

    try
        Screen('Preference', 'SkipSyncTests', 1);
        PsychDefaultSetup(2);
        bg = [124 124 124] / 255;              % same grey as the task

        screen_number = max(Screen('Screens'));
        [w, rect] = PsychImaging('OpenWindow', screen_number, bg);
        [cx, cy]  = RectCenter(rect);
        sw = rect(3); sh = rect(4);
        escapeKey = KbName('ESCAPE');

        for k = 1:numel(files)
            img = imread(fullfile(stim, files(k).name));
            tex = Screen('MakeTexture', w, img);

            scale = min(sw / size(img, 2), sh / size(img, 1));
            dest  = CenterRectOnPointd([0 0 size(img,2)*scale size(img,1)*scale], cx, cy);

            Screen('FillRect', w, bg);
            Screen('DrawTexture', w, tex, [], dest);
            Screen('Flip', w);

            fprintf('[%d/%d] %s  (%d x %d, scale %.3f) -- any key = next, ESC = quit\n', ...
                    k, numel(files), files(k).name, size(img,2), size(img,1), scale);

            % wait for a keypress; ESC quits the slideshow
            KbReleaseWait;
            quitNow = false;
            while true
                [down, ~, kc] = KbCheck;
                if down
                    if kc(escapeKey), quitNow = true; end
                    break;
                end
                WaitSecs(0.001);
            end

            Screen('Close', tex);
            if quitNow, break; end
        end

        sca;

    catch ME
        sca;
        rethrow(ME);
    end
end
