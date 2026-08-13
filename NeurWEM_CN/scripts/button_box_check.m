function button_box_check()
% BUTTON_BOX_CHECK  Confirm the Siemens Prisma button box + scanner trigger
% actually reach PsychToolbox on THIS stim computer, using the SAME device
% index the tasks use (p.keys.device = -3, "all keyboards/keypads").
%
% HOW TO USE (run before the real scan, box connected):
%   1) cd to this scripts folder, type:  button_box_check
%   2) Press each button box button: index finger, then middle finger.
%   3) If you can, fire/simulate the scanner trigger.
%   4) Confirm the printout shows:  index -> "1" | middle -> "2" | trigger -> "5"
%
% If a button shows a DIFFERENT name/code, update p.keys.same / p.keys.diff /
% p.keys.trigger in A_subject_setup.m (and the KbName lists in C/D/F_run_*.m)
% to match what the box actually sends. If a press does NOT register here at
% all, device -3 is not seeing the box -- check the device list printed below
% (the box must appear as a keyboard/keypad) and its HID/keyboard mode.
%
% Press ESC or q to finish.  Michelmann Lab @ NYU.

    KbName('UnifyKeyNames');
    device = -3;                       % MUST match p.keys.device used by the tasks
    escKey = KbName('ESCAPE');
    qKey   = KbName('q');

    fprintf('\n============================================================\n');
    fprintf(' Button-box / trigger check   (listening on device %d)\n', device);
    fprintf('============================================================\n');
    try
        [idx, names] = GetKeyboardIndices();
        fprintf('Detected keyboard/keypad devices (the box should be one of these):\n');
        for k = 1:numel(idx)
            fprintf('   device index %d : %s\n', idx(k), names{k});
        end
    catch
        fprintf('(could not enumerate devices)\n');
    end
    fprintf('\nPress: index finger, middle finger, then the trigger.\n');
    fprintf('Expect:  index -> "1"   middle -> "2"   trigger -> "5"\n');
    fprintf('ESC or q to finish.\n\n');

    try
        ListenChar(2);                 % keep the presses out of the command window
        prev = [];
        while true
            [down, secs, keyCode] = KbCheck(device);
            if down
                codes = find(keyCode);
                if ~isequal(codes, prev)
                    nm = KbName(codes);
                    if ischar(nm), nm = {nm}; end
                    line = sprintf('%8.3f s | ', secs);
                    for k = 1:numel(codes)
                        n = nm{k}; if isempty(n), n = '<unnamed>'; end
                        line = [line, sprintf('"%s" (code %d)   ', n, codes(k))]; %#ok<AGROW>
                    end
                    disp(line);
                    if any(codes == escKey) || any(codes == qKey), break; end
                    prev = codes;
                end
            else
                prev = [];
            end
            WaitSecs(0.005);
        end
    catch ME
        ListenChar(0);
        rethrow(ME);
    end
    ListenChar(0);

    fprintf('\nDone.\n');
    fprintf('You should have seen "1" and "2" from the two buttons, and "5" from the trigger.\n');
    fprintf('Anything different -> update p.keys.same / p.keys.diff / p.keys.trigger to match.\n');
end
