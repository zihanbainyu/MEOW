%==========================================================================
%              instructions & practice (before the MEG session)
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
% Same for every participant: no subject ID, no setup file, nothing saved.
%==========================================================================
function run_practice()
    addpath(genpath(fullfile('/Users/Shared/Psychtoolbox')));
    clear; clc; sca; Priority(0); ListenChar(0); ShowCursor;

    try
        %%%%%%%%%%%%%%%%%%%%%%%
        % setup
        %%%%%%%%%%%%%%%%%%%%%%%
        rng('shuffle');
        Screen('Preference', 'SkipSyncTests', 1);

        base_dir = '..';
        p.stim_dir = fullfile(base_dir, 'stimulus/stim_pool/');
        p.keys.quit = 'escape';
        p.timing.image_dur = 1.5;          % matches the MEG runs

        %%%%%%%%%%%%%%%%%%%%%%%
        % psychtoolbox
        %%%%%%%%%%%%%%%%%%%%%%%
        PsychDefaultSetup(2);
        screens = Screen('Screens');
        screen_number = max(screens);
        p.colors.white   = WhiteIndex(screen_number);
        p.colors.black   = BlackIndex(screen_number);
        p.colors.bgcolor = [124/255 124/255 124/255];
        p.text_size = 26;
        [p.window, p.windowRect] = PsychImaging('OpenWindow', screen_number, p.colors.bgcolor);
        Screen('TextSize', p.window, p.text_size);
        Screen('TextFont', p.window, 'Helvetica');
        Screen('BlendFunction', p.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
        Priority(MaxPriority(p.window));
        HideCursor(screen_number);
        p.screenX = p.windowRect(3);
        p.screenY = p.windowRect(4);
        p.centerX = p.screenX/2;
        p.centerY = p.screenY/2;
        [p.xCenter, p.yCenter] = RectCenter(p.windowRect);
        p.ifi = Screen('GetFlipInterval', p.window);
        p.fix_dot_d1    = 36;              % outer disc diameter (px), Thaler et al. (2013) ABC target
        p.fix_dot_d2    = 12;              % central dot diameter / crosshair width (px)
        p.fix_dot_color = p.colors.black;
        KbName('UnifyKeyNames');
        p.keys.device = -3;               % listen to all keyboards
        KbReleaseWait(p.keys.device);
        % Block keystrokes from reaching MATLAB (editor / command window) for the
        % whole session; KbQueue/KbCheck still read the keyboard. Restored by
        % ListenChar(0) in the clean-up below, also after an error.
        ListenChar(-1);

        %%%%%%%%%%%%%%%%%%%%%%%
        % 1-back: intro -> practice -> post-practice recap (with example images)
        %%%%%%%%%%%%%%%%%%%%%%%
        instructions(p, 'oneback');
        fprintf('   Run 1-back practice\n');
        C_run_1_back_practice(p);
        instructions(p, 'oneback_pp');

        %%%%%%%%%%%%%%%%%%%%%%%
        % 2-back: intro -> practice -> post-practice recap (with example images)
        %%%%%%%%%%%%%%%%%%%%%%%
        instructions(p, 'twoback');
        fprintf('   Run 2-back practice\n');
        D_run_2_back_practice(p);
        instructions(p, 'twoback_pp');

        instructions(p, 'final');

    catch ME
        fprintf(2, '\n! AN ERROR OCCURRED: %s !\n', ME.message);
        disp(getReport(ME, 'extended'));
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    %%%%%%%%%%%%%%%%%%%%%%%
    Priority(0);
    ListenChar(0);
    sca;
    ShowCursor;
    fprintf('\nPractice done.\n');
end
