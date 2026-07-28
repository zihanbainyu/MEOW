%==========================================================================
%          Recognition task (independent, run separately)
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%
% Old/new recognition of the A-B compared pairs. This used to be Part 2 of
% main.m; it is now a standalone task the experimenter launches on its own,
% whenever the recognition test is scheduled (e.g. outside the scanner).
% It reads the same sub###_setup.mat and saves sub###_rec.mat.
%
%   >> main_recognition
%==========================================================================
function main_recognition()
    clear;
    clc;
    sca;
    Priority(0);
    ListenChar(0);
    ShowCursor;

    try
        %%%%%%%%%%%%%%%%%%%%%%%
        % setup
        %%%%%%%%%%%%%%%%%%%%%%%
        rng('shuffle');
        Screen('Preference', 'SkipSyncTests', 1);

        p.subj_id = input('Enter subject ID (e.g., 101): ');
        base_dir = '..';
        addpath(genpath(fullfile(base_dir, 'functions')));
        p.stim_dir = fullfile(base_dir, 'stimulus/stim_pool/');
        p.instr_dir = fullfile(base_dir, 'stimulus');
        p.setup_dir = fullfile(base_dir, 'subj_setup');
        p.results_dir = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
        if ~exist(p.results_dir, 'dir'), mkdir(p.results_dir); end

        setup_filename = fullfile(base_dir, 'subj_setup', sprintf('sub%03d_setup.mat', p.subj_id));
        load(setup_filename, 'subject_data');

        p.keys = subject_data.parameters.keys;
        p.timing = subject_data.parameters.timing;
        p.stim_dir = subject_data.parameters.stim_dir;
        sequence_recognition = subject_data.sequence_recognition;

        %%%%%%%%%%%%%%%%%%%%%%%
        % psychtoolbox
        %%%%%%%%%%%%%%%%%%%%%%%
        PsychDefaultSetup(2);
        screens = Screen('Screens');
        screen_number = max(screens);
        p.colors.white = WhiteIndex(screen_number);
        p.colors.black = BlackIndex(screen_number);
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
        p.fix_dot_d1    = 36;
        p.fix_dot_d2    = 12;
        p.fix_dot_color = p.colors.black;
        KbName('UnifyKeyNames');
        p.keys.device = -3;
        KbReleaseWait(p.keys.device);

        %%%%%%%%%%%%%%%%%%%%%%%
        % run recognition
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('*** Recognition task begins\n\n');
        results_recognition = E_run_recognition(p, sequence_recognition);
        instructions(p, 'goodbye');

        %%%%%%%%%%%%%%%%%%%%%%%
        % save
        %%%%%%%%%%%%%%%%%%%%%%%
        try
            rec_filename = sprintf('sub%03d_rec.mat', p.subj_id);
            rec_filepath = fullfile(p.results_dir, rec_filename);
            save(rec_filepath, 'results_recognition');
            fprintf('Recognition data saved to:\n%s\n', rec_filepath);
        catch ME
            warning('SAVE_FAILED: Could not save recognition data. Reason: %s', ME.message);
        end

    catch ME
        fprintf(2, '\n! AN ERROR OCCURRED: %s !\n', ME.message);
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    %%%%%%%%%%%%%%%%%%%%%%%
    Priority(0);
    ListenChar(0);
    sca;
    ShowCursor;
    fprintf('\nThe End.\n');
end
