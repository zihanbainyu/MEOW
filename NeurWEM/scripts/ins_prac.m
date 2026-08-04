%==========================================================================
%                  Hybrid MST N-Back Task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================

function ins_prac()
    addpath(genpath('/Users/Shared/Psychtoolbox'));
    clear;
    clc;
    sca;
    Priority(0);
    ListenChar(0);
    ShowCursor;

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
    p.results_dir  = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
    final_data_filename = fullfile(p.results_dir, sprintf('sub%03d_concat.mat', p.subj_id));

    setup_filename = fullfile(base_dir, 'subj_setup', sprintf('sub%03d_setup.mat', p.subj_id));
    load(setup_filename, 'subject_data');

    % additional parameters
    p.keys = subject_data.parameters.keys;
    p.timing = subject_data.parameters.timing;
    p.stim_dir = subject_data.parameters.stim_dir;
    p.nBlocks = subject_data.parameters.nBlocks;
    sequence_1_back_all = subject_data.sequence_1_back;
    sequence_2_back_all = subject_data.sequence_2_back;
    sequence_recognition = subject_data.sequence_recognition;

    %%%%%%%%%%%%%%%%%%%%%%%
    % psychotoolbox
    %%%%%%%%%%%%%%%%%%%%%%%
    PsychDefaultSetup(2);
    screens = Screen('Screens');
    screen_number = max(screens);
    p.colors.white=WhiteIndex(screen_number);
    p.colors.black=BlackIndex(screen_number);
    p.colors.bgcolor=[124/255 124/255 124/255];
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
    p.fix_cross_size = 45;
    p.fix_cross_width = 6;
    % Thaler et al. (2013) combined fixation target, used as the
    % inter-trial fixation point in the 1-ba12fck and 2-back
    p.fix_dot_d1    = 36;              % outer disc diameter (px)
    p.fix_dot_d2    = 12;              % central dot diameter / crosshair width (px)
    p.fix_dot_color = p.colors.black;  % disc + central dot colour
    KbName('UnifyKeyNames');
    p.keys.device = -1; % listen to all keyboard (experimenter room + test room)
    fprintf('Listening to all keyboards\n');
    disp(' ');
    KbReleaseWait(p.keys.device);

    %%%%%%%%%%%%%%%%%%%%%%%
    % run experiment
    %%%%%%%%%%%%%%%%%%%%%%%
    fprintf('***Instructions begins\n\n\n');

    % overview + 1-back rules, then 1-back practice
    run_instructions(p, {'ins_start','ins_1','ins_2','ins_3','ins_4','ins_5','ins_prac_1back'});
    fprintf('   Run 1-back practice\n');
    C_run_1_back_practice(p);

    % 2-back rules, then 2-back practice
    run_instructions(p, {'ins_6','ins_7','ins_8','ins_9','ins_10','ins_prac_2back'});
    fprintf('   Run 2-back practice\n');
    D_run_2_back_practice(p);
        
    %%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    %%%%%%%%%%%%%%%%%%%%%%%
    Priority(0);
    ListenChar(0);
    sca;
    ShowCursor;
    fprintf('\n Practice End.\n');
end