%==========================================================================
%              final recognition test
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================
function run_recognition()
    addpath(genpath(fullfile('/Users/Shared/Psychtoolbox')));
    clear; clc; sca; Priority(0); ListenChar(0); ShowCursor;

    try
        %%%%%%%%%% 1111221122122221221111211222221112211122211111212211221211111212222222111212212122212122111112122111221111212112221111112111122122121122122222122121211221121222222%%%%%%%%%%%%
        % setup
        %%%%%%%%%%%%%%%%%%%%%%%
        rng('shuffle');
        Screen('Preference', 'SkipSyncTests', 1);

        p.subj_id = input('Enter subject ID (e.g., 102): ');
        base_dir = '..';
        addpath(genpath(fullfile(base_dir, 'functions')));
        p.stim_dir    = fullfile(base_dir, 'stimulus/stim_pool/');
        p.results_dir = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
        if ~exist(p.results_dir, 'dir'), mkdir(p.results_dir); end
        p.backup_dir  = fullfile(base_dir, 'data_backup', ...
            sprintf('sub%03d_recog_%s', p.subj_id, datestr(now, 'yyyymmdd_HHMMSS')));
        if ~exist(p.backup_dir, 'dir'), mkdir(p.backup_dir); end
        fprintf('Data folder:   %s\nBackup folder: %s\n', p.results_dir, p.backup_dir);

        setup_filename = fullfile(base_dir, 'subj_setup', sprintf('sub%03d_setup.mat', p.subj_id));
        if ~exist(setup_filename, 'file')
            error('Setup file not found: %s (run A_subject_setup.m first).', setup_filename);
        end
        load(setup_filename, 'subject_data');
        if ~isfield(subject_data, 'sequence_recognition')
            error('sequence_recognition not found in %s. Re-generate with A_subject_setup.m.', setup_filename);
        end
        p.keys   = subject_data.parameters.keys;
        p.timing = subject_data.parameters.timing;
        sequence_recognition = subject_data.sequence_recognition;

        %%%%%%%%%%%%%%%%%%%%%%%
        % psychtoolbox (laptop display; no external monitor)
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
        [p.xCenter, p.yCenter] = RectCenter(p.windowRect);
        p.ifi = Screen('GetFlipInterval', p.window);
        p.fix_dot_d1    = 36;              % outer disc diameter (px), Thaler et al. (2013) ABC target
        p.fix_dot_d2    = 12;              % central dot diameter / crosshair width (px)
        p.fix_dot_color = p.colors.black;
        KbName('UnifyKeyNames');
        p.keys.device = -3;               % listen to all keyboards (built-in laptop keyboard)
        KbReleaseWait(p.keys.device);

        %%%%%%%%%%%%%%%%%%%%%%%
        % run recognition
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('***Recognition test begins\n\n');
        results_recognition = E_run_recognition(p, sequence_recognition);

        %%%%%%%%%%%%%%%%%%%%%%%
        % save
        %%%%%%%%%%%%%%%%%%%%%%%
        rec_filepath = fullfile(p.results_dir, sprintf('sub%03d_recognition.mat', p.subj_id));
        robust_save(rec_filepath, 'results_recognition', results_recognition, p.backup_dir);
        fprintf('Recognition data saved to:\n%s\n', rec_filepath);

        % append to the concatenated file if it already exists
        concat_file = fullfile(p.results_dir, sprintf('sub%03d_concat.mat', p.subj_id));
        if exist(concat_file, 'file')
            S = load(concat_file);
            if isfield(S, 'final_data_output')
                final_data_output = S.final_data_output;
                final_data_output.results_recognition = results_recognition;
                robust_save(concat_file, 'final_data_output', final_data_output, p.backup_dir);
                fprintf('Appended recognition to:\n%s\n', concat_file);
            end
        end

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
    fprintf('\nThe End.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%
function robust_save(primary_path, varname, data, backup_dir)
% Save `data` under variable name `varname` to primary_path AND a mirror copy in
% backup_dir. A failure on either target warns but never aborts; a total failure
% is flagged CRITICAL.
    S = struct(); S.(varname) = data;
    [~, nm, ext] = fileparts(primary_path);
    targets = {primary_path};
    if nargin >= 4 && ~isempty(backup_dir)
        targets{end+1} = fullfile(backup_dir, [nm ext]);
    end
    saved_any = false;
    for i = 1:numel(targets)
        try
            save(targets{i}, '-struct', 'S');
            saved_any = true;
        catch ME
            warning('SAVE_FAILED: %s  (%s)', targets{i}, ME.message);
        end
    end
    if ~saved_any
        fprintf(2, '\n!! CRITICAL: %s could not be saved to ANY location !!\n\n', varname);
    end
end
