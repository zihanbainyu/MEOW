%==========================================================================
%                  MAIN EXPERIMENTAL SCRIPT
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================
function main()
clear;
clc;
sca; 
Priority(0);
ListenChar(0);
ShowCursor;

try
    %% Set up
    rng('shuffle');
    Screen('Preference', 'SkipSyncTests', 1);

    % directories
    base_dir = '..';
    addpath(genpath(fullfile(base_dir, 'functions')));

    % subj info
    p.subj_id = input('Enter subject ID (e.g., 101): ');
    p.eyetracking = input('Eyetracking? (1=yes, 0=no): ');
    p.setup_dir = fullfile(base_dir, 'subj_setup');
    p.results_dir  = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
    if ~exist(p.results_dir, 'dir'), mkdir(p.results_dir); end

    % load subj set up
    setup_filename = fullfile(base_dir, 'subj_setup', sprintf('sub%03d_setup.mat', p.subj_id));
    if ~exist(setup_filename, 'file')
        error('you have to run A_subject_setup.m first.');
    end
    load(setup_filename, 'subject_data');
    final_data_filename = fullfile(p.results_dir, sprintf('sub%03d_concat.mat', p.subj_id));
    if exist(final_data_filename, 'file')
        if ~strcmpi(input('final data file already exists. Overwrite? (y/n): ', 's'), 'y')
            fprintf('Experiment aborted to prevent overwriting.\n'); return;
        end
    end

    % additional parameters
    p.keys = subject_data.parameters.keys;
    p.timing = subject_data.parameters.timing;
    p.stim_dir = subject_data.parameters.stim_dir;
    p.nBlocks = subject_data.parameters.nBlocks;
    sequence_1_back_all = subject_data.sequence_1_back;
    sequence_2_back_all = subject_data.sequence_2_back;
    sequence_recognition = subject_data.sequence_recognition;

    %% Psychtoolbox
    % PTB screen
    PsychDefaultSetup(2);

    % project to secondary screen
    screens = Screen('Screens');
    screen_number = max(screens);
    p.colors.white=WhiteIndex(screen_number);
    p.colors.black=BlackIndex(screen_number);
    p.colors.bgcolor=[124/255 124/255 124/255];
    p.text_size = 26;
    [p.window, p.windowRect] = PsychImaging('OpenWindow', screen_number, p.colors.bgcolor);

    % font and blend settings
    Screen('TextSize', p.window, p.text_size);
    Screen('TextFont', p.window, 'Helvetica');
    Screen('BlendFunction', p.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    % set priority and hide cursor
    Priority(MaxPriority(p.window));
    HideCursor(screen_number); 
    % ListenChar(-1); 
    p.screenX = p.windowRect(3);
    p.screenY = p.windowRect(4);
    p.centerX = p.screenX/2;
    p.centerY = p.screenY/2;
    [p.xCenter, p.yCenter] = RectCenter(p.windowRect);
    p.ifi = Screen('GetFlipInterval', p.window);
    p.fix_cross_size = 30;
    p.fix_cross_width = 4;

    % keyboard connection
    KbName('UnifyKeyNames');
    p.keys.device = -3; % listen to all keyboard (experimenter room + test room)
    fprintf('Listening to all keyboards\n');
    disp(' ');
    KbReleaseWait(p.keys.device);

    %% Eyelink set up
    el = []; % in case no eyetracking
    if p.eyetracking == 1
        dummymode = 0; % set to 1 for debugging without a tracker
        el=EyelinkInitDefaults(p.window);
        if ~EyelinkInit(dummymode)
            fprintf('Eyelink Ignit aborted.\n');
            error('EYELINK_INIT_FAILED'); %
        end

        [v vs]=Eyelink('GetTrackerVersion');

        Eyelink('command', 'set_idle_mode');
        WaitSecs(0.05);

        % eyelink configuration
        Eyelink('command', 'sample_rate = 1000');
        Eyelink('command', 'pupil_model = ELLIPSE');
        Eyelink('command', 'pupil_size_diameter = YES');
        [width, height] = Screen('WindowSize', screen_number);
        Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
        Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
        Eyelink('command', 'calibration_type = HV9');
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
       
        if sscanf(vs(12:end),'%f') >= 4
            Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,PUPIL,HREF,VELOCITY,HTARGET,GAZERES,STATUS,INPUT');
            Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,PUPIL,VELOCITY,GAZERES,HTARGET,STATUS,INPUT');
        else
            Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,PUPIL,HREF,VELOCITY,GAZERES,STATUS,INPUT');
            Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,PUPIL,VELOCITY,GAZERES,STATUS,INPUT');
        end
        Eyelink('command', 'button_function 5 "accept_target_fixation"');

        % make sure we're still connected
        if Eyelink('IsConnected')~=1 && dummymode == 0
            error('EYELINK_CONNECTION_LOST');
        end

        % set up calibration
        el.backgroundcolour = p.colors.bgcolor;
        el.foregroundcolour = p.colors.black;
        el.calibrationtargetcolour = p.colors.black;
        EyelinkUpdateDefaults(el);
    end

    %% Run Experiment
    
    fprintf('***Experiment begins\n\n\n');

    % global instructions
    % instructions(p, 'welcome');
    % 
    % % initial eyetracker calibration
    % if p.eyetracking == 1
    %     instructions(p, 'calibration');
    %     fprintf('Performing initial calibration\n');
    %     EyelinkDoTrackerSetup(el);
    % end

    % run all or some?
        % if the task crashes for some reason (although it should not), set to
        % the block you want to run to continue
    b_to_run = [3 4];
    if b_to_run == 0, b_seq = 1:p.nBlocks; else, b_seq = b_to_run; end

    for b = b_seq
        fprintf('Start Block...%d\n\n', b);


        %% Phase 1: run 1-back
        fprintf('   Run 1-back\n');
        sequence_1_back_block = subject_data.sequence_1_back(subject_data.sequence_1_back.block == b, :);

        % instruction and practice: 1-back for only block 1
        if b == 1
            instructions(p, '1_back'); 
            fprintf('   Run practice\n');
            C_run_1_back_practice(p);
        end


        if p.eyetracking == 1
            edf_filename = sprintf('%d_1_%d.edf', p.subj_id, b);
            Eyelink('OpenFile', edf_filename);
            fprintf('EYELINK: opened edf file: %s\n', edf_filename);
            Eyelink('command', 'add_file_preamble_text ''1_Back, Block %d''', b);

            % Eye-tracking version
            results_1_back = C_run_1_back(p, el, sequence_1_back_block, b);

            % fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
            % Eyelink('CloseFile');
            % WaitSecs(0.01);
            % try
            %     Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
            % catch ME
            %     fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
            % end
        else

            % Behavior-only version
            results_1_back = C_run_1_back(p, el, sequence_1_back_block, b);
        end

        %% rest
        message = sprintf('Well done!\n\nPlease use the next 45 seconds to relax.\n\nThe screen will go blank shortly');
        DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
        task_end_flip = Screen('Flip', p.window);

        %% save data during rest
        if p.eyetracking == 1
            fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
            Eyelink('CloseFile');
            WaitSecs(0.1);
            try
                Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
            catch ME
                fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
            end
        end

        try
            block_filename = sprintf('sub%03d_1_back_b%d.mat', p.subj_id, b);
            block_filepath = fullfile(p.results_dir, block_filename);
            save(block_filepath, 'results_1_back');
            fprintf('1-back block %d data saved.\n', b);
        catch ME
            warning('SAVE_FAILED: Could not save 1-back data for block %d. Reason: %s', b, ME.message);
        end

        WaitSecs('UntilTime', task_end_flip); 

        fprintf('Rest started... (45 seconds)\n');
        Screen('Flip', p.window);
        WaitSecs(45);

        if p.eyetracking == 1
            fprintf('Checking Calibration\n');
            ask_for_recalibration(p, el);
        end

        % try
        %     block_filename = sprintf('sub%03d_1_back_b%d.mat', p.subj_id, b);
        %     block_filepath = fullfile(p.results_dir, block_filename);
        %     save(block_filepath, 'results_1_back');
        %     fprintf('1-back block %d data saved.\n', b);
        % catch ME
        %     warning('SAVE_FAILED: Could not save 1-back data for block %d. Reason: %s', b, ME.message);
        % end


        %% Phase 1: run 2-back
        fprintf('   Running 2-back\n\n');
        sequence_2_back_block = subject_data.sequence_2_back(subject_data.sequence_2_back.block == b, :);

        % instruction and practice: 2-back for only block 1
        if b == 1
            instructions(p, '2_back'); 
            fprintf('   Run practice\n');
            D_run_2_back_practice(p);
        end

        if p.eyetracking == 1
            edf_filename = sprintf('%d_2_%d.edf', p.subj_id, b);
            fprintf('EYELINK: opening edf file: %s\n', edf_filename);
            Eyelink('OpenFile', edf_filename);
            Eyelink('command', 'add_file_preamble_text ''2_Back, Block %d''', b);

            % Eye-tracking version
            results_2_back = D_run_2_back(p, el, sequence_2_back_block, b);
            % fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
            % Eyelink('CloseFile');
            % WaitSecs(0.1);
            % try
            %     Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
            % catch ME
            %     fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
            % end
        else
            % Behavior-only version
            results_2_back = D_run_2_back(p, el, sequence_2_back_block, b);
        end

        if b < p.nBlocks
            message = sprintf('Fantastic job!\n\nPlease use the next 1 minute to relax.\n\nThe screen will go blank shortly');
            DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
            task_end_flip = Screen('Flip', p.window);

            if p.eyetracking == 1
                fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
                Eyelink('CloseFile');
                WaitSecs(0.1);
                try
                    Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
                catch ME
                    fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
                end
            end
            
            try
                block_filename = sprintf('sub%03d_2_back_b%d.mat', p.subj_id, b);
                block_filepath = fullfile(p.results_dir, block_filename);
                save(block_filepath, 'results_2_back');
                fprintf('2-back block %d data saved.\n', b);
            catch ME
                warning('SAVE_FAILED: Could not save 2-back data for block %d. Reason: %s', b, ME.message);
            end

            WaitSecs('UntilTime', task_end_flip);
            fprintf('Rest started... (60 seconds)\n');
            Screen('Flip', p.window); % Blank screen
            WaitSecs(60);
            
            if p.eyetracking == 1
                fprintf('Checking Calibration\n');
                ask_for_recalibration(p, el);
            end
        else
            % LAST block, so just save the data without a rest period
            if p.eyetracking == 1
                fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
                Eyelink('CloseFile');
                WaitSecs(0.1);
                try
                    Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
                catch ME
                    fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
                end
            end
            
            try
                block_filename = sprintf('sub%03d_2_back_b%d.mat', p.subj_id, b);
                block_filepath = fullfile(p.results_dir, block_filename);
                save(block_filepath, 'results_2_back');
                fprintf('2-back block %d data saved.\n', b);
            catch ME
                warning('SAVE_FAILED: Could not save 2-back data for block %d. Reason: %s', b, ME.message);
            end
        end

        % save behavioral data
        % try
        %     block_filename = sprintf('sub%03d_2_back_b%d.mat', p.subj_id, b);
        %     block_filepath = fullfile(p.results_dir, block_filename);
        %     save(block_filepath, 'results_2_back');
        %     fprintf('2-back block %d data saved.\n', b);
        % catch ME
        %     warning('SAVE_FAILED: Could not save 2-back data for block %d. Reason: %s', b, ME.message);
        % end

        

    end % block loop ends

    %% Save Phase 1 data
    results_1_back_all = consolidate_data(p, '1_back');
    results_2_back_all = consolidate_data(p, '2_back');
    final_data_output.subj_id = p.subj_id;
    final_data_output.parameters = p;
    final_data_output.results_1_back_all = results_1_back_all;
    final_data_output.results_2_back_all = results_2_back_all;
    save(final_data_filename, 'final_data_output');
    fprintf('All Phase 1 data saved to:\n%s\n', final_data_filename);

    %% Phase 2: Recognition task
    fprintf('Running Recognition Task\n\n');    
    sequence_recognition = subject_data.sequence_recognition;

    if p.eyetracking == 1
            edf_filename = sprintf('%d_rec.edf', p.subj_id);
            Eyelink('OpenFile', edf_filename);
            fprintf('EYELINK: opened edf file: %s\n', edf_filename);
            Eyelink('command', 'add_file_preamble_text ''Recog''');

            % Eye-tracking version
            results_recognition = E_run_recognition(p, el, sequence_recognition);
            instructions(p, 'goodbye');
            fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
            Eyelink('CloseFile');
            WaitSecs(0.1);
            try
                Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
            catch ME
                fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
            end
        else
            % Behavior-only version
            results_recognition = E_run_recognition(p, el, sequence_recognition);
            instructions(p, 'goodbye');
    end
    
    % save behavioral data
    try
        rec_filename = sprintf('sub%03d_rec.mat', p.subj_id);
        rec_filepath = fullfile(p.results_dir, rec_filename);
        save(rec_filepath, 'results_recognition');
        fprintf('Recognition data saved.\n');
    catch ME
        warning('SAVE_FAILED: Could not save recognition data. Reason: %s', ME.message);
    end

    %% Final save
    final_data_output.results_recognition = results_recognition;
    save(final_data_filename, 'final_data_output');
    fprintf('All data (Phase 1 + Phase 2) saved to:\n%s\n', final_data_filename);
    
catch ME
    fprintf(2, '\n! AN ERROR OCCURRED: %s !\n', ME.message);
end

%% Clean up
if exist('p', 'var') && isfield(p, 'eyetracking') && p.eyetracking == 1
    if Eyelink('IsConnected') == 1
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        Eyelink('Shutdown');
    end
end
Priority(0);
ListenChar(0);
sca;
ShowCursor;
fprintf('\nThe End.\n');
end

%% Functions
function ask_for_recalibration(p, el)
    recal_text = 'Welcome back! Recalibrate eye tracker?\n\nWait for experimenter to press y for yes, n for skip';
    DrawFormattedText(p.window, recal_text, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    KbReleaseWait(p.keys.device);
    while true
        % Use KbCheck with the specific device index
        [~, ~, keyCode] = KbCheck(p.keys.device);
        if keyCode(KbName('y'))
            fprintf('Recalibrating...\n');
            EyelinkDoTrackerSetup(el);
            break;
        elseif keyCode(KbName('n'))
            fprintf('Skipping recalibration.\n');
            break;
        end
        WaitSecs(0.001); % Yield CPU
    end
    Screen('Flip', p.window);
    WaitSecs(0.5);
end


function [full_data_table] = consolidate_data(p, task_name)
    full_data_table = table();
    fprintf('appending %s data...\n', task_name);
    for b = 1:p.nBlocks
        block_filename = sprintf('sub%03d_%s_b%d.mat', p.subj_id, task_name, b);
        block_filepath = fullfile(p.results_dir, block_filename);
        
        if exist(block_filepath, 'file')
            loaded_data = load(block_filepath);
            
            if isfield(loaded_data, 'results_1_back')
                full_data_table = [full_data_table; loaded_data.results_1_back];
            elseif isfield(loaded_data, 'results_2_back')
                full_data_table = [full_data_table; loaded_data.results_2_back];
            else
                warning('Variable (results_1_back/results_2_back) not found in file: %s.', block_filename);
            end
        else
            warning('Data file not found: %s. Final data set will be incomplete.', block_filename);
        end
    end
end