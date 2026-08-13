%==========================================================================
%                  MST N-Back Task
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================

function main()
    addpath(genpath(fullfile('/Users/Shared/Psychtoolbox')));
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
        p.eyetracking = input('Eyetracking? (1=yes, 0=no): ');
        base_dir = '..';
        addpath(genpath(fullfile(base_dir, 'functions')));
        p.stim_dir = fullfile(base_dir, 'stimulus/stim_pool/');
        p.instr_dir = fullfile(base_dir, 'stimulus');   % instruction slides (run_instructions)
        p.setup_dir = fullfile(base_dir, 'subj_setup');
        p.results_dir  = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
        if ~exist(p.results_dir, 'dir'), mkdir(p.results_dir); end
        final_data_filename = fullfile(p.results_dir, sprintf('sub%03d_concat.mat', p.subj_id));
        % --- Robust saving: mirror every result into a per-session timestamped
        %     backup folder, so re-running a subject never overwrites or loses a
        %     previous session's data. ---
        p.backup_dir = fullfile(base_dir, 'data_backup', ...
            sprintf('sub%03d_%s', p.subj_id, datestr(now, 'yyyymmdd_HHMMSS')));
        if ~exist(p.backup_dir, 'dir'), mkdir(p.backup_dir); end
        fprintf('Data folder:   %s\nBackup folder: %s\n', p.results_dir, p.backup_dir);

        setup_filename = fullfile(base_dir, 'subj_setup', sprintf('sub%03d_setup.mat', p.subj_id));
        load(setup_filename, 'subject_data');

        % additional parameters
        p.keys = subject_data.parameters.keys;
        p.timing = subject_data.parameters.timing;
        % NB: do NOT inherit subject_data.parameters.stim_dir. A setup file
        % generated on Windows stores a backslash path (e.g. '..\stimulus\
        % stim_pool') that dir()/fullfile() cannot resolve on macOS/Linux,
        % which leaves the practice stimulus list empty and crashes with
        % "Index exceeds array bounds". Keep the local path set at line 27.
        p.nBlocks = subject_data.parameters.nBlocks;
        sequence_1_back_all = subject_data.sequence_1_back;
        sequence_2_back_all = subject_data.sequence_2_back;
        

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
        p.fix_cross_size = 30;
        p.fix_cross_width = 4;
        p.fix_dot_d1    = 36;              % outer disc diameter (px), Thaler et al. (2013) ABC target
        p.fix_dot_d2    = 12;              % central dot diameter / crosshair width (px)
        p.fix_dot_color = p.colors.black;  % disc + central dot colour
        KbName('UnifyKeyNames');
        p.keys.device = -3; % listen to all keyboard (experimenter room + test room)
        fprintf('Listening to all keyboards\n');
        disp(' ');
        KbReleaseWait(p.keys.device);
    
        %%%%%%%%%%%%%%%%%%%%%%%
        % eyelink
        %%%%%%%%%%%%%%%%%%%%%%%
        el = []; % in case no eyetracking
        if p.eyetracking == 1
            dummymode = 0; % set to 1 for debugging without a tracker
            % Point at this Host PC's IP (default is 100.1.1.2). Must be set
            % before EyelinkInit.
            Eyelink('SetAddress', '192.168.1.5');
            el=EyelinkInitDefaults(p.window);
            if ~EyelinkInit(dummymode)
                fprintf('Eyelink Ignit aborted.\n');
                error('EYELINK_INIT_FAILED'); %
            end
            [v vs]=Eyelink('GetTrackerVersion');
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
            el.backgroundcolour = p.colors.bgcolor;
            el.foregroundcolour = p.colors.black;
            el.calibrationtargetcolour = p.colors.black;
            EyelinkUpdateDefaults(el);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%
        % run experiment
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('***Experiment begins\n\n\n');

        % %%%%%%%%%%%%%%%%%%%%%%%
        % % instructions & practice (run during the anatomical scan)
        % %%%%%%%%%%%%%%%%%%%%%%%
        % do_practice = 1;   % set 0 to skip (e.g., if practised outside the scanner)
        % if do_practice
        %     run_instructions(p, {'ins_start','ins_1','ins_2','ins_3','ins_4','ins_5','ins_prac_1back'});
        %     fprintf('   Run 1-back practice\n');
        %     C_run_1_back_practice(p);
        %     run_instructions(p, {'ins_6','ins_7','ins_8','ins_9','ins_10','ins_prac_2back'});
        %     fprintf('   Run 2-back practice\n');
        %     D_run_2_back_practice(p);
        %     run_instructions(p, {'ins_11'});
        % end
        % 
        % %%%%%%%%%%%%%%%%%%%%%%%
        % % eye-tracker calibration
        % %%%%%%%%%%%%%%%%%%%%%%%
        % if p.eyetracking == 1
        %     instructions(p, 'calibration');
        %     fprintf('Performing initial calibration\n');
        %     EyelinkDoTrackerSetup(el);
        % end

        %%%%%%%%%%%%%%%%%%%%%%%
        % which blocks to run
        %%%%%%%%%%%%%%%%%%%%%%%
        b_to_run = 0; % 0 = all; [x] = specific block numbers [2 3 4]
        if b_to_run == 0, b_seq = 1:p.nBlocks; else, b_seq = b_to_run; end

        edf_to_transfer = {};

        for b = b_seq
            fprintf('Block...%d\n\n', b);

            %%%%%%%%%%%%%%%%%%%%%%%
            % Part 1: 1-back
            %%%%%%%%%%%%%%%%%%%%%%%
            fprintf('   Run 1-back\n');
            sequence_1_back_block = subject_data.sequence_1_back(subject_data.sequence_1_back.block == b, :);

            %%%%%%%%%%%%%%%%%%%%%%%
            % eyetracking version
            %%%%%%%%%%%%%%%%%%%%%%%
            if p.eyetracking == 1
                edf_filename = sprintf('%d_1_%d.edf', p.subj_id, b);
                Eyelink('OpenFile', edf_filename);
                fprintf('EYELINK: opened edf file: %s\n', edf_filename);
                Eyelink('command', 'add_file_preamble_text ''1_Back, Block %d''', b);

                results_1_back = C_run_1_back(p, el, sequence_1_back_block, b);
            else
            %%%%%%%%%%%%%%%%%%%%%%%
            % behavior-only version
            %%%%%%%%%%%%%%%%%%%%%%%
                results_1_back = C_run_1_back(p, el, sequence_1_back_block, b);
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            % save 1-back data
            %%%%%%%%%%%%%%%%%%%%%%%
            if p.eyetracking == 1
                Eyelink('CloseFile');
                edf_to_transfer{end+1} = edf_filename;
            end
            block_filepath = fullfile(p.results_dir, sprintf('sub%03d_1_back_b%d.mat', p.subj_id, b));
            robust_save(block_filepath, 'results_1_back', results_1_back, p.backup_dir);
            fprintf('1-back block %d data saved.\n', b);

            %%%%%%%%%%%%%%%%%%%%%%%
            % rest
            %%%%%%%%%%%%%%%%%%%%%%%
            rest_dur = 30;
            message = sprintf('Fantastic job!\n\nYou are halfway through this block.\n\nPlease use the next 30 seconds to relax.\n\nYou can have you eyes open or closed.');
            DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
            rest_onset = Screen('Flip', p.window);
            fprintf('Rest... (%d s)\n', rest_dur);
            WaitSecs('UntilTime', rest_onset + rest_dur);

            %%%%%%%%%%%%%%%%%%%%%%%
            % optional recalibration
            %%%%%%%%%%%%%%%%%%%%%%%
            if p.eyetracking == 1
                fprintf('Checking Calibration\n');
                ask_for_recalibration(p, el);
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            % Part 1: 2-back
            %%%%%%%%%%%%%%%%%%%%%%%
            fprintf('   Running 2-back\n\n');
            sequence_2_back_block = subject_data.sequence_2_back(subject_data.sequence_2_back.block == b, :);

            %%%%%%%%%%%%%%%%%%%%%%%
            % eyetracking version
            %%%%%%%%%%%%%%%%%%%%%%%
            if p.eyetracking == 1
                edf_filename = sprintf('%d_2_%d.edf', p.subj_id, b);
                fprintf('EYELINK: opening edf file: %s\n', edf_filename);
                Eyelink('OpenFile', edf_filename);
                Eyelink('command', 'add_file_preamble_text ''2_Back, Block %d''', b);
                results_2_back = D_run_2_back(p, el, sequence_2_back_block, b);
            else
            %%%%%%%%%%%%%%%%%%%%%%%
            % behavior-only version
            %%%%%%%%%%%%%%%%%%%%%%%
                results_2_back = D_run_2_back(p, el, sequence_2_back_block, b);
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            % rest
            %%%%%%%%%%%%%%%%%%%%%%%
            if b < p.nBlocks
                rest_dur = 60;  
                message = sprintf('Fantastic job!\n\nYou have completed this block.\n\nPlease use the next 1 minute to relax.\n\nYou can have you eyes open or closed.');
                DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
                rest_onset = Screen('Flip', p.window);

                %%%%%%%%%%%%%%%%%%%%%%%
                % save 2-back data
                %%%%%%%%%%%%%%%%%%%%%%%
                if p.eyetracking == 1
                    Eyelink('CloseFile');
                    edf_to_transfer{end+1} = edf_filename;
                end
                block_filepath = fullfile(p.results_dir, sprintf('sub%03d_2_back_b%d.mat', p.subj_id, b));
                robust_save(block_filepath, 'results_2_back', results_2_back, p.backup_dir);
                fprintf('2-back block %d data saved.\n', b);

                fprintf('Rest... (%d s)\n', rest_dur);
                WaitSecs('UntilTime', rest_onset + rest_dur);

                %%%%%%%%%%%%%%%%%%%%%%%
                % optional recalibration
                %%%%%%%%%%%%%%%%%%%%%%%
                if p.eyetracking == 1
                    fprintf('Checking Calibration\n');
                    ask_for_recalibration(p, el);
                end
            else

                if p.eyetracking == 1
                    Eyelink('CloseFile');
                    edf_to_transfer{end+1} = edf_filename;
                end
                block_filepath = fullfile(p.results_dir, sprintf('sub%03d_2_back_b%d.mat', p.subj_id, b));
                robust_save(block_filepath, 'results_2_back', results_2_back, p.backup_dir);
                fprintf('2-back block %d data saved.\n', b);
            end
        end % block loop ends

        %%%%%%%%%%%%%%%%%%%%%%%
        % save 1-back & 2-back data
        %%%%%%%%%%%%%%%%%%%%%%%
        results_1_back_all = consolidate_data(p, '1_back');
        results_2_back_all = consolidate_data(p, '2_back');
        final_data_output.subj_id = p.subj_id;
        final_data_output.parameters = p;
        final_data_output.results_1_back_all = results_1_back_all;
        final_data_output.results_2_back_all = results_2_back_all;
        robust_save(final_data_filename, 'final_data_output', final_data_output, p.backup_dir);
        fprintf('Part 1 (1-back & 2-back) data saved to:\n%s\n', final_data_filename);

        %%%%%%%%%%%%%%%%%%%%%%%
        % Part 2: post-task MST (old / similar / new)
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n   Running Part 2: MST\n\n');
        sequence_mst = subject_data.sequence_mst;
        if p.eyetracking == 1
            fprintf('Calibration check before the MST\n');
            ask_for_recalibration(p, el);
        end
        if p.eyetracking == 1
            edf_filename = sprintf('%d_mst.edf', p.subj_id);
            Eyelink('OpenFile', edf_filename);
            fprintf('EYELINK: opened edf file: %s\n', edf_filename);
            Eyelink('command', 'add_file_preamble_text ''Post-task MST''');
            results_mst = F_run_mst(p, el, sequence_mst);
            Eyelink('CloseFile');
            edf_to_transfer{end+1} = edf_filename;
        else
            results_mst = F_run_mst(p, el, sequence_mst);
        end
        mst_filepath = fullfile(p.results_dir, sprintf('sub%03d_mst.mat', p.subj_id));
        robust_save(mst_filepath, 'results_mst', results_mst, p.backup_dir);
        fprintf('MST data saved.\n');

        %%%%%%%%%%%%%%%%%%%%%%%
        % append the MST to the concatenated file
        %%%%%%%%%%%%%%%%%%%%%%%
        final_data_output.results_mst = results_mst;
        robust_save(final_data_filename, 'final_data_output', final_data_output, p.backup_dir);
        fprintf('All data saved to:\n%s\n', final_data_filename);

        %%%%%%%%%%%%%%%%%%%%%%%
        % transfer all EDF files
        %%%%%%%%%%%%%%%%%%%%%%%
        if p.eyetracking == 1 && ~isempty(edf_to_transfer)
            fprintf('EYELINK: transferring %d EDF file(s) from Host PC...\n', numel(edf_to_transfer));
            for k = 1:numel(edf_to_transfer)
                try
                    Eyelink('ReceiveFile', edf_to_transfer{k}, p.results_dir, 1);
                    fprintf('  received %s\n', edf_to_transfer{k});
                catch ME
                    fprintf(2, '  Problem receiving %s: %s\n', edf_to_transfer{k}, ME.message);
                end
            end
        end
    
    catch ME
        fprintf(2, '\n! AN ERROR OCCURRED: %s !\n', ME.message);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    %%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%
function ask_for_recalibration(p, el)
    recal_text = 'Welcome back! Recalibrate eye tracker?\n\nPlease wait for the experimenter';
    DrawFormattedText(p.window, recal_text, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    KbReleaseWait(p.keys.device);
    while true
        [~, ~, keyCode] = KbCheck(p.keys.device);
        if keyCode(KbName('y'))
            fprintf('Recalibrating...\n');
            EyelinkDoTrackerSetup(el);
            break;
        elseif keyCode(KbName('n'))
            fprintf('Skipping recalibration.\n');
            break;
        end
        WaitSecs(0.001);
    end
    Screen('Flip', p.window);
    WaitSecs(0.5);
end

function [full_data_table] = consolidate_data(p, task_name)
    full_data_table = table();
    for b = 1:p.nBlocks
        block_filename = sprintf('sub%03d_%s_b%d.mat', p.subj_id, task_name, b);
        block_filepath = fullfile(p.results_dir, block_filename);
        if ~exist(block_filepath, 'file') && isfield(p, 'backup_dir')
            block_filepath = fullfile(p.backup_dir, block_filename);   % fall back to the backup copy
        end
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
            warning('Data file not found in primary OR backup: %s. Final set incomplete.', block_filename);
        end
    end
end

function robust_save(primary_path, varname, data, backup_dir)
% Save `data` under variable name `varname` to primary_path AND a mirror copy in
% backup_dir (the per-session timestamped folder). A failure on either target
% warns loudly but never aborts the run; a total failure is flagged CRITICAL.
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

