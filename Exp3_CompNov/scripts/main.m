%==========================================================================
%                  Hybrid MST N-Back Task
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
        %%%%%%%%%%%%%%%%%%%%%%%
        % setup
        %%%%%%%%%%%%%%%%%%%%%%%
        rng('shuffle');
        Screen('Preference', 'SkipSyncTests', 1);

        p.subj_id = input('Enter subject ID (e.g., 101): ');
        base_dir = '..';
        addpath(genpath(fullfile(base_dir, 'functions')));
        p.stim_dir = fullfile(base_dir, 'stimulus/stim_final/');
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
        % inter-trial fixation point in the 1-back and 2-back
        p.fix_dot_d1    = 36;              % outer disc diameter (px)
        p.fix_dot_d2    = 12;              % central dot diameter / crosshair width (px)
        p.fix_dot_color = p.colors.black;  % disc + central dot colour
        KbName('UnifyKeyNames');
        p.keys.device = -3; % listen to all keyboard (experimenter room + test room)
        fprintf('Listening to all keyboards\n');
        disp(' ');
        KbReleaseWait(p.keys.device);

        %%%%%%%%%%%%%%%%%%%%%%%
        % run experiment
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('***Experiment begins\n\n\n');

        %%%%%%%%%%%%%%%%%%%%%%%
        % global instructions
        %%%%%%%%%%%%%%%%%%%%%%%
        instructions(p, 'welcome');

        %%%%%%%%%%%%%%%%%%%%%%%
        % % which blocks to run
        % %%%%%%%%%%%%%%%%%%%%%%%
        % b_to_run = 0; % 0 = all; [x] = specific block numbers [2 3 4]
        % if b_to_run == 0, b_seq = 1:p.nBlocks; else, b_seq = b_to_run; end
        % 
        % for b = b_seq
        %     fprintf('Block...%d\n\n', b);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     % Phase 1: 1-back
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     fprintf('   Run 1-back\n');
        %     sequence_1_back_block = subject_data.sequence_1_back(subject_data.sequence_1_back.block == b, :);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     % instruction and practice for block 1 only
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     if b == 1
        %         instructions(p, '1_back');
        %         fprintf('   Run practice\n');
        %         C_run_1_back_practice(p);
        %     end
        % 
        %     results_1_back = C_run_1_back(p, sequence_1_back_block, b);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     % rest
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     message = sprintf('Well done. Please take the next 30 seconds to rest.\n\nThe screen will go blank shortly.');
        %     DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
        %     task_end_flip = Screen('Flip', p.window);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     % save 1-back data
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     try
        %         block_filename = sprintf('sub%03d_1_back_b%d.mat', p.subj_id, b);
        %         block_filepath = fullfile(p.results_dir, block_filename);
        %         save(block_filepath, 'results_1_back');
        %         fprintf('1-back block %d data saved.\n', b);
        %     catch ME
        %         warning('Could not save 1-back data for block %d. Reason: %s', b, ME.message);
        %     end
        % 
        %     WaitSecs('UntilTime', task_end_flip + 2);  % keep rest message up for >= 2 s
        % 
        %     fprintf('Rest started... (30 seconds)\n');
        %     Screen('Flip', p.window);
        %     WaitSecs(30);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     % Phase 1: 2-back
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     fprintf('   Running 2-back\n\n');
        %     sequence_2_back_block = subject_data.sequence_2_back(subject_data.sequence_2_back.block == b, :);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     % instruction and practice for only block 1
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     if b == 1
        %         instructions(p, '2_back');
        %         fprintf('   Run practice\n');
        %         D_run_2_back_practice(p);
        %     end
        % 
        %     results_2_back = D_run_2_back(p, sequence_2_back_block, b);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     % rest
        %     %%%%%%%%%%%%%%%%%%%%%%%
        %     if b < p.nBlocks
        %         message = sprintf('Well done. Please take the next 30 seconds to rest.\n\nThe screen will go blank shortly.');
        %         DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
        %         task_end_flip = Screen('Flip', p.window);
        % 
        %         %%%%%%%%%%%%%%%%%%%%%%%
        %         % save 2-back data
        %         %%%%%%%%%%%%%%%%%%%%%%%
        %         try
        %             block_filename = sprintf('sub%03d_2_back_b%d.mat', p.subj_id, b);
        %             block_filepath = fullfile(p.results_dir, block_filename);
        %             save(block_filepath, 'results_2_back');
        %             fprintf('2-back block %d data saved.\n', b);
        %         catch ME
        %             warning('SAVE_FAILED: Could not save 2-back data for block %d. Reason: %s', b, ME.message);
        %         end
        % 
        %         WaitSecs('UntilTime', task_end_flip + 2);  % keep rest message up for >= 2 s
        %         fprintf('Rest started... (30 seconds)\n');
        %         Screen('Flip', p.window); % Blank screen
        %         WaitSecs(30);
        %     else
        %         % LAST block, so just save the data without a rest
        %         try
        %             block_filename = sprintf('sub%03d_2_back_b%d.mat', p.subj_id, b);
        %             block_filepath = fullfile(p.results_dir, block_filename);
        %             save(block_filepath, 'results_2_back');
        %             fprintf('2-back block %d data saved.\n', b);
        %         catch ME
        %             warning('Could not save 2-back data for block %d. Reason: %s', b, ME.message);
        %         end
        %     end
        % end % block loop ends
        % 
        %  %%%%%%%%%%%%%%%%%%%%%%%
        %  % save 1-back and 2-back data
        %  %%%%%%%%%%%%%%%%%%%%%%%
        % results_1_back_all = consolidate_data(p, '1_back');
        % results_2_back_all = consolidate_data(p, '2_back');
        % final_data_output.subj_id = p.subj_id;
        % final_data_output.parameters = p;
        % final_data_output.results_1_back_all = results_1_back_all;
        % final_data_output.results_2_back_all = results_2_back_all;
        % save(final_data_filename, 'final_data_output');
        % fprintf('All Phase 1 data saved to:\n%s\n', final_data_filename);

        %%%%%%%%%%%%%%%%%%%%%%%
        % Phase 2: Recognition
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Running Recognition Task\n\n');
        sequence_recognition = subject_data.sequence_recognition;

        results_recognition = E_run_recognition(p, sequence_recognition);
        instructions(p, 'goodbye');

        %%%%%%%%%%%%%%%%%%%%%%%
        % save recognitoin data
        %%%%%%%%%%%%%%%%%%%%%%%
        try
            rec_filename = sprintf('sub%03d_rec.mat', p.subj_id);
            rec_filepath = fullfile(p.results_dir, rec_filename);
            save(rec_filepath, 'results_recognition');
            fprintf('Recognition data saved.\n');
        catch ME
            warning('SAVE_FAILED: Could not save recognition data. Reason: %s', ME.message);
        end

        %%%%%%%%%%%%%%%%%%%%%%%
        % save full data
        %%%%%%%%%%%%%%%%%%%%%%%
        final_data_output.results_recognition = results_recognition;
        save(final_data_filename, 'final_data_output');
        fprintf('All data (Phase 1 + Phase 2) saved to:\n%s\n', final_data_filename);

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

%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%
function [full_data_table] = consolidate_data(p, task_name)
    full_data_table = table();
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
