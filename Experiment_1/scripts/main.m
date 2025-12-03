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
        p.eyetracking = input('Eyetracking? (1=yes, 0=no): ');
        base_dir = '..';
        addpath(genpath(fullfile(base_dir, 'functions')));
        p.stim_dir = fullfile(base_dir, 'stimulus/stim_final/');
        p.setup_dir = fullfile(base_dir, 'subj_setup');
        p.results_dir  = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
        final_data_filename = fullfile(p.results_dir, sprintf('sub%03d_concat.mat', p.subj_id));
    
        % generate subj set up
        fprintf('Generating subject sequences internally...\n');
        subject_data = generate_sequences(p);
    
        % setup_filename = fullfile(base_dir, 'subj_setup', sprintf('sub%03d_setup.mat', p.subj_id));
        % load(setup_filename, 'subject_data');
        
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
        p.fix_cross_size = 30;
        p.fix_cross_width = 4;
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

        %%%%%%%%%%%%%%%%%%%%%%%
        % global instructions
        %%%%%%%%%%%%%%%%%%%%%%%
        instructions(p, 'welcome');

        %%%%%%%%%%%%%%%%%%%%%%%
        % eye-tracker calibration
        %%%%%%%%%%%%%%%%%%%%%%%
        if p.eyetracking == 1
            instructions(p, 'calibration');
            fprintf('Performing initial calibration\n');
            EyelinkDoTrackerSetup(el);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%
        % which blocks to run
        %%%%%%%%%%%%%%%%%%%%%%%
        b_to_run = 0; % 0 = all; [x] = specific block number
        if b_to_run == 0, b_seq = 1:p.nBlocks; else, b_seq = b_to_run; end
    
        for b = b_seq
            fprintf('Block...%d\n\n', b);

            %%%%%%%%%%%%%%%%%%%%%%%
            % Phase 1: 1-back
            %%%%%%%%%%%%%%%%%%%%%%%
            fprintf('   Run 1-back\n');
            sequence_1_back_block = subject_data.sequence_1_back(subject_data.sequence_1_back.block == b, :);

            %%%%%%%%%%%%%%%%%%%%%%%
            % instruction and practice for block 1 only
            %%%%%%%%%%%%%%%%%%%%%%%
            if b == 1
                instructions(p, '1_back');
                fprintf('   Run practice\n');
                C_run_1_back_practice(p);
            end

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
            % rest
            %%%%%%%%%%%%%%%%%%%%%%%
            message = sprintf('Well done!\n\nPlease use the next 45 seconds to relax.\n\nThe screen will go blank shortly');
            DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
            task_end_flip = Screen('Flip', p.window);

            %%%%%%%%%%%%%%%%%%%%%%%
            % save 1-back data
            %%%%%%%%%%%%%%%%%%%%%%%
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
                warning('Could not save 1-back data for block %d. Reason: %s', b, ME.message);
            end

            WaitSecs('UntilTime', task_end_flip);

            fprintf('Rest started... (45 seconds)\n');
            Screen('Flip', p.window);
            WaitSecs(45);

            %%%%%%%%%%%%%%%%%%%%%%%
            % optional recalibration
            %%%%%%%%%%%%%%%%%%%%%%%
            if p.eyetracking == 1
                fprintf('Checking Calibration\n');
                ask_for_recalibration(p, el);
            end
    
            %%%%%%%%%%%%%%%%%%%%%%%
            % Phase 1: 2-back
            %%%%%%%%%%%%%%%%%%%%%%%
            fprintf('   Running 2-back\n\n');
            sequence_2_back_block = subject_data.sequence_2_back(subject_data.sequence_2_back.block == b, :);
    
            %%%%%%%%%%%%%%%%%%%%%%%
            % instruction and practice for only block 1
            %%%%%%%%%%%%%%%%%%%%%%%
            if b == 1
                instructions(p, '2_back');
                fprintf('   Run practice\n');
                D_run_2_back_practice(p);
            end
    
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
                message = sprintf('Fantastic job!\n\nPlease use the next 1 minute to relax.\n\nThe screen will go blank shortly');
                DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
                task_end_flip = Screen('Flip', p.window);
    
                %%%%%%%%%%%%%%%%%%%%%%%
                % save 2-back data
                %%%%%%%%%%%%%%%%%%%%%%%
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
    
                %%%%%%%%%%%%%%%%%%%%%%%
                % optional recalibration
                %%%%%%%%%%%%%%%%%%%%%%%
                if p.eyetracking == 1
                    fprintf('Checking Calibration\n');
                    ask_for_recalibration(p, el);
                end
            else
                % LAST block, so just save the data without a rest
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
                    warning('Could not save 2-back data for block %d. Reason: %s', b, ME.message);
                end
            end
        end % block loop ends
    
         %%%%%%%%%%%%%%%%%%%%%%%
         % save 1-back and 2-back data
         %%%%%%%%%%%%%%%%%%%%%%%
        results_1_back_all = consolidate_data(p, '1_back');
        results_2_back_all = consolidate_data(p, '2_back');
        final_data_output.subj_id = p.subj_id;
        final_data_output.parameters = p;
        final_data_output.results_1_back_all = results_1_back_all;
        final_data_output.results_2_back_all = results_2_back_all;
        save(final_data_filename, 'final_data_output');
        fprintf('All Phase 1 data saved to:\n%s\n', final_data_filename);
    
        %%%%%%%%%%%%%%%%%%%%%%%
        % Phase 2: Recognition
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Running Recognition Task\n\n');
        sequence_recognition = subject_data.sequence_recognition;
    
        %%%%%%%%%%%%%%%%%%%%%%%
        % Phase 1: 2-back
        %%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%
        % eyetracking version
        %%%%%%%%%%%%%%%%%%%%%%%
        if p.eyetracking == 1
            edf_filename = sprintf('%d_rec.edf', p.subj_id);
            Eyelink('OpenFile', edf_filename);
            fprintf('EYELINK: opened edf file: %s\n', edf_filename);
            Eyelink('command', 'add_file_preamble_text ''Recog''');
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
        %%%%%%%%%%%%%%%%%%%%%%%
        % behavior-only version
        %%%%%%%%%%%%%%%%%%%%%%%
            results_recognition = E_run_recognition(p, el, sequence_recognition);
            instructions(p, 'goodbye');
        end
    
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

function subject_data = generate_sequences(p)
    subj_id = p.subj_id;
    sequence_1_back = []; sequence_2_back = []; goal_list_full = [];
    p.counts.oneback = struct('resp_same', 0, 'resp_similar', 0, 'resp_new', 0);
    p.counts.twoback = struct('resp_same', 0, 'resp_similar', 0, 'resp_new', 0);
    p.nComparison = 120; p.nIsolated_Both = 120; p.nNovel = 120;
    p.nTotalPairs = p.nComparison + p.nIsolated_Both + p.nNovel;
    p.nBlocks = 4; p.keys.same = 'j'; p.keys.diff = 'k'; p.keys.quit = 'escape';
    p.timing.image_dur = 1.5; p.timing.fix_dur = 0.75; p.timing.fix_jitter = 0.25;
    
    fprintf('loading stimuli...\n');
    all_A_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_A_l1.png'));
    all_B_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_B_l1.png'));
    A_names_l1 = string({all_A_files_l1.name}'); B_names_l1 = string({all_B_files_l1.name}');
    master_pair_list_l1 = table(sort(A_names_l1), sort(B_names_l1), 'VariableNames', {'A', 'B'});
    
    all_A_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_A_l2.png'));
    all_B_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_B_l2.png'));
    A_names_l2 = string({all_A_files_l2.name}'); B_names_l2 = string({all_B_files_l2.name}');
    master_pair_list_l2 = table(sort(A_names_l2), sort(B_names_l2), 'VariableNames', {'A', 'B'});
    
    all_foil_A_files = dir(fullfile(p.stim_dir, 'mst_*_A_foil.png'));
    all_foil_B_files = dir(fullfile(p.stim_dir, 'mst_*_B_foil.png'));
    all_foil_pairs = table(string({all_foil_A_files.name}'), string({all_foil_B_files.name}'), 'VariableNames', {'A_foil', 'B_foil'});
    all_foil_pairs = all_foil_pairs(randperm(height(all_foil_pairs)), :);
    
    n_cond_l1 = height(master_pair_list_l1)/3; n_cond_l2 = height(master_pair_list_l2)/3;
    final_list_l1 = master_pair_list_l1(randperm(height(master_pair_list_l1)), :);
    comp_l1 = final_list_l1(1:n_cond_l1, :); iso_l1 = final_list_l1(n_cond_l1+1:2*n_cond_l1, :); nov_l1 = final_list_l1(2*n_cond_l1+1:end, :);
    final_list_l2 = master_pair_list_l2(randperm(height(master_pair_list_l2)), :);
    comp_l2 = final_list_l2(1:n_cond_l2, :); iso_l2 = final_list_l2(n_cond_l2+1:2*n_cond_l2, :); nov_l2 = final_list_l2(2*n_cond_l2+1:end, :);
    
    comp_pairs = [comp_l1; comp_l2]; comp_pairs = comp_pairs(randperm(height(comp_pairs)), :);
    iso_pairs = [iso_l1; iso_l2]; iso_pairs = iso_pairs(randperm(height(iso_pairs)), :);
    novel_pairs = [nov_l1; nov_l2]; novel_pairs = novel_pairs(randperm(height(novel_pairs)), :);
    
    p.stim.comp_l1=comp_l1; p.stim.comp_l2=comp_l2; p.stim.iso_l1=iso_l1; p.stim.iso_l2=iso_l2; p.stim.novel_l1=nov_l1; p.stim.novel_l2=nov_l2;
    p.stim.compared=comp_pairs; p.stim.isolated=iso_pairs; p.stim.novel=novel_pairs;
    
    partition_idx = @(N, nb) arrayfun(@(k) ((floor((k-1)*N/nb)+1):floor(k*N/nb)), 1:nb, 'UniformOutput', false);
    p.block_indices.comp_l1 = partition_idx(height(p.stim.comp_l1), p.nBlocks);
    p.block_indices.iso_l1 = partition_idx(height(p.stim.iso_l1), p.nBlocks);
    p.block_indices.nov_l1 = partition_idx(height(p.stim.novel_l1), p.nBlocks);
    p.block_indices.comp_l2 = partition_idx(height(p.stim.comp_l2), p.nBlocks);
    p.block_indices.iso_l2 = partition_idx(height(p.stim.iso_l2), p.nBlocks);
    p.block_indices.nov_l2 = partition_idx(height(p.stim.novel_l2), p.nBlocks);
    
    all_foils_remain = all_foil_pairs;
    for b = 1:p.nBlocks
        fprintf('Block %d \n', b);
        comp_l1_b = p.stim.comp_l1(p.block_indices.comp_l1{b}, :);
        comp_l2_b = p.stim.comp_l2(p.block_indices.comp_l2{b}, :);
        comp_pairs_b = [comp_l1_b; comp_l2_b]; comp_pairs_b = comp_pairs_b(randperm(height(comp_pairs_b)), :);
        iso_l1_b = p.stim.iso_l1(p.block_indices.iso_l1{b}, :);
        iso_l2_b = p.stim.iso_l2(p.block_indices.iso_l2{b}, :);
        iso_pairs_b = [iso_l1_b; iso_l2_b]; iso_pairs_b = iso_pairs_b(randperm(height(iso_pairs_b)), :);
        nov_l1_b = p.stim.novel_l1(p.block_indices.nov_l1{b}, :);
        nov_l2_b = p.stim.novel_l2(p.block_indices.nov_l2{b}, :);
        novel_pairs_b = [nov_l1_b; nov_l2_b]; novel_pairs_b = novel_pairs_b(randperm(height(novel_pairs_b)), :);
        
        nComp_b = height(comp_pairs_b); nIso_b = height(iso_pairs_b); nNov_b = height(novel_pairs_b);
        comp_miniblocks = {}; repeat_miniblocks = {}; iso_trials = {};
        for i = 1:nComp_b
            comp_miniblocks{end+1} = {comp_pairs_b.A(i), "compared", "A", "none"; comp_pairs_b.B(i), "compared", "B", "k"};
        end
        n_foil_repeats = nComp_b;
        if height(all_foils_remain) < n_foil_repeats, error('Not enough foils 1-back'); end
        repeat_foil_pairs = all_foils_remain(1:n_foil_repeats, :); all_foils_remain(1:n_foil_repeats, :) = [];
        for i = 1:n_foil_repeats
            foil_item = repeat_foil_pairs.A_foil(i);
            repeat_miniblocks{end+1} = {foil_item, "repeat", "A", "none"; foil_item, "repeat", "A", "j"};
        end
        for i = 1:nIso_b
            iso_trials{end+1} = {iso_pairs_b.A(i), "isolated", "A", "none"}; iso_trials{end+1} = {iso_pairs_b.B(i), "isolated", "B", "none"};
        end
        miniblocks_1_back = [comp_miniblocks, repeat_miniblocks, iso_trials];
        final_1_back_list = vertcat(miniblocks_1_back{randperm(numel(miniblocks_1_back))});
        sequence_1_back_block = cell2table(final_1_back_list, 'VariableNames', {'stim_id', 'condition', 'identity', 'corr_resp'});
        sequence_1_back_block.block = repmat(b, height(sequence_1_back_block), 1);
        
        for i = 1:height(iso_pairs_b)
            tgt = iso_pairs_b.A(i); lur = iso_pairs_b.B(i);
            idx_tgt = find(sequence_1_back_block.stim_id == tgt, 1); idx_lur = find(sequence_1_back_block.stim_id == lur, 1);
            if ~isempty(idx_tgt) && ~isempty(idx_lur) && idx_lur < idx_tgt
                [sequence_1_back_block(idx_tgt, :), sequence_1_back_block(idx_lur, :)] = deal(sequence_1_back_block(idx_lur, :), sequence_1_back_block(idx_tgt, :));
            end
        end
        sequence_1_back = [sequence_1_back; sequence_1_back_block];
        p.counts.oneback.resp_same = p.counts.oneback.resp_same + sum(strcmp(sequence_1_back_block.corr_resp, 'j'));
        p.counts.oneback.resp_similar = p.counts.oneback.resp_similar + sum(strcmp(sequence_1_back_block.corr_resp, 'k'));
        p.counts.oneback.resp_new = p.counts.oneback.resp_new + sum(strcmp(sequence_1_back_block.corr_resp, 'none'));
    
        block_pairs = [comp_pairs_b; iso_pairs_b; novel_pairs_b];
        block_conditions = [repmat("compared", nComp_b, 1); repmat("isolated", nIso_b, 1); repmat("novel", nNov_b, 1)];
        goal_list = table(block_pairs.A, block_pairs.B, block_conditions, repmat("", 90, 1), repmat("", 90, 1), 'VariableNames', {'A', 'B', 'condition', 'goal_type', 'X'});
        goal_types = [repmat("A-B", 10, 1); repmat("A-A", 10, 1); repmat("A-N", 10, 1)];
        goal_list.goal_type(1:30) = goal_types(randperm(30)); goal_list.goal_type(31:60) = goal_types(randperm(30)); goal_list.goal_type(61:90) = goal_types(randperm(30));
        for i = 1:90
            if goal_list.goal_type(i)=="A-B", goal_list.X(i)=goal_list.B(i); elseif goal_list.goal_type(i)=="A-A", goal_list.X(i)=goal_list.A(i); end
        end
        goal_list = goal_list(randperm(height(goal_list)), :);
        goal_list_b = goal_list; goal_list_b.block = repmat(b, height(goal_list_b), 1);
        goal_list_full = [goal_list_full; goal_list_b];
    
        sequence = cell(300, 5); row_idx = 1;
        junk1 = all_foils_remain(1,:); junk2 = all_foils_remain(2,:); all_foils_remain(1:2,:) = [];
        sequence(row_idx,:) = {junk1.A_foil, "init_junk", "J", "JUNK", "none"}; row_idx=row_idx+1;
        sequence(row_idx,:) = {junk2.A_foil, "init_junk", "J", "JUNK", "none"}; row_idx=row_idx+1;
        active_goals = []; goals_started = false(height(goal_list), 1); goal_pointer = 1;
        
        while goal_pointer <= height(goal_list) || ~isempty(active_goals)
            to_complete = [];
            for k = 1:size(active_goals, 1)
                if row_idx - active_goals(k, 2) == 2, to_complete(end+1) = k; end
            end
            if ~isempty(to_complete)
                k = to_complete(1); goal_idx = active_goals(k, 1);
                goal = goal_list(goal_idx, :); X_type = goal.goal_type;
                if strcmp(X_type, "A-N")
                    next_unstarted = [];
                    for check_idx = 1:height(goal_list)
                        if ~goals_started(check_idx) && ~ismember(check_idx, active_goals(:, 1)), next_unstarted = check_idx; break; end
                    end
                    if ~isempty(next_unstarted), X = goal_list.A(next_unstarted); else, X = all_foils_remain.A_foil(1); all_foils_remain(1,:) = []; end
                    X_props = goal_list(strcmp(goal_list.A, X), :);
                    if isempty(X_props), log_cond="junk_foil"; log_goal="JUNK"; X_id="J"; else, log_cond=X_props.condition(1); log_goal=X_props.goal_type(1); X_id="A"; end
                    X_resp = "none";
                else
                    X = goal.X; log_cond = goal.condition; log_goal = X_type;
                    if strcmp(X_type, "A-B"), X_resp="k"; X_id="B"; else, X_resp="j"; X_id="A"; end
                end
                sequence(row_idx,:) = {X, log_cond, X_id, log_goal, X_resp}; row_idx=row_idx+1; active_goals(k, :) = [];
                if strcmp(X_type, "A-N")
                    X_goal_idx = find(strcmp(goal_list.A, X), 1);
                    if ~isempty(X_goal_idx) && ~goals_started(X_goal_idx), goals_started(X_goal_idx)=true; active_goals(end+1, :) = [X_goal_idx, row_idx-1]; end
                end
            else
                while goal_pointer <= height(goal_list) && goals_started(goal_pointer), goal_pointer = goal_pointer + 1; end
                if goal_pointer <= height(goal_list)
                    goal = goal_list(goal_pointer, :); sequence(row_idx,:) = {goal.A, goal.condition, "A", goal.goal_type, "none"};
                    goals_started(goal_pointer) = true; active_goals(end+1, :) = [goal_pointer, row_idx]; row_idx=row_idx+1; goal_pointer=goal_pointer+1;
                else
                    if isempty(active_goals), break; end
                    if height(all_foils_remain)<1, error('Out of foils'); end
                    sequence(row_idx,:) = {all_foils_remain.A_foil(1), "padding_junk", "J", "JUNK", "none"}; all_foils_remain(1,:) = []; row_idx=row_idx+1;
                end
            end
        end
        if height(all_foils_remain) >= 5
            for jj = 1:5, sequence(row_idx,:) = {all_foils_remain.A_foil(jj), "end_junk", "J", "JUNK", "none"}; row_idx=row_idx+1; end
            all_foils_remain(1:5,:) = [];
        end
        sequence_2_back_block = cell2table(sequence(1:row_idx-1, :), 'VariableNames', {'stim_id', 'condition', 'identity', 'goal','corr_resp'});
        sequence_2_back_block.block = repmat(b, height(sequence_2_back_block), 1);
        sequence_2_back = [sequence_2_back; sequence_2_back_block];
        p.counts.twoback.resp_same = p.counts.twoback.resp_same + sum(strcmp(sequence_2_back_block.corr_resp, 'j'));
        p.counts.twoback.resp_similar = p.counts.twoback.resp_similar + sum(strcmp(sequence_2_back_block.corr_resp, 'k'));
        p.counts.twoback.resp_new = p.counts.twoback.resp_new + sum(strcmp(sequence_2_back_block.corr_resp, 'none'));
    end
    
    tested_goals = goal_list_full(goal_list_full.goal_type ~= "A-N" & goal_list_full.condition ~= "novel", :);
    n_tested = height(tested_goals); selected_items = strings(n_tested, 1); selected_identity = strings(n_tested, 1);
    for i = 1:n_tested
        if rand() > 0.5, selected_items(i) = tested_goals.A(i); selected_identity(i) = "A"; else, selected_items(i) = tested_goals.B(i); selected_identity(i) = "B"; end
    end
    all_old_items = table(selected_items, tested_goals.condition, selected_identity, repmat("old", n_tested, 1), repmat(p.keys.same, n_tested, 1), ...
        'VariableNames', {'stim_id', 'condition', 'identity', 'trial_type', 'corr_resp'});
    n_rec_foils = 160;
    assert(height(all_foils_remain) >= n_rec_foils, 'Not enough foils for recognition');
    new_foils_list = all_foils_remain.A_foil(1:n_rec_foils); all_foils_remain(1:n_rec_foils, :) = [];
    all_new_foils = table(new_foils_list, repmat("foil", n_rec_foils, 1), repmat("N", n_rec_foils, 1), repmat("new", n_rec_foils, 1), repmat(p.keys.diff, n_rec_foils, 1), ...
        'VariableNames', {'stim_id', 'condition', 'identity', 'trial_type', 'corr_resp'});
    sequence_recognition = [all_old_items; all_new_foils]; sequence_recognition = sequence_recognition(randperm(height(sequence_recognition)), :);
    
    sequence_1_back.fix_duration = p.timing.fix_dur + (rand(height(sequence_1_back), 1) * 2 - 1) * p.timing.fix_jitter;
    sequence_2_back.fix_duration = p.timing.fix_dur + (rand(height(sequence_2_back), 1) * 2 - 1) * p.timing.fix_jitter;
    sequence_recognition.fix_duration = repmat(0.5, height(sequence_recognition), 1);
    sequence_1_back.subj_id = repmat(subj_id, height(sequence_1_back), 1);
    sequence_2_back.subj_id = repmat(subj_id, height(sequence_2_back), 1);
    sequence_recognition.subj_id = repmat(subj_id, height(sequence_recognition), 1);
    
    subject_data.subj_id = subj_id; subject_data.parameters = p;
    subject_data.sequence_1_back = sequence_1_back; subject_data.sequence_2_back = sequence_2_back; subject_data.sequence_recognition = sequence_recognition;
    output_filename = fullfile(p.setup_dir, sprintf('sub%03d_setup.mat', subj_id));
    save(output_filename, 'subject_data');
end