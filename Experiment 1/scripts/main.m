%======================================================================= ===
%                  MAIN EXPERIMENTAL SCRIPT
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================

function main()
    clear; 
    clc; 
    sca;        
    
    %% ========================================================================
    %  SET UP
    %  ========================================================================
    rng('shuffle');
    Screen('Preference', 'SkipSyncTests', 1);
    
    % --- define Directories ---
    base_dir = '..';
    addpath(genpath(fullfile(base_dir, 'functions')));

    % --- get Subject Info ---
    p.subj_id = input('Enter subject ID (e.g., 101): ');
    p.eyetracking = input('Eyetracking? (1=yes, 0=no): ');
    
    % --- create data directory for this subject ---
    p.setup_dir = fullfile(base_dir, 'subj_setup');
    p.results_dir  = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
    if ~exist(p.results_dir, 'dir'), mkdir(p.results_dir); end

    % --- load pre-generated schedule ---
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

    %% ========================================================================
    %  PSYCHTOOLBOX & EYELINK SETUP
    %  ========================================================================
    try
        % --- PTB Screen Setup ---
        PsychDefaultSetup(2);
        screen_number = max(Screen('Screens'));
        [p.window, p.windowRect] = PsychImaging('OpenWindow', screen_number, WhiteIndex(screen_number) / 2);
                 [p.xCenter, p.yCenter] = RectCenter(p.windowRect);
        p.ifi = Screen('GetFlipInterval', p.window);
        p.text_size = 32;
        Screen('TextSize', p.window, p.text_size);
        p.colors.white=WhiteIndex(screen_number);
        p.colors.black=BlackIndex(screen_number);
        p.fix_cross_size = 40;
        p.fix_cross_width = 4;
        p.keys = subject_data.parameters.keys;
        p.timing = subject_data.parameters.timing;
        p.stim_dir = subject_data.parameters.stim_dir;
        p.nBlocks = subject_data.parameters.nBlocks;

        % --- Eyelink Setup ---
        el = []; % initialize empty in case of no eyetracking
        if p.eyetracking == 1
            el = Eyelink_setup(p);
        end

        %% ========================================================================
        %   MAIN EXPERIMENT
        %  ========================================================================
        HideCursor(screen_number);
        
        % --- global instructions ---
        instructions(p, 'welcome');

        blocks_to_run = [1];  % 'full' / [#] for testing
        if strcmp(blocks_to_run, 'all')
            block_sequence = 1:p.nBlocks;
        else
            block_sequence = blocks_to_run;
        end

        for b = block_sequence %
            fprintf('\n\n===== STARTING BLOCK %d of %d =====\n', b, p.nBlocks);

            %==================================================================
            % Encoding
            %==================================================================
            
            if b == 1
                instructions(p, 'encoding');
                fprintf('--- Running Encoding Practice ---\n');
                load(fullfile(p.setup_dir, 'practice_encoding_schedule.mat'), 'practice_schedule');
                C_run_encoding_practice(p, practice_schedule);
                fprintf('\n\n===== Completed Encoding Practice =====\n');   
            end

            fprintf('--- Running Encoding ---\n');
            encoding_schedule_block = subject_data.encoding_schedule(subject_data.encoding_schedule.block == b, :);
            C_run_encoding(p, encoding_schedule_block);
            
            %==================================================================
            % Rest
            %==================================================================
            
            fprintf('--- Starting Rest ---\n');
            rest(p, b, p.nBlocks);

            %==================================================================
            % Recalibration (if needed)
            %==================================================================
            
            if p.eyetracking == 1
                fprintf('--- Checking Calibration ---\n');
                ask_for_recalibration(p, el);
            end
            
            %==================================================================
            % Memory
            %==================================================================
            
            if b == 1
                instructions(p, 'memory');
                fprintf('--- Running Memory Practice ---\n');
                load(fullfile(p.setup_dir, 'practice_memory_schedule.mat'), 'practice_schedule');
                D_run_memory_practice(p, practice_schedule);
                fprintf('\n\n===== Completed Memory Practice =====\n');
            end

            fprintf('--- Running Memory ---\n');
            test_schedule_block = subject_data.test_schedule(subject_data.test_schedule.block == b, :);
            D_run_memory(p, test_schedule_block);

        end % blook loop ends

        % 
        %% ========================================================================
        %   SAVE DATA
        %  ========================================================================
        encoding_results_all = consolidate_data(p, 'enc');
        memory_results_all = consolidate_data(p, 'mem');

        final_data_output.subj_id = p.subj_id;
        final_data_output.parameters = subject_data.parameters;
        final_data_output.encoding_results = encoding_results_all;
        final_data_output.test_results = memory_results_all;

        save(final_data_filename, 'final_data_output');
        fprintf('All data saved to:\n%s\n', final_data_filename);

        instructions(p, 'goodbye');
        if p.eyetracking == 1, Eyelink('Shutdown'); end
        sca; ShowCursor;
        
        fprintf('\nThe End.\n');
        
    catch ME %
        sca; ShowCursor;
        if p.eyetracking == 1, Eyelink('Shutdown'); end
        fprintf(2, '\n! AN ERROR OCCURRED !\n'); % Print in red
        
        if strcmp(ME.identifier, 'USER_ABORT')
             fprintf('Experiment gracefully aborted by user.\n');
        else
            rethrow(ME); %
        end
    end
end

%% ========================================================================
%  HELPER FUNCTIONS FOR MAIN SCRIPT
%  ========================================================================

function el = Eyelink_setup(p)
    dummymode = 0; % set to 1 to run without a tracker connected
    if Eyelink('Initialize') ~= 0, error('Eyelink initialization failed.'); end
    
    el = EyelinkInitDefaults(p.window);
    if ~EyelinkInit(dummymode), error('Eyelink Init aborted.'); end
    
    edf_filename = sprintf('s%03d.edf', p.subj_id);
    Eyelink('OpenFile', edf_filename);
    
    EyelinkDoTrackerSetup(el); % Run calibration
    Eyelink('StartRecording');
end

function rest(p, current_block, total_blocks)
    rest_text = sprintf(['Block %d of %d is complete.\n\n' ...
                         'You may take a 2-minute rest.\n\n' ...
                         'The screen will go blank shortly.'], current_block, total_blocks);
    DrawFormattedText(p.window, rest_text, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    WaitSecs(10);

    fprintf('Rest period started... (120 seconds)\n');
    Screen('Flip', p.window); %
    WaitSecs(10);
    fprintf('Rest period finished.\n');
end

function ask_for_recalibration(p, el)
    recal_text = ['The experimenter may now check eye tracker calibration.\n\n' ...
                  'Press ''C'' to re-calibrate, or SPACE to continue.'];
    DrawFormattedText(p.window, recal_text, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);

    KbName('UnifyKeyNames');
    while true
        [~, ~, keyCode] = KbCheck;
        if keyCode(KbName('space')), break; end
        if keyCode(KbName('c'))
            Eyelink('StopRecording');
            EyelinkDoTrackerSetup(el);
            Eyelink('StartRecording');
            DrawFormattedText(p.window, recal_text, 'center', 'center', p.colors.black);
            Screen('Flip', p.window);
        end
    end
end

function [full_data_table] = consolidate_data(p, task_name)
    full_data_table = table();
    fprintf('appending %s data...\n', task_name);
    
    for b = 1:p.nBlocks
        block_filename = sprintf('sub%03d_%s_b%d.mat', p.subj_id, task_name, b);
        block_filepath = fullfile(p.results_dir, block_filename);
        
        if exist(block_filepath, 'file')
            loaded_data = load(block_filepath);
            if isfield(loaded_data, 'block_results')
                full_data_table = [full_data_table; loaded_data.block_results];
            elseif isfield(loaded_data, 'results_table')
                 full_data_table = [full_data_table; loaded_data.results_table];
            end
        else
            warning('Data file not found: %s. Final data set will be incomplete.', block_filename);
        end
    end
end