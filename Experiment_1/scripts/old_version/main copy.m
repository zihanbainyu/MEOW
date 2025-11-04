%==========================================================================
%                  MAIN EXPERIMENTAL SCRIPT
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================
function main()
    % --- Initial Cleanup ---
    clear; 
    clc;  
    sca; % Close all screens
    Priority(0); % Ensure priority is reset
    ListenChar(0); % Ensure keyboard is reset
    ShowCursor; % Ensure cursor is visible
    
    % --- Use a robust try-catch block for guaranteed cleanup ---
    try
        %% ========================================================================
        %  SETUP: DIRECTORIES AND SUBJECT INFO
        %  ========================================================================
        rng('shuffle');
        % Set to 0 (or remove) for a real experiment to enforce timing checks!
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
        
        p.keys = subject_data.parameters.keys;
        p.timing = subject_data.parameters.timing;
        p.stim_dir = subject_data.parameters.stim_dir;
        p.nBlocks = subject_data.parameters.nBlocks;
        
        %% ========================================================================
        %  PSYCHTOOLBOX & EYELINK INITIALIZATION
        %  ========================================================================
        
        % --- PTB Screen Setup ---
        PsychDefaultSetup(2);
        
        % project to secondary screen
        screens = Screen('Screens');
        screen_number = max(screens);
        p.colors.white=WhiteIndex(screen_number);
        p.colors.black=BlackIndex(screen_number);
        p.colors.bgcolor=[124/255 124/255 124/255];
        p.text_size = 26;
        
        [p.window, p.windowRect] = PsychImaging('OpenWindow', screen_number, p.colors.bgcolor);
        
        % Apply font and blend settings once
        Screen('TextSize', p.window, p.text_size);
        Screen('TextFont', p.window, 'Helvetica');
        Screen('BlendFunction', p.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
        
        % Set Priority and Hide Cursor
        Priority(MaxPriority(p.window));
        HideCursor(screen_number); % Hide mouse pointer
        ListenChar(-1); % Suppress keyboard output to command window
        
        p.screenX = p.windowRect(3);
        p.screenY = p.windowRect(4);
        p.centerX = p.screenX/2;
        p.centerY = p.screenY/2;
        [p.xCenter, p.yCenter] = RectCenter(p.windowRect);
        p.ifi = Screen('GetFlipInterval', p.window);
        
        p.fix_cross_size = 30;
        p.fix_cross_width = 4;

        % --- Keyboard Identification ---
        KbName('UnifyKeyNames');
        
        % By setting the device index to -3, Psychtoolbox will listen to ALL
        % connected keyboard devices simultaneously. This is the simplest way
        % to allow both the experimenter and participant keyboards to work.
        p.keys.device = -3; 
        
        fprintf('SUCCESS: Listening to ALL keyboards (device index -3)\n');
        disp(' ');
        KbReleaseWait(p.keys.device); % This line is optional but good practice to wait for key release
               
        % --- eyelink setup ---
        el = []; % in case no eyetracking
        if p.eyetracking == 1
            dummymode = 0; % set to 1 for debugging without a tracker
            el=EyelinkInitDefaults(p.window);
    
            if ~EyelinkInit(dummymode)
                fprintf('Eyelink Ignit aborted.\n');
                error('EYELINK_INIT_FAILED'); % Throw error to trigger cleanup
            end
        
            [v vs]=Eyelink('GetTrackerVersion');
            fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        
            % --- eyelink configuration ---
            Eyelink('command', 'sample_rate = 1000');
            Eyelink('command', 'pupil_model = ELLIPSE');
            Eyelink('command', 'pupil_size_diameter = YES'); % Set mode to Diameter
            
            [width, height] = Screen('WindowSize', screen_number);
            Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
            Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);                
            
            Eyelink('command', 'calibration_type = HV13'); 
            Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
            Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');

            % Check tracker version and set data fields accordingly
            if sscanf(vs(12:end),'%f') >= 4
                % Request PUPIL for diameter, GAZE, VELOCITY and other useful measures
                Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,PUPIL,HREF,VELOCITY,HTARGET,GAZERES,STATUS,INPUT');
                Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,PUPIL,VELOCITY,GAZERES,HTARGET,STATUS,INPUT');
            else
                % Request PUPIL for diameter, GAZE, VELOCITY and other useful measures
                Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,PUPIL,HREF,VELOCITY,GAZERES,STATUS,INPUT');
                Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,PUPIL,VELOCITY,GAZERES,STATUS,INPUT');
            end
            
            Eyelink('command', 'button_function 5 "accept_target_fixation"');
           
            % make sure we're still connected
            if Eyelink('IsConnected')~=1 && dummymode == 0
                error('EYELINK_CONNECTION_LOST'); % Throw error to trigger cleanup
            end
           % --- set up calibration ---
            el.backgroundcolour = p.colors.bgcolor;
            el.foregroundcolour = p.colors.black;
            el.calibrationtargetcolour = p.colors.black;
            EyelinkUpdateDefaults(el);
        end
        
        %% ========================================================================
        %   MAIN EXPERIMENT LOOP
        %  ========================================================================

        % --- global instructions ---
        instructions(p, 'welcome');
        
        % --- initial calibration here ---
        if p.eyetracking == 1
            instructions(p, 'calibration');
            fprintf(' --- Performing initial calibration ---\n');
            EyelinkDoTrackerSetup(el);
            
        end
        
        blocks_to_run = 'all';  % 'full' / [#] for testing
        if strcmp(blocks_to_run, 'all')
            block_sequence = 1:p.nBlocks;
        else
            block_sequence = blocks_to_run;
        end 
        
        for b = block_sequence 
            fprintf('\n\n===== STARTING BLOCK %d of %d =====\n', b, p.nBlocks);
            
            %==================================================================
            % Part A: Encoding
            %==================================================================
            if b == 1
                instructions(p, 'encoding');
                fprintf('--- Running Encoding Practice ---\n');
                load(fullfile(p.setup_dir, 'practice_encoding_schedule.mat'), 'practice_schedule');
                C_run_encoding_practice(p, practice_schedule);
            end

            fprintf('--- Running Encoding ---\n');
            encoding_schedule_block = subject_data.encoding_schedule(subject_data.encoding_schedule.block == b, :);

            if p.eyetracking == 1
                edf_filename = sprintf('%d_e_%d.edf', p.subj_id, b);
                Eyelink('OpenFile', edf_filename);
                fprintf('EYELINK: opened edf file: %s\n', edf_filename);
                Eyelink('command', 'add_file_preamble_text ''Encoding, Block %d''', b);

                C_run_encoding(p, el, encoding_schedule_block);

                fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
                Eyelink('CloseFile');
                try
                    Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
                catch ME
                    fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
                end
            else
                C_run_encoding(p, el, encoding_schedule_block);
            end

            %==================================================================
            % Mid rest
            %==================================================================
            post_enc_rest(p, el, b);

            % ==================================================================
            % Memory
            % ==================================================================

            if b == 1
                instructions(p, 'memory');
                fprintf('--- Running Memory Practice ---\n');
                load(fullfile(p.setup_dir, 'practice_memory_schedule.mat'), 'practice_schedule');
                D_run_memory_practice(p, practice_schedule);
            end
            
            fprintf('--- Running Memory ---\n');
            test_schedule_block = subject_data.test_schedule(subject_data.test_schedule.block == b, :);
            
            if p.eyetracking == 1
                edf_filename = sprintf('%d_m_%d.edf', p.subj_id, b);
                fprintf('EYELINK: opening edf file: %s\n', edf_filename);
                Eyelink('OpenFile', edf_filename);
                Eyelink('command', 'add_file_preamble_text ''Memory, Block %d''', b);

                D_run_memory(p, el, test_schedule_block); % Corrected: passed 'el'

                fprintf('EYELINK: receiving edf file: %s\n', edf_filename);
                Eyelink('CloseFile');
                try
                    Eyelink('ReceiveFile', edf_filename, p.results_dir, 1);
                catch ME
                    fprintf('Problem receiving data file ''%s'': %s\n', edf_filename, ME.message);
                end
            else 
                D_run_memory(p, el, test_schedule_block);
            end
            
            %==================================================================
            % Post Memory Rest
            %==================================================================
            if b < p.nBlocks
                post_mem_rest(p, el, b);
            end
        end % block loop ends
        
        %% ========================================================================
        %   SAVE DATA & FINAL CLEANUP
        %  ========================================================================
        instructions(p, 'goodbye');

        encoding_results_all = consolidate_data(p, 'enc');
        memory_results_all = consolidate_data(p, 'mem');
        final_data_output.subj_id = p.subj_id;
        final_data_output.parameters = subject_data.parameters;
        final_data_output.encoding_results = encoding_results_all;
        final_data_output.test_results = memory_results_all;
        save(final_data_filename, 'final_data_output');
        fprintf('All data saved to:\n%s\n', final_data_filename);

    % --- Cleanup outside the try block ---
    catch ME
        fprintf(2, '\n! AN ERROR OCCURRED: %s !\n', ME.message); % Print error in red
    end

    % --- Guaranteed Final Cleanup ---
    if exist('p', 'var') && isfield(p, 'eyetracking') && p.eyetracking == 1
        Eyelink('Shutdown'); 
    end
    Priority(0);
    ListenChar(0);
    sca; 
    ShowCursor;
    fprintf('\nThe End.\n');
end

%% ========================================================================
%  FUNCTIONS FOR MAIN SCRIPT
%  ========================================================================
function ask_for_recalibration(p, el)
    recal_text = 'Welcome back! Recalibrate eye tracker?\n\nPress y for yes, n for skip';
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

function post_enc_rest(p, el, b)
    fprintf('--- Starting Mid Rest ---\n');
    message = sprintf('You have completed Part A of Block %d.\n\nPlease use the next 1 minute to relax.\n\nThe screen will go blank shortly', b);
    DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    WaitSecs(10);
    fprintf('Rest started... (60 seconds)\n');
    Screen('Flip', p.window);
    WaitSecs(60); 
    fprintf('Rest finished.\n');
    if p.eyetracking == 1
        fprintf('--- Checking Calibration ---\n');
        ask_for_recalibration(p, el);
    end
end

function post_mem_rest(p, el, b)
    fprintf('--- Starting Rest ---\n');
    message = sprintf('You have completed Part B of Block %d.\n\nPlease use the next 2 minutes to relax.\n\nThe screen will go blank shortly', b);
    DrawFormattedText(p.window, message, 'center', 'center', p.colors.black);
    Screen('Flip', p.window);
    WaitSecs(10);
    fprintf('Rest started... (120 seconds)\n');
    Screen('Flip', p.window);
    WaitSecs(120); 
    fprintf('Rest finished.\n');
    if p.eyetracking == 1
        fprintf('--- Checking Calibration ---\n');
        ask_for_recalibration(p, el);
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