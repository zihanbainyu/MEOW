%==========================================================================
%                  MAIN EXPERIMENTAL SCRIPT
%==========================================================================
% author: Zihan Bai

function main()

    clear; 
    clc;
    sca;
    
    %% Initial setup
    rng('shuffle');

    Screen('Preference', 'SkipSyncTests', 1);
    
    % define directory
    base_dir = '..';
    setup_dir = fullfile(base_dir, 'subj_setup/');
    data_dir = fullfile(base_dir, 'data/'); 
    addpath(genpath('functions'));
    
    if ~exist(data_dir, 'dir')
        mkdir(data_dir); 
    end

    %% Prompt of subject information and eyetracking
    subj_id = [];
    while isempty(subj_id)
        try
            subj_id_str = input('subject id? (e.g., 101): ', 's');
            subj_id = str2double(subj_id_str);
            if isnan(subj_id), subj_id = []; error('not a number'); end
        catch
            fprintf('invalid subject ID. please enter a number.\n');
        end
    end

    p.eyetracking = input('eyetracking? (1=yes, 0=no): ');

    %% Pre-load subject setup file
    schedule_filename = fullfile(setup_dir, sprintf('sub%03d.mat', subj_id));
    if ~exist(schedule_filename, 'file')
        error('must run subject_setup first', subj_id);
    end

    load(schedule_filename);
    
    % prevent overwrite
    data_filename = fullfile(data_dir, sprintf('sub%03d.mat', subj_id));
    if exist(data_filename, 'file')
        overwrite = input('data already exists. overwrite? (y/n): ', 's');
        if ~strcmpi(overwrite, 'y')
            fprintf('experiment aborted.\n');
            return;
        end
    end                                      

    encoding_data = [];
    test_data = [];

    %% Setup psychotoolbox
    try
        % general setup
        PsychDefaultSetup(2); % sets color range to 0-1, unifies key names
        
        % find screen and set parameters
        % ALWAYS SET TO SCREEN 1
        screen_number = max(Screen('Screens'));
        white = WhiteIndex(screen_number);
        black = BlackIndex(screen_number);
        grey = white / 2;
        
        % open an on-screen window
        [w, wRect] = PsychImaging('OpenWindow', screen_number, grey);
        
        % get screen center and flip interval (for timing)
        [p.xCenter, p.yCenter] = RectCenter(wRect);
        p.ifi = Screen('GetFlipInterval', w);
        
        % set text parameters for instructions
        Screen('TextFont', w, 'Open Sans');
        Screen('TextSize', w, 26);
        
        % put key PTB parameters into our main struct 'p' for easy passing ---
        p.window = w;
        p.windowRect = wRect;
        p.colors.white = white;
        p.colors.black = black;
        p.colors.grey = grey;
        
        % combine parameters from setup file with PTB parameters ---
        p.subj_id = subject_data.subj_id;
        p.keys = subject_data.parameters.keys;
        p.timing = subject_data.parameters.timing;
        p.stim_dir = subject_data.parameters.stim_dir;

        %% Set up eyelink 
        if p.eyetracking ==1
            dummymode=0; 
            
            % note - need to update IP address 
            Eyelink('SetAddress', '100.1.1.1');
        
            % initialize eye tracker
            if Eyelink('initialize') ~= 0
                fprintf('Error in connecting to the eye tracker\n\n');
                Eyelink('Shutdown');
                sca;
                ListenChar(0);
                return;
            end
            
            el=EyelinkInitDefaults(window);
            el.debugPrint = true;
            el.displayCalResults = 1;
            el.callback = [];
            if ~EyelinkInit(dummymode)
                fprintf('Eyelink Init aborted.\n');
                Eyelink('Shutdown');
                sca;
                ListenChar(0);
                return;
            end
            
            % Ensure eye tracker is connected
            if Eyelink('IsConnected') ~= 1
                error('Eyetracker lost!!!')
            end
        
        end

        
        %% Main experiment
        HideCursor(screen_number);
        
        %------------------------------------------------------------------
        % Instructions
        %------------------------------------------------------------------
        % instructions(p, 'welcome');
        
        %------------------------------------------------------------------
        % Encoding Phase
        %------------------------------------------------------------------
        % creat an edf file
        if p.eyetracking == 1
            if Eyelink('IsConnected') ~= 1
                error('Eyetracker lost!!!')
            end
            % create data file
            edf_filename_enc = sprintf('s%03de.edf', subj_id); % e for encoding
            Eyelink('OpenFile', edf_filename_enc);
            EyelinkDoTrackerSetup(el);
        end
        
        % run behavioral task
        [encoding_data] = run_encoding(p, subject_data.encoding_schedule);

        if p.eyetracking == 1 
            % stop, close and safe the file
            if Eyelink('IsConnected') ~= 1
                error('Eyetracker lost!!!')
            end
            Eyelink('StopRecording');
            Eyelink('CloseFile');
            edf_dir_enc = fullfile(data_dir, edf_filename_enc);
            try Eyelink('ReceiveFile', edf_filename_enc, edf_dir_enc);
                fprintf('saved encoding edf data to: %s\n', edf_dir_enc);
            catch
                fprintf('failed to save encdoing edf data.\n');
            end
        end

        fprintf('\nphase 1 complete. saving intermediate data...\n');
        save(data_filename, 'subject_data', 'encoding_data', 'test_data');

        %------------------------------------------------------------------
        % Rest
        %------------------------------------------------------------------
        % instructions(p, 'rest');
        
        %------------------------------------------------------------------
        % Test
        %------------------------------------------------------------------
        % instructions(p, 'start_test');

        % if p.eyetracking
        %     fprintf('preparing test eye tracking...\n');
        %     % option to re-calibrate before the test phase
        %     DrawFormattedText(p.window, 'The experimenter may now check the eye tracker calibration.\n\nPlease stay still.', 'center', 'center', p.colors.white);
        %     Screen('Flip', p.window);
        %     
        %     edf_filename_test = sprintf('s%03dt.edf', subj_id); % t for test phase
        %     Eyelink('OpenFile', edf_filename_test);
        %     EyelinkDoTrackerSetup(el);
        %     Eyelink('StartRecording');
        % end

        % run test phase
        [test_data] = run_test(p, subject_data.test_schedule);

        % if p.eyetracking
        %     % stop, close and safe the file
        %     Eyelink('StopRecording')   jkj k  b            ;
        %     Eyelink('CloseFile');
        %     edf_dir_test = fullfile(data_dir, edf_filename_test);
        %     try Eyelink('ReceiveFile', edf_filename_enc, edf_dir_test);
        %         fprintf('saved encoding edf data to: %s\n, edf_dir_test');
        %     catch
        %         fprintf('failed to save test edf data.\n');
        %     end
        % end

        fprintf('\nphase 2 complete. Saving final data...\n');
        save(data_filename, 'subject_data', 'encoding_data', 'test_data');
        
        %------------------------------------------------------------------
        % End
        %------------------------------------------------------------------
        if p.eyetracking == 1
            if Eyelink('IsConnected') ~= 1
                error('Eyetracker lost!!!')
            end
            Eyelink('Shutdown');
        end
        
        % instructions(p, 'goodbye');
        
        %% Save data and clean up
        % gather all data into a final structure
        final_data_output.subj_id = p.subj_id;
        final_data_output.parameters = subject_data.parameters; % save original parameters
        final_data_output.encoding_results = encoding_data;
        final_data_output.test_results = test_data;
        
        % save the final data to a .mat file
        save(data_filename, 'final_data_output');
        fprintf('\saved all data to:\n%s\n', data_filename);
        
        % show mouse, close psychtoolbox screen, and close script
        ShowCursor;
        sca;
        fprintf('\nThe End.\n');
        
    catch ME % this 'catch' block will run if anything goes wrong
        
        fprintf('\nERROR !!!\n');
        sca;
        ShowCursor;
        
        % check if this was an intentional abort by the user
        if strcmp(ME.identifier, 'USER_ABORT:ExperimentAborted')
            fprintf('Experiment was gracefully aborted by the user.\n');
   
            save(data_filename, 'subject_data', 'encoding_data', 'test_data');
            fprintf('Partial data (if any) and schedule saved to: %s\n', data_filename);
            
        else % This was an unexpected crash
            fprintf('An unexpected error occurred.\n');
            % also try to save whatever data we might have
            save(data_filename, 'subject_data', 'encoding_data', 'test_data');
            fprintf('Partial data (if any) saved to: %s\n', data_filename);
            
            % see the full bug report
            psychrethrow(ME);
        end
        
    end
end