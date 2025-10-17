function p = F_run_survey(p)
% F_run_survey - Post-task feedback survey (Likert scale)
% -------------------------------------------------------------------------
% Presents a short diagnostic survey to collect feedback on task difficulty,
% clarity, and engagement. Results are saved in p.survey_results.
% -------------------------------------------------------------------------
% Author: Zihan Bai, Michelmann Lab, NYU

    % ========================== SETUP ====================================
    black   = p.colors.black;
    gray    = [0.3 0.3 0.3];
    bgcolor = p.colors.bgcolor;

    start_key  = KbName('g');
    escape_key = KbName(p.keys.quit);
    keys_1_5   = KbName({'1!', '2@', '3#', '4$', '5%'}); 

    Screen('TextSize', p.window, p.text_size);
    Screen('TextFont', p.window, 'Helvetica');
    Screen('TextStyle', p.window, 0);
    Screen('FillRect', p.window, bgcolor);

    % ========================== SURVEY ITEMS =============================
    survey.questions = {
    'How clear and helpful were the instructions provided before the task began?'
    'Across the entire experiment, how mentally demanding was the task overall?'
    'Overall, how engaging was the task (did it hold your attention)?'
    'For Part A (1-Back), how difficult did you find the comparison task?'
    'For Part A (1-Back), how accurate do you feel your responses were?'
    'For Part A (1-Back), how confident were you in the majority of your judgments?'
    'For Part B (2-Back), how difficult did you find the comparison task?'
    'For Part B (2-Back), how accurate do you feel your responses were?'
    'For Part B (2-Back), how confident were you in the majority of your judgments?'
    'Overall, how difficult was it to identify a SIMILAR item rather than an identical one?'
    };

    survey.results = nan(length(survey.questions), 1);

    % ========================== WELCOME SCREEN ===========================
    welcome_text = [
        'You''re almost done!\n\n' ...
        'As this is a pilot study, your honest feedback will directly help us refine the design.\n\n' ...
        'On the next screens, please answer several short questions using the number keys (1–5).\n' ...
        'Each question will show what the numbers mean for that item.\n\n\n' ...
        'Press "G" to begin.'
    ];

    DrawFormattedText(p.window, welcome_text, 'center', 'center', black, 70, [], 1.5);
    Screen('Flip', p.window);
    waitForKeyPress(p, start_key, escape_key);

    % ========================== QUESTION LOOP ============================
    Y_QUESTION = p.yCenter - 150;
    Y_SCALE    = p.yCenter + 150;

    for q = 1:length(survey.questions)
        Screen('FillRect', p.window, bgcolor);

        % --- Section Title ---
        if q <= 3
            section = 'Part 1: General Experience';
        elseif q <= 6
            section = 'Part 2: 1-Back Task (Part A)';
        elseif q <= 9
            section = 'Part 3: 2-Back Task (Part B)';
        else
            section = 'Final Question';
        end
        DrawFormattedText(p.window, section, 'center', 70, gray, 60);

        % --- Question Text ---
        DrawFormattedText(p.window, survey.questions{q}, 'center', Y_QUESTION, black, 70, [], 1.6);

        % --- Scale Anchors ---
        if ismember(q, [2, 4, 7, 10])
            scale_text = '1 = Extremely Easy / Low Effort        5 = Extremely Difficult / High Effort';
        else
            scale_text = '1 = Not at all (Clear / Accurate / Confident / Engaging)        5 = Extremely';
        end
        DrawFormattedText(p.window, scale_text, 'center', Y_SCALE, gray, 70, [], 1.3);

        % --- Progress Prompt ---
        progress_text = sprintf('Question %d of %d   |   Press 1–5 to continue', q, length(survey.questions));
        DrawFormattedText(p.window, progress_text, 'center', p.windowRect(4) - 80, gray, 70);

        Screen('Flip', p.window);
        KbReleaseWait(p.keys.device);

        % --- Wait for Response ---
        response_valid = false;
        while ~response_valid
            [keyIsDown, ~, keyCode] = KbCheck(p.keys.device);
            if keyIsDown
                if keyCode(escape_key)
                    error('USER_ABORT: Survey terminated early.');
                end
                key_index = find(keyCode(keys_1_5), 1);
                if ~isempty(key_index)
                    survey.results(q) = key_index;
                    response_valid = true;
                end
            end
            WaitSecs(0.001);
        end
        WaitSecs(0.25);
    end

    % ========================== SAVE RESULTS =============================
    p.survey_results = survey;
    try
        survey_filename = sprintf('sub%03d_survey.mat', p.subj_id);
        survey_filepath = fullfile(p.results_dir, survey_filename);
        save(survey_filepath, 'p');
        fprintf('\nSurvey data saved to: %s\n', survey_filepath);
    catch ME
        warning('Survey save failed: %s', ME.message);
    end

    % ========================== FINAL SCREEN =============================
    final_text = [
        'Survey complete!\n\n' ...
        'Thank you for your time and thoughtful feedback.\n' ...
        'Your input helps us make the study clearer and more engaging for future participants.\n\n' ...
        'If you have any additional comments, please share them with the experimenter now.\n\n' ...
        'Press "G" to finish.'
    ];

    Screen('FillRect', p.window, bgcolor);
    DrawFormattedText(p.window, final_text, 'center', 'center', black, 70, [], 1.5);
    Screen('Flip', p.window);
    waitForKeyPress(p, start_key, escape_key);

    % ========================== PRINT TO CONSOLE ==========================
    fprintf('\n--- SURVEY RESULTS ---\n');
    for i = 1:length(survey.questions)
        fprintf('Q%02d: %s\n → Response: %d\n', i, survey.questions{i}, survey.results(i));
    end
    fprintf('----------------------\n\n');
end
