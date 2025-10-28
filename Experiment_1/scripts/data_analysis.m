subj_ids = [901, 902, 903, 904];
base_dir = '..';
data_dir = fullfile(base_dir, 'data');


AllSubsEncoding = table();
AllSubsTest = table();

% j for same, k for similar, no response for new

for s = 1:length(subj_ids)
    subj_id = subj_ids(s);
    data_filename_template = 'sub%03d/sub%03d_concat.mat';
    load_filename = sprintf(data_filename_template, subj_id, subj_id);
    load_filepath = fullfile(data_dir, load_filename);

    if ~exist(load_filepath, 'file')
        warning('Cannot find data file for subject %d: %s. Skipping.', subj_id, load_filepath);
        continue;
    end

    fprintf('Loading file: %s\n', load_filename);
    S = load(load_filepath, 'final_data_output');

    if ~isfield(S, 'final_data_output')
        warning('File for subject %d does not contain FINAL_DATA_OUTPUT. Skipping.', subj_id);
        continue;
    end

    T_subject = S.final_data_output;

    % Extract subj_id value from T_subject (for reliability)
    if isfield(T_subject, 'subj_id')
        true_subj_id = T_subject.subj_id;
    else
        true_subj_id = subj_id; % fallback
    end

    % Process encoding results
    if isfield(T_subject, 'encoding_results') && istable(T_subject.encoding_results)
        T_enc = T_subject.encoding_results;
        T_enc.subj_id = repmat(true_subj_id, height(T_enc), 1);
        AllSubsEncoding = [AllSubsEncoding; T_enc];
    else
        warning('Subject %d has no encoding_results.', true_subj_id);
    end

    % Process test results
    if isfield(T_subject, 'test_results') && istable(T_subject.test_results)
        T_test = T_subject.test_results;
        T_test.subj_id = repmat(true_subj_id, height(T_test), 1);
        AllSubsTest = [AllSubsTest; T_test];
    else
        warning('Subject %d has no test_results.', true_subj_id);
    end
end

fprintf('Finished loading all subjects.\n');


%% DATA PRE-PROCESSING & RECODING (BASIC)
T.test  = AllSubsTest;
T.encod = AllSubsEncoding;

strVars = {'response_key','correct_response','condition','trial_type_final'};
for v = strVars
    if ~isstring(T.test.(v{1}))
        T.test.(v{1}) = string(T.test.(v{1}));
    end
end

%% BASIC SANITY CHECKS
T_correct = T.test(T.test.correct == 1 & ~isnan(T.test.rt), :);
Stats_Acc = groupsummary(T.test, {'subj_id','trial_type_final'}, 'mean', 'correct');
% 
% T_correct = T.test(T.test.correct == 1, :);
Stats_RT = groupsummary(T_correct, {'subj_id','trial_type_final'}, 'mean', 'rt');

Stats_Acc_Overall = groupsummary(T.test, 'subj_id', 'mean', 'correct');

fprintf('\n## Overall Accuracy by Subject ##\n');
disp(Stats_Acc_Overall);

% --- 2. Accuracy by Trial Type ---
Stats_Acc_byType_Overall = groupsummary(T.test, 'trial_type_final', 'mean', 'correct');

fprintf('\n## Overall Accuracy by Trial Type (Across All Subjects) ##\n');
disp(Stats_Acc_byType_Overall);

% --- 3. Overall RT per Subject (Correct Trials Only) ---
Stats_RT_Overall = groupsummary(T_correct, 'subj_id', 'mean', 'rt');

fprintf('\n## Overall RT by Subject (Correct Trials Only) ##\n');
disp(Stats_RT_Overall);

% --- 4. RT by Trial Type (Correct Trials Only) ---
% a) Overall RT for each trial type (averaged across all subjects)
Stats_RT_byType_Overall = groupsummary(T_correct, 'trial_type_final', 'mean', 'rt');

fprintf('\n## Overall RT by Trial Type (Across All Subjects, Correct Only) ##\n');
disp(Stats_RT_byType_Overall);

% b) RT for each trial type *for each subject*
fprintf('\n## RT by Trial Type & Subject (Correct Only) ##\n');
disp(Stats_RT); % 'Stats_RT' was calculated in your step 4


%% Plot: accuracy and RT
% RT Distribution (with Mean Line) ---
figure('Name', 'RT Distribution');
histogram(T_correct.rt, 50);
hold on;

mean_rt_overall = mean(T_correct.rt);
xline(mean_rt_overall, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean RT = %.3fs', mean_rt_overall));

title('Distribution of Reaction Times (Correct Trials Only)');
xlabel('Reaction Time (s)');
ylabel('Frequency (Count)');
grid on;
hold off;

% Speed-Accuracy Tradeoff
figure('Name', 'Speed-Accuracy Tradeoff');
Stats_Acc_Overall = groupsummary(T.test, 'subj_id', 'mean', 'correct');
Stats_RT_Overall = groupsummary(T_correct, 'subj_id', 'mean', 'rt');
Stats_SpeedAcc = outerjoin(Stats_Acc_Overall, Stats_RT_Overall, 'Keys', 'subj_id', 'MergeKeys', true);
Stats_SpeedAcc.Properties.VariableNames = {'subj_id', 'GroupCount_Acc', 'mean_correct', 'GroupCount_RT', 'mean_rt'};

plot(Stats_SpeedAcc.mean_rt, Stats_SpeedAcc.mean_correct, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
lsline;
hold off;
title('Speed-Accuracy Tradeoff (Each point = 1 subject)');
xlabel('Mean RT (s) - Correct Trials');
ylabel('Mean Accuracy (All Trials)');
grid on;

% accuracy and RT
[G_acc, T_acc] = findgroups(Stats_Acc.trial_type_final);
mean_acc_by_type = splitapply(@mean, Stats_Acc.mean_correct, G_acc);
sem_acc_by_type = splitapply(@(x) std(x)/sqrt(length(x)), Stats_Acc.mean_correct, G_acc);
T_acc = categorical(T_acc);

% Bar chart
figure('Name', 'Accuracy by Trial Type');
b_acc = bar(mean_acc_by_type);
hold on;
x = 1:length(mean_acc_by_type);
errorbar(x, mean_acc_by_type, sem_acc_by_type, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTick', x, 'XTickLabel', cellstr(T_acc));

title('Mean Accuracy by Trial Type');
ylabel('Mean Accuracy');
xlabel('Trial Type');
ylim([0 1.1]);
grid on;
hold off;


[G_rt, T_rt] = findgroups(Stats_RT.trial_type_final);
mean_rt_by_type = splitapply(@mean, Stats_RT.mean_rt, G_rt);
sem_rt_by_type = splitapply(@(x) std(x)/sqrt(length(x)), Stats_RT.mean_rt, G_rt);
T_rt = categorical(T_rt);

figure('Name', 'RT by Trial Type');
b_rt = bar(mean_rt_by_type); 
hold on;
x = 1:length(mean_rt_by_type);
errorbar(x, mean_rt_by_type, sem_rt_by_type, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTick', x, 'XTickLabel', cellstr(T_rt));

title('Mean RT by Trial Type (Correct Trials Only)');
ylabel('Mean Reaction Time (s)');
xlabel('Trial Type');
grid on;
hold off;





%% 8. LURE ANALYSIS BY CONDITION
% This section specifically analyzes performance on "lure" trials,
% broken down by the encoding "condition" ("isolated_both", "comparison", "novel").

fprintf('Running analysis on Lure trials by condition...\n');

% --- 1. Define the desired order of conditions ---
% (Based on your prompt: "isolated_both", "comparison", and "novel")
condition_order = ["isolated_both", "comparison", "novel"];

% --- 2. Filter data for ONLY lure trials ---
% T.test has all trials, for accuracy analysis
T_lures = T.test(T.test.trial_type_final == "lure", :);
% T_correct has only correct+responded trials, for RT analysis
T_lures_correct = T_correct(T_correct.trial_type_final == "lure", :);

% --- 3. Convert the 'condition' column to categorical with specific order ---
% This ensures the plots are in your desired order, not alphabetical.
% We do this on our filtered tables so we don't alter T.test.
T_lures.condition = categorical(T_lures.condition, condition_order, 'Ordinal', false);
T_lures_correct.condition = categorical(T_lures_correct.condition, condition_order, 'Ordinal', false);

% --- 4. Get per-subject stats for lures ---
% Accuracy for lures, by subject and condition
Stats_Lure_Acc = groupsummary(T_lures, {'subj_id', 'condition'}, 'mean', 'correct');
% RT for correct lures, by subject and condition
Stats_Lure_RT = groupsummary(T_lures_correct, {'subj_id', 'condition'}, 'mean', 'rt');

% --- 5. Get across-subject stats (Mean & SEM) for plotting ---
% Accuracy
[G_lure_acc, T_lure_acc] = findgroups(Stats_Lure_Acc.condition);
mean_lure_acc = splitapply(@mean, Stats_Lure_Acc.mean_correct, G_lure_acc);
sem_lure_acc = splitapply(@(x) std(x)/sqrt(length(x)), Stats_Lure_Acc.mean_correct, G_lure_acc);

% RT
[G_lure_rt, T_lure_rt] = findgroups(Stats_Lure_RT.condition);
mean_lure_rt = splitapply(@mean, Stats_Lure_RT.mean_rt, G_lure_rt);
sem_lure_rt = splitapply(@(x) std(x)/sqrt(length(x)), Stats_Lure_RT.mean_rt, G_lure_rt);

% --- 6. Plot Lure Accuracy by Condition ---
figure('Name', 'Lure Accuracy by Condition');
b_acc = bar(mean_lure_acc);
hold on;
x_acc = 1:length(mean_lure_acc);
errorbar(x_acc, mean_lure_acc, sem_lure_acc, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTick', x_acc, 'XTickLabel', cellstr(T_lure_acc));
title('Lure Identification Accuracy by Condition');
% A more descriptive Y-label for this specific analysis
ylabel('Mean Accuracy (p(Correct "Similar" Response | Lure))');
xlabel('Encoding Condition');
ylim([0 1.1]); % Accuracy plot
grid on;
hold off;

% --- 7. Plot Lure RT by Condition ---
figure('Name', 'Lure RT by Condition');
b_rt = bar(mean_lure_rt);
hold on;
x_rt = 1:length(mean_lure_rt);
errorbar(x_rt, mean_lure_rt, sem_lure_rt, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTick', x_rt, 'XTickLabel', cellstr(T_lure_rt));
title('Lure Identification RT by Condition (Correct Trials)');
ylabel('Mean Reaction Time (s)');
xlabel('Encoding Condition');
% Let Y-axis be automatic for RT, but add grid
grid on; 
hold off;

% ------------------------------------------------------------------------



%% 9. ENCODING: LURE ANALYSIS BY CONDITION
encStrVars = {'response_key', 'correct_response', 'condition'};
for v = encStrVars
    if ismember(v{1}, T.encod.Properties.VariableNames)
        if ~isstring(T.encod.(v{1}))
            T.encod.(v{1}) = string(T.encod.(v{1}));
        end
    else
        % It's ok if 'condition' isn't there, but we need the response keys
        if v{1} ~= "condition"
            warning('Encoding table T.encod is missing expected variable: %s', v{1});
        end
    end
end


% Define the logical conditions for a correct response (will be true/false)
cond1 = (T.encod.correct_response == "none" & T.encod.response_key == "NA");
cond2 = (T.encod.correct_response == "j" & T.encod.response_key == "j");
cond3 = (T.encod.correct_response == "k" & T.encod.response_key == "k");

% Combine all "correct" conditions with an element-wise OR (|)
is_correct_logical = cond1 | cond2 | cond3;

% Convert the logical (true/false) to double (1/0) for calculations
T.encod.correct = double(is_correct_logical);

