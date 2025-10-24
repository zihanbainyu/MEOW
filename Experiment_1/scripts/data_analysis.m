clear; clc;

%% --- Configuration ---
data_dir = fullfile('..','data');
participant_ids = [901,902,903,904];
KEY_SIMILAR = 'd'; % for lure rejection

%% --- Load and combine data ---
all_results_table = table();

for p_idx = 1:length(participant_ids)
    subj_id = participant_ids(p_idx);
    filename = fullfile(data_dir, sprintf('sub%03d',subj_id), sprintf('sub%03d_concat.mat',subj_id));
    
    if ~exist(filename,'file')
        warning('File not found: %s. Skipping...', filename);
        continue;
    end
    
    load(filename,'final_data_output');
    
    if ~isfield(final_data_output,'test_results')
        warning('No test_results field for sub%03d. Skipping...', subj_id);
        continue;
    end
    
    tbl = final_data_output.test_results;
    tbl.subj_id = repmat(subj_id,height(tbl),1);
    
    % Lure rejection correctness
    lure_idx = (tbl.correct_response == KEY_SIMILAR);
    tbl.lure_rejection_correct = NaN(height(tbl),1);
    tbl.lure_rejection_correct(lure_idx) = (tbl.response_key(lure_idx) == KEY_SIMILAR);
    
    % Correct RT only
    tbl.rt_correct = tbl.rt;
    tbl.rt_correct(~tbl.correct) = NaN;  % <-- use your "correct" column
    
    all_results_table = [all_results_table; tbl];
end

%% --- Prepare data for plotting ---
if ~ismember('trial_type_final', all_results_table.Properties.VariableNames)
    error('Column "trial_type_final" not found in data table.');
end

all_results_table.lure_type = categorical(all_results_table.trial_type_final);
valid_lure_types = {'ComparativeLure','IsolatedLure','NovelLure'};
analysis_table = all_results_table(ismember(all_results_table.lure_type, valid_lure_types), :);

if isempty(analysis_table)
    error('No valid lure trials found. Check your trial type column or data.');
end

%% --- Aggregate mean and SEM per lure type ---
mean_acc = grpstats(analysis_table,'lure_type',{'mean','sem'},'DataVars','lure_rejection_correct');
mean_rt = grpstats(analysis_table,'lure_type',{'mean','sem'},'DataVars','rt_correct');

%% --- Plot Lure Rejection Accuracy ---
figure;
bar(mean_acc.mean_lure_rejection_correct); hold on;
errorbar(1:height(mean_acc), mean_acc.mean_lure_rejection_correct, mean_acc.sem_lure_rejection_correct,'k','linestyle','none');
set(gca,'XTickLabel',cellstr(mean_acc.lure_type));
ylabel('Lure Rejection Accuracy');
ylim([0 1]);
title('Preliminary Lure Rejection Accuracy');

%% --- Plot Reaction Time ---
figure;
bar(mean_rt.mean_rt_correct); hold on;
errorbar(1:height(mean_rt), mean_rt.mean_rt_correct, mean_rt.sem_rt_correct,'k','linestyle','none');
set(gca,'XTickLabel',cellstr(mean_rt.lure_type));
ylabel('Reaction Time (ms)');
title('Preliminary Reaction Time (Correct Lures)');
