clear;
clc;
close all;

%% 1. Group Setup
% -------------------------------------------------------------------------
subj_ids = [602, 601]; % <--- ADD YOUR SUBJECT LIST HERE
% -------------------------------------------------------------------------

base_dir = '..'; 
MIN_RT_CUTOFF = 0.150;

% Define standard colors and fonts
plot_font_name = 'Helvetica';
plot_font_size_axis = 12;
plot_font_size_title = 14;
plot_font_size_sgtitle = 16;
plot_color_comp = [180 174 211]/255; 
plot_color_iso  = [176 230 255]/255;
plot_color_nov  = [183 210 205]/255; 
plot_color_similar = [255 191 205]/255; 
plot_color_same = [97 125 184]/255; 
plot_color_new = [219 219 219]/255;

% Initialize a structure array to hold every subject's stats
all_subjs = struct(); 

%% 2. Analysis Loop
for s = 1:length(subj_ids)
    curr_id = subj_ids(s);
    fprintf('Processing Subject %d...\n', curr_id);
    
    % --- Load Data ---
    subj_folder = sprintf('sub%03d', curr_id); 
    results_dir = fullfile(base_dir, 'data', subj_folder);
    concat_file = fullfile(results_dir, sprintf('example.mat'));
    % concat_file = fullfile(results_dir, sprintf('sub%03d_concat.mat', curr_id));
    rec_file = fullfile(results_dir, sprintf('sub%03d_rec.mat', curr_id)); 
    
    if ~exist(concat_file, 'file') || ~exist(rec_file, 'file')
        warning('Files not found for subject %d. Skipping.', curr_id);
        continue;
    end
    
    load(concat_file, 'final_data_output');
    load(rec_file, 'results_recognition');
    
    results_1_back = final_data_output.results_1_back_all;
    results_2_back = final_data_output.results_2_back_all;
    
    % Initialize temp stats for this subject
    subj_stats = [];

    % =====================================================================
    % 1-Back Math
    % =====================================================================
    results_1_back.corr_resp = cellstr(results_1_back.corr_resp);
    results_1_back.resp_key = cellstr(results_1_back.resp_key);
    na_idx = strcmp(results_1_back.resp_key, 'NA');
    results_1_back.resp_key(na_idx) = {'none'};
    results_1_back.correct = strcmp(results_1_back.corr_resp, results_1_back.resp_key);
    
    % Indices
    comp_B_idx = (results_1_back.condition == "compared" & results_1_back.identity == "B");
    rep_A2_idx = (results_1_back.condition == "repeat" & strcmp(results_1_back.corr_resp, 'j'));
    % Correct Rejection Indices
    comp_A_idx = (results_1_back.condition == "compared" & results_1_back.identity == "A");
    iso_A_idx = (results_1_back.condition == "isolated" & results_1_back.identity == "A");
    iso_B_idx = (results_1_back.condition == "isolated" & results_1_back.identity == "B");
    rep_A1_idx = (results_1_back.condition == "repeat" & strcmp(results_1_back.corr_resp, 'none'));
    all_cr_trials_idx = comp_A_idx | iso_A_idx | iso_B_idx | rep_A1_idx;
    
    % Accuracy
    subj_stats.oneback.acc_hit_similar = mean(results_1_back.correct(comp_B_idx));
    subj_stats.oneback.acc_hit_same = mean(results_1_back.correct(rep_A2_idx));
    subj_stats.oneback.acc_cr_new = mean(results_1_back.correct(all_cr_trials_idx));
    
    % RT
    valid_rt_idx = results_1_back.rt > MIN_RT_CUTOFF;
    subj_stats.oneback.rt_hit_similar = median(results_1_back.rt(comp_B_idx & results_1_back.correct == 1 & valid_rt_idx), 'omitnan');
    subj_stats.oneback.rt_hit_same = median(results_1_back.rt(rep_A2_idx & results_1_back.correct == 1 & valid_rt_idx), 'omitnan');
    
    % Errors (Proportions)
    n_comp_B = sum(comp_B_idx);
    subj_stats.oneback.err_comp_as_same = sum(strcmp(results_1_back.resp_key(comp_B_idx), 'j')) / n_comp_B;
    subj_stats.oneback.err_comp_as_new = sum(strcmp(results_1_back.resp_key(comp_B_idx), 'none')) / n_comp_B;
    
    n_rep_A2 = sum(rep_A2_idx);
    subj_stats.oneback.err_rep_as_similar = sum(strcmp(results_1_back.resp_key(rep_A2_idx), 'k')) / n_rep_A2;
    subj_stats.oneback.err_rep_as_new = sum(strcmp(results_1_back.resp_key(rep_A2_idx), 'none')) / n_rep_A2;
    
    n_cr_new = sum(all_cr_trials_idx);
    subj_stats.oneback.err_new_as_same = sum(strcmp(results_1_back.resp_key(all_cr_trials_idx), 'j')) / n_cr_new;
    subj_stats.oneback.err_new_as_similar = sum(strcmp(results_1_back.resp_key(all_cr_trials_idx), 'k')) / n_cr_new;

    % =====================================================================
    % 2-Back Math
    % =====================================================================
    results_2_back.corr_resp = cellstr(results_2_back.corr_resp);
    results_2_back.resp_key = cellstr(results_2_back.resp_key);
    na_idx = strcmp(results_2_back.resp_key, 'NA');
    results_2_back.resp_key(na_idx) = {'none'};
    results_2_back.correct = strcmp(results_2_back.corr_resp, results_2_back.resp_key);
    
    % Post A-N logic
    results_2_back.post_AN = zeros(height(results_2_back), 1);
    for i = 1:(height(results_2_back) - 2)
        if strcmp(results_2_back.goal(i), 'A-N')
            results_2_back.post_AN(i + 2) = 1;
        end
    end

    % Indices
    real_trials_idx = ~contains(results_2_back.goal, "JUNK");
    valid_rt_idx = results_2_back.rt > MIN_RT_CUTOFF;
    aa_idx = real_trials_idx & strcmp(results_2_back.goal, 'A-A');
    ab_idx = real_trials_idx & strcmp(results_2_back.goal, 'A-B'); 
    an_idx = (results_2_back.post_AN == 1);
    comp_idx = real_trials_idx & strcmp(results_2_back.condition, 'compared');
    iso_idx = real_trials_idx & strcmp(results_2_back.condition, 'isolated');
    nov_idx = real_trials_idx & strcmp(results_2_back.condition, 'novel');
    j_idx = real_trials_idx & strcmp(results_2_back.corr_resp, 'j'); 
    k_idx = real_trials_idx & strcmp(results_2_back.corr_resp, 'k'); 
    
    % Accuracy
    subj_stats.twoback.acc_AA_comp = mean(results_2_back.correct(aa_idx & comp_idx & j_idx));
    subj_stats.twoback.acc_AA_iso = mean(results_2_back.correct(aa_idx & iso_idx & j_idx));
    subj_stats.twoback.acc_AA_nov = mean(results_2_back.correct(aa_idx & nov_idx & j_idx));
    
    subj_stats.twoback.acc_AB_comp = mean(results_2_back.correct(ab_idx & comp_idx & k_idx));
    subj_stats.twoback.acc_AB_iso = mean(results_2_back.correct(ab_idx & iso_idx & k_idx));
    subj_stats.twoback.acc_AB_nov = mean(results_2_back.correct(ab_idx & nov_idx & k_idx));
    
    subj_stats.twoback.acc_AN_comp = mean(results_2_back.correct(an_idx & comp_idx));
    subj_stats.twoback.acc_AN_iso = mean(results_2_back.correct(an_idx & iso_idx));
    subj_stats.twoback.acc_AN_nov = mean(results_2_back.correct(an_idx & nov_idx));
    
    % RT
    subj_stats.twoback.rt_AA_comp = median(results_2_back.rt(aa_idx & comp_idx & j_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
    subj_stats.twoback.rt_AA_iso = median(results_2_back.rt(aa_idx & iso_idx & j_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
    subj_stats.twoback.rt_AA_nov = median(results_2_back.rt(aa_idx & nov_idx & j_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
    
    subj_stats.twoback.rt_AB_comp = median(results_2_back.rt(ab_idx & comp_idx & k_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
    subj_stats.twoback.rt_AB_iso = median(results_2_back.rt(ab_idx & iso_idx & k_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
    subj_stats.twoback.rt_AB_nov = median(results_2_back.rt(ab_idx & nov_idx & k_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');

    % LDI (Lure Discrimination Index)
    subj_stats.twoback.err_AN_comp_as_k = sum(strcmp(results_2_back.resp_key(an_idx & comp_idx), 'k')) / sum(an_idx & comp_idx);
    subj_stats.twoback.err_AN_iso_as_k = sum(strcmp(results_2_back.resp_key(an_idx & iso_idx), 'k')) / sum(an_idx & iso_idx);
    subj_stats.twoback.err_AN_nov_as_k = sum(strcmp(results_2_back.resp_key(an_idx & nov_idx), 'k')) / sum(an_idx & nov_idx);
    
    subj_stats.twoback.ldi_comp = subj_stats.twoback.acc_AB_comp - subj_stats.twoback.err_AN_comp_as_k;
    subj_stats.twoback.ldi_iso = subj_stats.twoback.acc_AB_iso - subj_stats.twoback.err_AN_iso_as_k;
    subj_stats.twoback.ldi_nov = subj_stats.twoback.acc_AB_nov - subj_stats.twoback.err_AN_nov_as_k;

    % d-prime (Identity Recognition)
    % Note: Simplified d' calc for brevity, assumes standard formula as per original code
    n_AA_comp = sum(aa_idx & comp_idx & j_idx);
    n_AN_comp = sum(an_idx & comp_idx);
    hr = subj_stats.twoback.acc_AA_comp; far = sum(strcmp(results_2_back.resp_key(an_idx & comp_idx), 'j')) / n_AN_comp;
    if hr==1, hr=1-(1/(2*n_AA_comp)); elseif hr==0, hr=1/(2*n_AA_comp); end
    if far==0, far=1/(2*n_AN_comp); elseif far==1, far=1-(1/(2*n_AN_comp)); end
    subj_stats.twoback.d_AA_comp = norminv(hr) - norminv(far);

    n_AA_iso = sum(aa_idx & iso_idx & j_idx);
    n_AN_iso = sum(an_idx & iso_idx);
    hr = subj_stats.twoback.acc_AA_iso; far = sum(strcmp(results_2_back.resp_key(an_idx & iso_idx), 'j')) / n_AN_iso;
    if hr==1, hr=1-(1/(2*n_AA_iso)); elseif hr==0, hr=1/(2*n_AA_iso); end
    if far==0, far=1/(2*n_AN_iso); elseif far==1, far=1-(1/(2*n_AN_iso)); end
    subj_stats.twoback.d_AA_iso = norminv(hr) - norminv(far);

    n_AA_nov = sum(aa_idx & nov_idx & j_idx);
    n_AN_nov = sum(an_idx & nov_idx);
    hr = subj_stats.twoback.acc_AA_nov; far = sum(strcmp(results_2_back.resp_key(an_idx & nov_idx), 'j')) / n_AN_nov;
    if hr==1, hr=1-(1/(2*n_AA_nov)); elseif hr==0, hr=1/(2*n_AA_nov); end
    if far==0, far=1/(2*n_AN_nov); elseif far==1, far=1-(1/(2*n_AN_nov)); end
    subj_stats.twoback.d_AA_nov = norminv(hr) - norminv(far);

    % =====================================================================
    % Recognition Math
    % =====================================================================
    results_recognition.corr_resp = cellstr(results_recognition.corr_resp);
    results_recognition.resp_key = cellstr(results_recognition.resp_key);
    results_recognition.correct = strcmp(results_recognition.corr_resp, results_recognition.resp_key);
    
    old_trials = results_recognition(results_recognition.trial_type == "old", :);
    new_trials = results_recognition(results_recognition.trial_type ~= "old", :);
    valid_rt_idx = results_recognition.rt > MIN_RT_CUTOFF;

    % Overall FAR
    n_new = height(new_trials);
    n_fa = sum(strcmp(new_trials.resp_key, 'j') & new_trials.rt > MIN_RT_CUTOFF);
    FAR = n_fa / n_new;
    if FAR == 0, FAR = 1/(2*n_new); elseif FAR == 1, FAR = 1-(1/(2*n_new)); end
    subj_stats.recog.FAR_overall = FAR;
    
    % Comp
    comp_trials = old_trials(old_trials.condition == "compared", :);
    n_comp = height(comp_trials);
    n_hits_comp = sum(comp_trials.correct == 1 & comp_trials.rt > MIN_RT_CUTOFF);
    HR_comp = n_hits_comp / n_comp;
    if HR_comp==1, HR_comp=1-(1/(2*n_comp)); elseif HR_comp==0, HR_comp=1/(2*n_comp); end
    subj_stats.recog.HR_comp = HR_comp;
    subj_stats.recog.d_prime_comp = norminv(HR_comp) - norminv(FAR);
    subj_stats.recog.rt_comp_hit = mean(comp_trials.rt(comp_trials.correct==1 & comp_trials.rt>MIN_RT_CUTOFF), 'omitnan');

    % Iso
    iso_trials = old_trials(old_trials.condition == "isolated", :);
    n_iso = height(iso_trials);
    n_hits_iso = sum(iso_trials.correct == 1 & iso_trials.rt > MIN_RT_CUTOFF);
    HR_iso = n_hits_iso / n_iso;
    if HR_iso==1, HR_iso=1-(1/(2*n_iso)); elseif HR_iso==0, HR_iso=1/(2*n_iso); end
    subj_stats.recog.HR_iso = HR_iso;
    subj_stats.recog.d_prime_iso = norminv(HR_iso) - norminv(FAR);
    subj_stats.recog.rt_iso_hit = mean(iso_trials.rt(iso_trials.correct==1 & iso_trials.rt>MIN_RT_CUTOFF), 'omitnan');

    % A vs B
    A_trials = old_trials(old_trials.identity == "A", :);
    n_A = height(A_trials);
    HR_A = sum(A_trials.correct==1 & A_trials.rt>MIN_RT_CUTOFF)/n_A;
    if HR_A==1, HR_A=1-(1/(2*n_A)); elseif HR_A==0, HR_A=1/(2*n_A); end
    subj_stats.recog.d_prime_A = norminv(HR_A) - norminv(FAR);
    subj_stats.recog.rt_A_hit = mean(A_trials.rt(A_trials.correct==1 & A_trials.rt>MIN_RT_CUTOFF), 'omitnan');

    B_trials = old_trials(old_trials.identity == "B", :);
    n_B = height(B_trials);
    HR_B = sum(B_trials.correct==1 & B_trials.rt>MIN_RT_CUTOFF)/n_B;
    if HR_B==1, HR_B=1-(1/(2*n_B)); elseif HR_B==0, HR_B=1/(2*n_B); end
    subj_stats.recog.d_prime_B = norminv(HR_B) - norminv(FAR);
    subj_stats.recog.rt_B_hit = mean(B_trials.rt(B_trials.correct==1 & B_trials.rt>MIN_RT_CUTOFF), 'omitnan');
    
    subj_stats.recog.rt_new_cr = mean(new_trials.rt(new_trials.correct==1 & new_trials.rt>MIN_RT_CUTOFF), 'omitnan');

    % --- STORE IN MASTER STRUCT ---
    all_subjs(s).stats = subj_stats;
end

%% 3. Helper Functions for Aggregation
% We create anonymous functions to extract vector of data from the struct array
% usage: get_vec('oneback', 'acc_hit_same') returns [subj1_acc, subj2_acc...]
get_vec = @(field1, field2) arrayfun(@(x) x.stats.(field1).(field2), all_subjs);
calc_sem = @(x) std(x, 'omitnan') / sqrt(sum(~isnan(x)));

%% 4. Plots (Using Group Means & SEM)

% --- 1-Back Figure ---
figure('color', 'white', 'Position', [100 100 1200 400]); 

% 1. Accuracy
h_ax1 = subplot(1, 3, 1);
data_mat = [get_vec('oneback', 'acc_hit_same'); get_vec('oneback', 'acc_hit_similar'); get_vec('oneback', 'acc_cr_new')]';
M = mean(data_mat, 1);
E = std(data_mat, 0, 1) ./ sqrt(size(data_mat, 1));

b1 = bar(M, 'FaceColor', 'flat');
hold on;
errorbar(1:3, M, E, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
b1.CData(1,:) = plot_color_same; b1.CData(2,:) = plot_color_similar; b1.CData(3,:) = plot_color_new;
ylabel('Group Accuracy', 'FontSize', plot_font_size_axis);
set(h_ax1, 'XTickLabel', {'Hit(Same)', 'Hit(Sim)', 'CR(New)'}, 'FontSize', plot_font_size_axis);
xtickangle(45); ylim([0, 1.05]);

% 2. RT
h_ax2 = subplot(1, 3, 2);
rt_mat = [get_vec('oneback', 'rt_hit_same'); get_vec('oneback', 'rt_hit_similar')]';
M = mean(rt_mat, 1, 'omitnan');
E = std(rt_mat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(rt_mat)));

b2 = bar(M, 'FaceColor', 'flat');
hold on;
errorbar(1:2, M, E, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
b2.CData(1,:) = plot_color_same; b2.CData(2,:) = plot_color_similar;
ylabel('Group RT (s)', 'FontSize', plot_font_size_axis);
set(h_ax2, 'XTickLabel', {'Hit(Same)', 'Hit(Sim)'}, 'FontSize', plot_font_size_axis);
ylim([0, max(M)*1.3]);

% 3. Errors (Stacked) - Note: Error bars are messy on stacked, plotting Means only
h_ax3 = subplot(1, 3, 3);
% Get means for every component
m_rep_sim = mean(get_vec('oneback', 'err_rep_as_similar')); m_rep_new = mean(get_vec('oneback', 'err_rep_as_new'));
m_comp_same = mean(get_vec('oneback', 'err_comp_as_same')); m_comp_new = mean(get_vec('oneback', 'err_comp_as_new'));
m_new_same = mean(get_vec('oneback', 'err_new_as_same')); m_new_sim = mean(get_vec('oneback', 'err_new_as_similar'));

% Use the calculated group means for the stack
% Col 1: Exp Same (Correct, ErrSim, ErrNew)
% Col 2: Exp Sim  (ErrSame, Correct, ErrNew)
% Col 3: Exp New  (ErrSame, ErrSim, Correct)
Y_stack = [
    mean(get_vec('oneback','acc_hit_same')), m_rep_sim, m_rep_new;
    m_comp_same, mean(get_vec('oneback','acc_hit_similar')), m_comp_new;
    m_new_same, m_new_sim, mean(get_vec('oneback','acc_cr_new'))
];
% Reorder for plotting (Bottom, Middle, Top segments)
plot_data = [Y_stack(1,1), Y_stack(1,2), Y_stack(1,3); 
             Y_stack(2,2), Y_stack(2,1), Y_stack(2,3); 
             Y_stack(3,3), Y_stack(3,1), Y_stack(3,2)];
         
b_stack = bar(plot_data, 'stacked', 'FaceColor', 'flat');
% Apply colors (Same logic as original: diagonal logic)
b_stack(1).CData(1,:) = plot_color_same; b_stack(1).CData(2,:) = plot_color_similar; b_stack(1).CData(3,:) = plot_color_new;
b_stack(2).CData(1,:) = plot_color_similar; b_stack(2).CData(2,:) = plot_color_same; b_stack(2).CData(3,:) = plot_color_same;
b_stack(3).CData(1,:) = plot_color_new; b_stack(3).CData(2,:) = plot_color_new; b_stack(3).CData(3,:) = plot_color_similar;

set(h_ax3, 'XTickLabel', {'Exp Same', 'Exp Sim', 'Exp New'}, 'FontSize', plot_font_size_axis);
ylabel('Mean Response Prop', 'FontSize', plot_font_size_axis);
xtickangle(45); sgtitle('Group 1-Back');

% --- 2-Back LDI & D-Prime ---
figure('color', 'white', 'Position', [100 100 800 400]);

% LDI
h_ax1 = subplot(1, 2, 1);
ldi_mat = [get_vec('twoback', 'ldi_comp'); get_vec('twoback', 'ldi_iso'); get_vec('twoback', 'ldi_nov')]';
M = mean(ldi_mat, 1); E = std(ldi_mat,0,1)./sqrt(size(ldi_mat,1));

b1 = bar(M, 'FaceColor', 'flat'); hold on;
errorbar(1:3, M, E, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
b1.CData(1,:) = plot_color_comp; b1.CData(2,:) = plot_color_iso; b1.CData(3,:) = plot_color_nov;
ylabel('Group LDI'); title('Lure Discrimination');
set(h_ax1, 'XTickLabel', {'Comp', 'Iso', 'Nov'});

% D-Prime
h_ax2 = subplot(1, 2, 2);
d_mat = [get_vec('twoback', 'd_AA_comp'); get_vec('twoback', 'd_AA_iso'); get_vec('twoback', 'd_AA_nov')]';
M = mean(d_mat, 1); E = std(d_mat,0,1)./sqrt(size(d_mat,1));

b2 = bar(M, 'FaceColor', 'flat'); hold on;
errorbar(1:3, M, E, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
b2.CData(1,:) = plot_color_comp; b2.CData(2,:) = plot_color_iso; b2.CData(3,:) = plot_color_nov;
ylabel('Group d'''); title('Identity Recognition');
set(h_ax2, 'XTickLabel', {'Comp', 'Iso', 'Nov'});
sgtitle('Group 2-Back');

% --- Recognition Figures ---
figure('color', 'white', 'Position', [100 100 1000 400]);

% 1. d' by Condition
subplot(1,3,1);
d_rec_mat = [get_vec('recog', 'd_prime_comp'); get_vec('recog', 'd_prime_iso')]';
M = mean(d_rec_mat,1); E = std(d_rec_mat,0,1)./sqrt(size(d_rec_mat,1));
b = bar(M, 'FaceColor', 'flat'); hold on;
errorbar(1:2, M, E, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
b.CData(1,:) = plot_color_comp; b.CData(2,:) = plot_color_iso;
title('Group d'' (Condition)'); set(gca, 'XTickLabel', {'Comp', 'Iso'});

% 2. d' by Identity
subplot(1,3,2);
d_id_mat = [get_vec('recog', 'd_prime_A'); get_vec('recog', 'd_prime_B')]';
M = mean(d_id_mat,1); E = std(d_id_mat,0,1)./sqrt(size(d_id_mat,1));
b = bar(M, 'FaceColor', 'flat'); hold on;
errorbar(1:2, M, E, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
b.CData(1,:) = plot_color_same; b.CData(2,:) = plot_color_similar;
title('Group d'' (Identity)'); set(gca, 'XTickLabel', {'A', 'B'});

% 3. HR vs FAR (The one you asked to fix earlier)
subplot(1,3,3);
hold on;
% Calculate means
hr_comp_m = mean(get_vec('recog', 'HR_comp'));
hr_comp_e = calc_sem(get_vec('recog', 'HR_comp'));
hr_iso_m = mean(get_vec('recog', 'HR_iso'));
hr_iso_e = calc_sem(get_vec('recog', 'HR_iso'));
far_m = mean(get_vec('recog', 'FAR_overall')); % Note: No error bar usually plotted for baseline line, but calculated here

% Plot Bars
b1 = bar(1, hr_comp_m); set(b1, 'FaceColor', plot_color_comp);
errorbar(1, hr_comp_m, hr_comp_e, 'k', 'LineWidth', 1.5);

b2 = bar(2, hr_iso_m); set(b2, 'FaceColor', plot_color_iso);
errorbar(2, hr_iso_m, hr_iso_e, 'k', 'LineWidth', 1.5);

% Plot Line
l1 = yline(far_m, 'Color', plot_color_new, 'LineStyle', '--', 'LineWidth', 2);

title('Group HR vs FAR'); ylabel('Proportion');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Comp', 'Iso'});
legend([b1, b2, l1], {'Compared', 'Isolated', 'FAR'}, 'Location', 'southoutside');
ylim([0, 1.0]);

sgtitle('Group Recognition');

disp('Group Analysis Complete.');