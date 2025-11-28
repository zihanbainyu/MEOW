
clear;
clc;
close all;


%% Set up
subj_id = 607;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id); 
results_dir = fullfile(base_dir, 'data', subj_folder);
concat_file = fullfile(results_dir, sprintf('sub%03d_concat.mat', subj_id));
% concat_file = fullfile(results_dir, sprintf('example.mat', subj_id));
rec_file = fullfile(results_dir, sprintf('sub%03d_rec.mat', subj_id)); 

load(concat_file, 'final_data_output');
load(rec_file, 'results_recognition');

MIN_RT_CUTOFF = 0.150;

subj_stats = [];
results_1_back = final_data_output.results_1_back_all;
results_2_back = final_data_output.results_2_back_all;

% define standard colors and fonts
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



%% Math. 1-Back

results_1_back.corr_resp = cellstr(results_1_back.corr_resp);
results_1_back.resp_key = cellstr(results_1_back.resp_key);
na_idx = strcmp(results_1_back.resp_key, 'NA');
results_1_back.resp_key(na_idx) = {'none'};
results_1_back.correct = strcmp(results_1_back.corr_resp, results_1_back.resp_key);

% overall accuracy
subj_stats.oneback.overall_acc = mean(results_1_back.correct);

% 'compared' trials
comp_A_idx = (results_1_back.condition == "compared" & results_1_back.identity == "A");
comp_B_idx = (results_1_back.condition == "compared" & results_1_back.identity == "B");

% 'isolated' trials
iso_A_idx = (results_1_back.condition == "isolated" & results_1_back.identity == "A");
iso_B_idx = (results_1_back.condition == "isolated" & results_1_back.identity == "B");

% 'repeat' trials
rep_A1_idx = (results_1_back.condition == "repeat" & strcmp(results_1_back.corr_resp, 'none'));
rep_A2_idx = (results_1_back.condition == "repeat" & strcmp(results_1_back.corr_resp, 'j'));

subj_stats.oneback.acc_hit_similar = mean(results_1_back.correct(comp_B_idx));

% accuracy: "same" hit (correctly press 'j' for 2nd repeat 'A')
subj_stats.oneback.acc_hit_same = mean(results_1_back.correct(rep_A2_idx));

% accuracy: "correct rejection" trials (where 'none' is required)
all_cr_trials_idx = comp_A_idx | iso_A_idx | iso_B_idx | rep_A1_idx;
subj_stats.oneback.acc_cr_new = mean(results_1_back.correct(all_cr_trials_idx));

subj_stats.oneback.fa_new = 1 - subj_stats.oneback.acc_cr_new;

% RT
valid_rt_idx = results_1_back.rt > MIN_RT_CUTOFF;

% median rt for correct "similar" hits
subj_stats.oneback.rt_hit_similar = median(results_1_back.rt(...
    comp_B_idx & ...
    results_1_back.correct == 1 & ...
    valid_rt_idx), 'omitnan');

% median rt for correct "same" hits
subj_stats.oneback.rt_hit_same = median(results_1_back.rt(...
    rep_A2_idx & ...
    results_1_back.correct == 1 & ...
    valid_rt_idx), 'omitnan');

% 'compared B' trials (correct = 'k')
n_comp_B = sum(comp_B_idx);
subj_stats.oneback.err_comp_as_same = sum(strcmp(results_1_back.resp_key(comp_B_idx), 'j')) / n_comp_B;
subj_stats.oneback.err_comp_as_new = sum(strcmp(results_1_back.resp_key(comp_B_idx), 'none')) / n_comp_B;

% 'repeat A2' trials (correct = 'j')
n_rep_A2 = sum(rep_A2_idx);
subj_stats.oneback.err_rep_as_similar = sum(strcmp(results_1_back.resp_key(rep_A2_idx), 'k')) / n_rep_A2;
subj_stats.oneback.err_rep_as_new = sum(strcmp(results_1_back.resp_key(rep_A2_idx), 'none')) / n_rep_A2;

% 'new' trials (correct = 'none')
n_cr_new = sum(all_cr_trials_idx);
subj_stats.oneback.err_new_as_same = sum(strcmp(results_1_back.resp_key(all_cr_trials_idx), 'j')) / n_cr_new;
subj_stats.oneback.err_new_as_similar = sum(strcmp(results_1_back.resp_key(all_cr_trials_idx), 'k')) / n_cr_new;

% print
fprintf('overall accuracy: %f\n', subj_stats.oneback.overall_acc);
fprintf('accuracy hit-same: %f\n', subj_stats.oneback.acc_hit_same);
fprintf('accuracy hit-similar: %f\n', subj_stats.oneback.acc_hit_similar);
fprintf('accuracy cr: %f\n', subj_stats.oneback.acc_cr_new);

%% Plots. 1-Back
figure('color', 'white'); 

% 1-Back accuracy
h_ax1 = subplot(1, 3, 1);
acc_data = [subj_stats.oneback.acc_hit_same, subj_stats.oneback.acc_hit_similar, subj_stats.oneback.acc_cr_new];
b1 = bar(acc_data, 'FaceColor', 'flat');
b1.CData(1,:) = plot_color_same;
b1.CData(2,:) = plot_color_similar;
b1.CData(3,:) = plot_color_new;
grid off;
box off;
ylabel('Accuracy', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
set(h_ax1, 'XTickLabel', {'Hit (Same)', 'Hit (Similar)', 'CR (New)'}, 'FontName', 'Helvetica', 'FontSize', plot_font_size_axis);
xtickangle(h_ax1, 45);
ylim([0, 1.05]);

% 1-Back RT
h_ax2 = subplot(1, 3, 2);
rt_data = [subj_stats.oneback.rt_hit_same, subj_stats.oneback.rt_hit_similar];
b2 = bar(rt_data, 'FaceColor', 'flat');
b2.CData(1,:) = plot_color_same;
b2.CData(2,:) = plot_color_similar;
grid off;
box off;
ylabel('RT (s)', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
set(h_ax2, 'XTickLabel', {'Hit (Same)', 'Hit (Similar)'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
xtickangle(h_ax2, 45);
ylim([0, 1.05]);
if all(~isnan(rt_data)), ylim(h_ax2, [0, max(rt_data) * 1.2]); end

% 1-Back response proportions (incorrect trials)
h_ax3 = subplot(1, 3, 3);
Y_full = [
    subj_stats.oneback.acc_hit_same, subj_stats.oneback.err_rep_as_similar, subj_stats.oneback.err_rep_as_new; ...
    subj_stats.oneback.err_comp_as_same, subj_stats.oneback.acc_hit_similar, subj_stats.oneback.err_comp_as_new; ...
    subj_stats.oneback.err_new_as_same, subj_stats.oneback.err_new_as_similar, subj_stats.oneback.acc_cr_new];
plot_data = [Y_full(1,1), Y_full(1,2), Y_full(1,3); ...
    Y_full(2,2), Y_full(2,1), Y_full(2,3); ...
    Y_full(3,3), Y_full(3,1), Y_full(3,2)];
b_stack = bar(plot_data, 'stacked', 'FaceColor', 'flat');
b_stack(1).CData(1,:) = plot_color_same; 
b_stack(1).CData(2,:) = plot_color_similar;
b_stack(1).CData(3,:) = plot_color_new;
b_stack(2).CData(1,:) = plot_color_similar;
b_stack(2).CData(2,:) = plot_color_same;
b_stack(2).CData(3,:) = plot_color_same;
b_stack(3).CData(1,:) = plot_color_new;
b_stack(3).CData(2,:) = plot_color_new;
b_stack(3).CData(3,:) = plot_color_similar;
box off;
grid off;
set(h_ax3, 'XTickLabel', {'Expected Same', 'Expected Similar', 'Expected New'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
ylabel('Response Proportions', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
xtickangle(h_ax3, 45);
ylim([0, 1.07]);
hold(h_ax3, 'on');
h_same    = bar(NaN, 'FaceColor', plot_color_same);
h_similar = bar(NaN, 'FaceColor', plot_color_similar);
h_new     = bar(NaN, 'FaceColor', plot_color_new);
hold(h_ax3, 'off');

% Apply one main title to the entire figure
sgtitle('1-Back Task', 'FontSize', plot_font_size_sgtitle, 'FontName', plot_font_name);

% print(gcf, '1Back_Figure_HighRes.tiff', '-dtiff', '-r300'); 
% disp('Saved figure as 1Back_Figure_HighRes.tiff (300 DPI)');












%% Math. 2-Back

results_2_back.corr_resp = cellstr(results_2_back.corr_resp);
results_2_back.resp_key = cellstr(results_2_back.resp_key);
na_idx = strcmp(results_2_back.resp_key, 'NA');
results_2_back.resp_key(na_idx) = {'none'};
results_2_back.correct = strcmp(results_2_back.corr_resp, results_2_back.resp_key);
results_2_back.post_AN = zeros(height(results_2_back), 1);
for i = 1:(height(results_2_back) - 2)
    if strcmp(results_2_back.goal(i), 'A-N')
        results_2_back.post_AN(i + 2) = 1;
    end
end

% overall accuracy
subj_stats.twoback.overall_acc = mean(results_2_back.correct);

% find trials
real_trials_idx = ~contains(results_2_back.goal, "JUNK");
valid_rt_idx = results_2_back.rt > MIN_RT_CUTOFF;
all_valid_rts = results_2_back.rt(valid_rt_idx);

figure('color', 'white');
histogram(all_valid_rts, 50, 'FaceColor', [0.5 0.5 0.5]);
title(sprintf('2-Back: RT Distribution (RTs > %.3fs)', MIN_RT_CUTOFF));
xlabel('RT (s)');
ylabel('Count');
grid on;
box off;

aa_idx = real_trials_idx & strcmp(results_2_back.goal, 'A-A');
ab_idx = real_trials_idx & strcmp(results_2_back.goal, 'A-B'); 
% an_idx = real_trials_idx & strcmp(results_2_back.goal, 'A-N');
j_idx = real_trials_idx & strcmp(results_2_back.corr_resp, 'j'); 
k_idx = real_trials_idx & strcmp(results_2_back.corr_resp, 'k'); 
an_idx = (results_2_back.post_AN == 1);  

comp_idx = real_trials_idx & strcmp(results_2_back.condition, 'compared');
iso_idx = real_trials_idx & strcmp(results_2_back.condition, 'isolated');
nov_idx = real_trials_idx & strcmp(results_2_back.condition, 'novel');

% (A-A) accuracy
subj_stats.twoback.acc_AA_comp = mean(results_2_back.correct(aa_idx & comp_idx & j_idx));
subj_stats.twoback.acc_AA_iso = mean(results_2_back.correct(aa_idx & iso_idx & j_idx));
subj_stats.twoback.acc_AA_nov = mean(results_2_back.correct(aa_idx & nov_idx & j_idx));

% (A-B) accuracy
subj_stats.twoback.acc_AB_comp = mean(results_2_back.correct(ab_idx & comp_idx & k_idx));
subj_stats.twoback.acc_AB_iso = mean(results_2_back.correct(ab_idx & iso_idx & k_idx));
subj_stats.twoback.acc_AB_nov = mean(results_2_back.correct(ab_idx & nov_idx & k_idx));

% (A-N) accuracy
subj_stats.twoback.acc_AN_comp = mean(results_2_back.correct(an_idx & comp_idx));
subj_stats.twoback.acc_AN_iso = mean(results_2_back.correct(an_idx & iso_idx));
subj_stats.twoback.acc_AN_nov = mean(results_2_back.correct(an_idx & nov_idx));

% (A-A) rt
subj_stats.twoback.rt_AA_comp = median(results_2_back.rt(aa_idx & comp_idx & j_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
subj_stats.twoback.rt_AA_iso = median(results_2_back.rt(aa_idx & iso_idx & j_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
subj_stats.twoback.rt_AA_nov = median(results_2_back.rt(aa_idx & nov_idx & j_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');

% (A-B) rt
subj_stats.twoback.rt_AB_comp = median(results_2_back.rt(ab_idx & comp_idx & k_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
subj_stats.twoback.rt_AB_iso = median(results_2_back.rt(ab_idx & iso_idx & k_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');
subj_stats.twoback.rt_AB_nov = median(results_2_back.rt(ab_idx & nov_idx & k_idx & results_2_back.correct == 1 & valid_rt_idx), 'omitnan');

% for 'A-A' goals (correct = 'j')
n_AA_comp = sum(aa_idx & comp_idx & j_idx);
n_AA_iso = sum(aa_idx & iso_idx & j_idx);
n_AA_nov = sum(aa_idx & nov_idx & j_idx);
subj_stats.twoback.err_AA_comp_as_k = sum(strcmp(results_2_back.resp_key(aa_idx & comp_idx & j_idx), 'k')) / n_AA_comp;
subj_stats.twoback.err_AA_iso_as_k = sum(strcmp(results_2_back.resp_key(aa_idx & iso_idx & j_idx), 'k')) / n_AA_iso;
subj_stats.twoback.err_AA_nov_as_k = sum(strcmp(results_2_back.resp_key(aa_idx & nov_idx & j_idx), 'k')) / n_AA_nov;
subj_stats.twoback.err_AA_comp_as_none = sum(strcmp(results_2_back.resp_key(aa_idx & comp_idx & j_idx), 'none')) / n_AA_comp;
subj_stats.twoback.err_AA_iso_as_none = sum(strcmp(results_2_back.resp_key(aa_idx & iso_idx & j_idx), 'none')) / n_AA_iso;
subj_stats.twoback.err_AA_nov_as_none = sum(strcmp(results_2_back.resp_key(aa_idx & nov_idx & j_idx), 'none')) / n_AA_nov;

% for 'A-B' goals (correct = 'k')
n_AB_comp = sum(ab_idx & comp_idx & k_idx);
n_AB_iso = sum(ab_idx & iso_idx & k_idx);
n_AB_nov = sum(ab_idx & nov_idx & k_idx);
subj_stats.twoback.err_AB_comp_as_j = sum(strcmp(results_2_back.resp_key(ab_idx & comp_idx & k_idx), 'j')) / n_AB_comp;
subj_stats.twoback.err_AB_iso_as_j = sum(strcmp(results_2_back.resp_key(ab_idx & iso_idx & k_idx), 'j')) / n_AB_iso;
subj_stats.twoback.err_AB_nov_as_j = sum(strcmp(results_2_back.resp_key(ab_idx & nov_idx & k_idx), 'j')) / n_AB_nov;
subj_stats.twoback.err_AB_comp_as_none = sum(strcmp(results_2_back.resp_key(ab_idx & comp_idx & k_idx), 'none')) / n_AB_comp;
subj_stats.twoback.err_AB_iso_as_none = sum(strcmp(results_2_back.resp_key(ab_idx & iso_idx & k_idx), 'none')) / n_AB_iso;
subj_stats.twoback.err_AB_nov_as_none = sum(strcmp(results_2_back.resp_key(ab_idx & nov_idx & k_idx), 'none')) / n_AB_nov;

% for 'A-N' goals (correct = 'none')
n_AN_comp = sum(an_idx & comp_idx);
n_AN_iso = sum(an_idx & iso_idx);
n_AN_nov = sum(an_idx & nov_idx);
subj_stats.twoback.err_AN_comp_as_j = sum(strcmp(results_2_back.resp_key(an_idx & comp_idx), 'j')) / n_AN_comp;
subj_stats.twoback.err_AN_iso_as_j = sum(strcmp(results_2_back.resp_key(an_idx & iso_idx), 'j')) / n_AN_iso;
subj_stats.twoback.err_AN_nov_as_j = sum(strcmp(results_2_back.resp_key(an_idx & nov_idx), 'j')) / n_AN_nov;
subj_stats.twoback.err_AN_comp_as_k = sum(strcmp(results_2_back.resp_key(an_idx & comp_idx), 'k')) / n_AN_comp;
subj_stats.twoback.err_AN_iso_as_k = sum(strcmp(results_2_back.resp_key(an_idx & iso_idx), 'k')) / n_AN_iso;
subj_stats.twoback.err_AN_nov_as_k = sum(strcmp(results_2_back.resp_key(an_idx & nov_idx), 'k')) / n_AN_nov;


%% Plots. 2-Back

figure('color', 'white');
h_ax1 = subplot(1, 2, 1);
acc_data = [
    subj_stats.twoback.acc_AA_comp, subj_stats.twoback.acc_AA_iso, subj_stats.twoback.acc_AA_nov;
    subj_stats.twoback.acc_AB_comp, subj_stats.twoback.acc_AB_iso, subj_stats.twoback.acc_AB_nov;
    subj_stats.twoback.acc_AN_comp, subj_stats.twoback.acc_AN_iso, subj_stats.twoback.acc_AN_nov];
b1 = bar(acc_data, 'grouped');
b1(1).FaceColor = plot_color_comp;
b1(2).FaceColor = plot_color_iso; 
b1(3).FaceColor = plot_color_nov;
grid off;
box off;
ylabel('Accuracy', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
set(h_ax1, 'XTickLabel', {'A-A', 'A-B', 'A-N'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
ylim([0, 1.05]);

h_ax2 = subplot(1, 2, 2);
rt_data = [
    subj_stats.twoback.rt_AA_comp, subj_stats.twoback.rt_AA_iso, subj_stats.twoback.rt_AA_nov;
    subj_stats.twoback.rt_AB_comp, subj_stats.twoback.rt_AB_iso, subj_stats.twoback.rt_AB_nov
];
b2 = bar(rt_data, 'grouped');
b2(1).FaceColor = plot_color_comp;
b2(2).FaceColor = plot_color_iso;
b2(3).FaceColor = plot_color_nov;
grid off;
box off;
ylabel('RT (s)', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
set(h_ax2, 'XTickLabel', {'A-A', 'A-B'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
if all(~isnan(rt_data(:))), ylim([0, max(rt_data(:)) * 1.2]); end

sgtitle('2-Back Task', 'FontSize', plot_font_size_sgtitle, 'FontName', plot_font_name);
legend({'Compared', 'Isolated', 'Novel'}, 'FontName', plot_font_name, 'Box', 'off');

% print(gcf, '2Back_Figure_Acc_RT.tiff', '-dtiff', '-r300'); 





%% 2-Back response proportions (A-B trials)
figure('color', 'white');

% Subplot 1: A-A trials
h_ax1 = subplot(1, 3, 1);
Y_AA = [
    subj_stats.twoback.acc_AA_comp, subj_stats.twoback.err_AA_comp_as_k, subj_stats.twoback.err_AA_comp_as_none; ...
    subj_stats.twoback.acc_AA_iso, subj_stats.twoback.err_AA_iso_as_k, subj_stats.twoback.err_AA_iso_as_none; ...
    subj_stats.twoback.acc_AA_nov, subj_stats.twoback.err_AA_nov_as_k, subj_stats.twoback.err_AA_nov_as_none
];
b_AA = bar(Y_AA, 'stacked', 'FaceColor', 'flat');
b_AA(1).FaceColor = plot_color_same;      % Resp: same (bottom)
b_AA(2).FaceColor = plot_color_similar;   % Resp: similar
b_AA(3).FaceColor = plot_color_new;       % Resp: new
box off;
grid off;
title('A-A', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
ylabel('Response Proportions', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
set(h_ax1, 'XTickLabel', {'Compared', 'Isolated', 'Novel'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
xtickangle(h_ax1, 45);
ylim([0, 1.05]);

% Subplot 2: A-B trials
h_ax2 = subplot(1, 3, 2);
Y_AB = [
    subj_stats.twoback.acc_AB_comp, subj_stats.twoback.err_AB_comp_as_j, subj_stats.twoback.err_AB_comp_as_none; ...
    subj_stats.twoback.acc_AB_iso, subj_stats.twoback.err_AB_iso_as_j, subj_stats.twoback.err_AB_iso_as_none; ...
    subj_stats.twoback.acc_AB_nov, subj_stats.twoback.err_AB_nov_as_j, subj_stats.twoback.err_AB_nov_as_none
];
b_AB = bar(Y_AB, 'stacked', 'FaceColor', 'flat');
b_AB(1).FaceColor = plot_color_similar;   % Resp: similar (bottom)
b_AB(2).FaceColor = plot_color_same;      % Resp: same
b_AB(3).FaceColor = plot_color_new;       % Resp: new
box off;
grid off;
title('A-B', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
set(h_ax2, 'XTickLabel', {'Compared', 'Isolated', 'Novel'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
xtickangle(h_ax2, 45);
ylim([0, 1.05]);

% Subplot 3: A-N trials
h_ax3 = subplot(1, 3, 3);
Y_AN = [
    subj_stats.twoback.acc_AN_comp, subj_stats.twoback.err_AN_comp_as_j, subj_stats.twoback.err_AN_comp_as_k; ...
    subj_stats.twoback.acc_AN_iso, subj_stats.twoback.err_AN_iso_as_j, subj_stats.twoback.err_AN_iso_as_k; ...
    subj_stats.twoback.acc_AN_nov, subj_stats.twoback.err_AN_nov_as_j, subj_stats.twoback.err_AN_nov_as_k
];
b_AN = bar(Y_AN, 'stacked', 'FaceColor', 'flat');
b_AN(1).FaceColor = plot_color_new;       % Resp: new (bottom)
b_AN(2).FaceColor = plot_color_same;      % Resp: same
b_AN(3).FaceColor = plot_color_similar;   % Resp: similar
box off;
grid off;
title('A-N', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
set(h_ax3, 'XTickLabel', {'Compared', 'Isolated', 'Novel'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
xtickangle(h_ax3, 45);
ylim([0, 1.05]);

% Create a unified legend
hold(h_ax3, 'on');

h_same_resp = bar(NaN, 'FaceColor', plot_color_same);
h_similar_resp = bar(NaN, 'FaceColor', plot_color_similar);
h_new_resp = bar(NaN, 'FaceColor', plot_color_new);
legend([h_same_resp, h_similar_resp, h_new_resp], ...
    {'Same', 'Similar', 'New'}, 'Location', 'east', 'FontName', plot_font_name, 'Box', 'off');
hold(h_ax3, 'off');

% print(gcf, '2Back_Figure_Resp.tiff', '-dtiff', '-r300'); 





%% 2-Back SDT
% LDI for A-B trials (similar lures vs new items)
% LDI = P(respond "similar" | similar item) - P(respond "similar" | new item)

% 'compared'
hr_AB_comp = subj_stats.twoback.acc_AB_comp;  % P(k | A-B)
far_AB_comp = subj_stats.twoback.err_AN_comp_as_k;  % P(k | A-N)
subj_stats.twoback.ldi_comp = hr_AB_comp - far_AB_comp;

% 'isolated'
hr_AB_iso = subj_stats.twoback.acc_AB_iso;
far_AB_iso = subj_stats.twoback.err_AN_iso_as_k;
subj_stats.twoback.ldi_iso = hr_AB_iso - far_AB_iso;

% 'novel'
hr_AB_nov = subj_stats.twoback.acc_AB_nov;
far_AB_nov = subj_stats.twoback.err_AN_nov_as_k;
subj_stats.twoback.ldi_nov = hr_AB_nov - far_AB_nov;

% Compute d-prime for A-A (identity discrimination)
% 'compared'
hr_AA_comp = subj_stats.twoback.acc_AA_comp;
far_AA_comp = subj_stats.twoback.err_AN_comp_as_j;
if hr_AA_comp == 1, hr_AA_comp = 1 - (1 / (2 * n_AA_comp)); end
if hr_AA_comp == 0, hr_AA_comp = 1 / (2 * n_AA_comp); end
if far_AA_comp == 0, far_AA_comp = 1 / (2 * n_AN_comp); end
if far_AA_comp == 1, far_AA_comp = 1 - (1 / (2 * n_AN_comp)); end
subj_stats.twoback.d_AA_comp = norminv(hr_AA_comp) - norminv(far_AA_comp);

% 'isolated'
hr_AA_iso = subj_stats.twoback.acc_AA_iso;
far_AA_iso = subj_stats.twoback.err_AN_iso_as_j;
if hr_AA_iso == 1, hr_AA_iso = 1 - (1 / (2 * n_AA_iso)); end
if hr_AA_iso == 0, hr_AA_iso = 1 / (2 * n_AA_iso); end
if far_AA_iso == 0, far_AA_iso = 1 / (2 * n_AN_iso); end
if far_AA_iso == 1, far_AA_iso = 1 - (1 / (2 * n_AN_iso)); end
subj_stats.twoback.d_AA_iso = norminv(hr_AA_iso) - norminv(far_AA_iso);

% 'novel'
hr_AA_nov = subj_stats.twoback.acc_AA_nov;
far_AA_nov = subj_stats.twoback.err_AN_nov_as_j;
if hr_AA_nov == 1, hr_AA_nov = 1 - (1 / (2 * n_AA_nov)); end
if hr_AA_nov == 0, hr_AA_nov = 1 / (2 * n_AA_nov); end
if far_AA_nov == 0, far_AA_nov = 1 / (2 * n_AN_nov); end
if far_AA_nov == 1, far_AA_nov = 1 - (1 / (2 * n_AN_nov)); end
subj_stats.twoback.d_AA_nov = norminv(hr_AA_nov) - norminv(far_AA_nov);

% print discrimination results
fprintf('\n--- 2-Back Discrimination Results ---\n');
fprintf('  LDI (Similar Discrimination) | Comp: %.3f, Iso: %.3f, Nov: %.3f\n', ...
    subj_stats.twoback.ldi_comp, subj_stats.twoback.ldi_iso, subj_stats.twoback.ldi_nov);
fprintf('  d'' (Identity Recognition)   | Comp: %.3f, Iso: %.3f, Nov: %.3f\n', ...
    subj_stats.twoback.d_AA_comp, subj_stats.twoback.d_AA_iso, subj_stats.twoback.d_AA_nov);

% Plot LDI and d-prime
figure('color', 'white');

% plot LDI (primary measure for similar lure discrimination)
h_ax1 = subplot(1, 2, 1);
ldi_data = [subj_stats.twoback.ldi_comp, subj_stats.twoback.ldi_iso, subj_stats.twoback.ldi_nov];
b1 = bar(ldi_data, 'FaceColor', 'flat');
b1.CData(1,:) = plot_color_comp;
b1.CData(2,:) = plot_color_iso;
b1.CData(3,:) = plot_color_nov;
grid off;
box off;
title('Similar Lure Discrimination', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
ylabel('LDI', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
set(h_ax1, 'XTickLabel', {'Compared', 'Isolated', 'Novel'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
xtickangle(h_ax1, 45);
ylim([-0.2, 1.0]);
yline(0, '--k', 'LineWidth', 1);

% plot A-A d-prime (identity recognition)
h_ax2 = subplot(1, 2, 2);
d_AA_data = [subj_stats.twoback.d_AA_comp, subj_stats.twoback.d_AA_iso, subj_stats.twoback.d_AA_nov];
b2 = bar(d_AA_data, 'FaceColor', 'flat');
b2.CData(1,:) = plot_color_comp;
b2.CData(2,:) = plot_color_iso;
b2.CData(3,:) = plot_color_nov;
grid off;
box off;
title('Same Item Recognition', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
ylabel('d-prime', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
set(h_ax2, 'XTickLabel', {'Compared', 'Isolated', 'Novel'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
xtickangle(h_ax2, 45);


% print(gcf, '2Back_Figure_LDI_Dprime.tiff', '-dtiff', '-r300'); 




% 
% %%Recognition
% % code correctness
% results_recognition.corr_resp = cellstr(results_recognition.corr_resp);
% results_recognition.resp_key = cellstr(results_recognition.resp_key);
% resp_idx = ~strcmp(results_recognition.resp_key, 'NA');
% num_resp = sum(resp_idx);
% num_total = height(results_recognition);
% results_recognition.correct = strcmp(results_recognition.corr_resp, results_recognition.resp_key);
% 
% % how much did they respond?
% subj_stats.recog.per_resp = num_resp / num_total;
% 
% % how accurately was their response?
% subj_stats.recog.overall_acc = mean(results_recognition.correct(resp_idx));
% 
% % how fast did they respond?
% valid_rt_idx = results_recognition.rt > MIN_RT_CUTOFF;
% all_valid_rts = results_recognition.rt(valid_rt_idx);
% 
% figure('color', 'white');
% histogram(all_valid_rts, 50, 'FaceColor', [0.5 0.5 0.5]);
% title(sprintf('RT Distribution (RTs > %.3fs)', MIN_RT_CUTOFF));
% xlabel('RT (s)');
% ylabel('Count');
% grid on;
% box off;
% 
% % store rt stats
% subj_stats.recog.median_rt = median(all_valid_rts);
% 
% % print sanity checks
% fprintf('percent responded: %f\n', subj_stats.recog.per_resp);
% fprintf('overall accuracy: %f\n', subj_stats.recog.overall_acc);
% fprintf('median reaction times: %f\n', subj_stats.recog.median_rt);
% 
% %% SDT
% % overall hit rate and false alarm rate
% old_trials = results_recognition(results_recognition.trial_type == "old", :);
% new_trials = results_recognition(results_recognition.trial_type ~= "old", :);
% n_old = height(old_trials);
% n_new = height(new_trials);
% n_hits = sum(old_trials.correct == 1 & old_trials.rt > MIN_RT_CUTOFF);
% n_fa = sum(strcmp(new_trials.resp_key, 'j') & new_trials.rt > MIN_RT_CUTOFF);
% 
% HR = n_hits / n_old;
% FAR = n_fa / n_new;
% 
% % correction
% if HR == 1, HR = 1 - (1 / (2 * n_old)); end
% if HR == 0, HR = 1 / (2 * n_old); end 
% if FAR == 0, FAR = 1 / (2 * n_new); end
% if FAR == 1, FAR = 1 - (1 / (2 * n_new)); end 
% 
% % store overall sdt
% subj_stats.recog.HR_overall = HR;
% subj_stats.recog.FAR_overall = FAR;
% subj_stats.recog.d_prime_overall = norminv(HR) - norminv(FAR);
% 
% fprintf('overall hit rate: %f\n', subj_stats.recog.HR_overall);
% fprintf('overall false alarm rate: %f\n', subj_stats.recog.FAR_overall);
% fprintf('overall d''prime: %f\n', subj_stats.recog.d_prime_overall);
% %% Compared vs Isolated
% comp_trials = old_trials(old_trials.condition == "compared", :);
% iso_trials = old_trials(old_trials.condition == "isolated", :);
% n_comp = height(comp_trials);
% n_iso = height(iso_trials);
% 
% % comp hits
% n_hits_comp = sum(comp_trials.correct == 1 & comp_trials.rt > MIN_RT_CUTOFF);
% HR_comp = n_hits_comp / n_comp;
% if HR_comp == 1, HR_comp = 1 - (1 / (2 * n_comp)); end
% if HR_comp == 0, HR_comp = 1 / (2 * n_comp); end
% 
% % iso hits
% n_hits_iso = sum(iso_trials.correct == 1 & iso_trials.rt > MIN_RT_CUTOFF);
% HR_iso = n_hits_iso / n_iso;
% if HR_iso == 1, HR_iso = 1 - (1 / (2 * n_iso)); end
% if HR_iso == 0, HR_iso = 1 / (2 * n_iso); end
% 
% % store conditional sdt
% subj_stats.recog.HR_comp = HR_comp;
% subj_stats.recog.d_prime_comp = norminv(HR_comp) - norminv(FAR);
% subj_stats.recog.HR_iso = HR_iso;
% subj_stats.recog.d_prime_iso = norminv(HR_iso) - norminv(FAR);
% 
% fprintf('compared HR: %f\n', subj_stats.recog.HR_comp);
% fprintf('compared d prime: %f\n', subj_stats.recog.d_prime_comp);
% fprintf('isolated HR: %f\n', subj_stats.recog.HR_iso);
% fprintf('isolated d prime: %f\n', subj_stats.recog.d_prime_iso);
% %% A vs B
% A_trials = old_trials(old_trials.identity == "A", :);
% B_trials = old_trials(old_trials.identity == "B", :);
% n_A = height(A_trials);
% n_B = height(B_trials);
% 
% % A its
% n_hits_A = sum(A_trials.correct == 1 & A_trials.rt > MIN_RT_CUTOFF);
% HR_A = n_hits_A / n_A;
% if HR_A == 1, HR_A = 1 - (1 / (2 * n_A)); end
% if HR_A == 0, HR_A = 1 / (2 * n_A); end
% 
% % B hits
% n_hits_B = sum(B_trials.correct == 1 & B_trials.rt > MIN_RT_CUTOFF);
% HR_B = n_hits_B / n_B;
% if HR_B == 1, HR_B = 1 - (1 / (2 * n_B)); end
% if HR_B == 0, HR_B = 1 / (2 * n_B); end
% 
% % store A vs B sdt
% subj_stats.recog.HR_A = HR_A;
% subj_stats.recog.d_prime_A = norminv(HR_A) - norminv(FAR);
% subj_stats.recog.HR_B = HR_B;
% subj_stats.recog.d_prime_B = norminv(HR_B) - norminv(FAR);
% 
% fprintf('A HR: %f\n', subj_stats.recog.HR_A);
% fprintf('A d prime: %f\n', subj_stats.recog.d_prime_A);
% fprintf('B HR: %f\n', subj_stats.recog.HR_B);
% fprintf('B d prime: %f\n', subj_stats.recog.d_prime_B);
% %% RT
% % Compared vs Isolated vs Foil
% rt_comp_hits = mean(results_recognition.rt(...
%     results_recognition.condition == "compared" & ...
%     results_recognition.correct == 1 & valid_rt_idx), 'omitnan');
% 
% rt_iso_hits = mean(results_recognition.rt(...
%     results_recognition.condition == "isolated" & ...
%     results_recognition.correct == 1 & valid_rt_idx), 'omitnan');
% 
% rt_new_crs = mean(results_recognition.rt(...
%     results_recognition.trial_type ~= "old" & ...
%     results_recognition.correct == 1 & valid_rt_idx), 'omitnan');
% 
% % by item identity
% rt_A_hits = mean(results_recognition.rt(...
%     results_recognition.identity == "A" & ...
%     results_recognition.correct == 1 & valid_rt_idx), 'omitnan');
% 
% rt_B_hits = mean(results_recognition.rt(...
%     results_recognition.identity == "B" & ...
%     results_recognition.correct == 1 & valid_rt_idx), 'omitnan');
% 
% % store rt
% subj_stats.recog.rt_comp_hit = rt_comp_hits;
% subj_stats.recog.rt_iso_hit = rt_iso_hits;
% subj_stats.recog.rt_new_cr = rt_new_crs;
% subj_stats.recog.rt_A_hit = rt_A_hits;
% subj_stats.recog.rt_B_hit = rt_B_hits;
% 
% fprintf('RT compared hit: %f\n', subj_stats.recog.rt_comp_hit);
% fprintf('RT isolated hit: %f\n', subj_stats.recog.rt_iso_hit);
% fprintf('RT new CR: %f\n', subj_stats.recog.rt_new_cr);
% fprintf('RT A hit: %f\n', subj_stats.recog.rt_A_hit);
% fprintf('RT B hit: %f\n', subj_stats.recog.rt_B_hit);
% 
% %% plots: sdt
% figure('color', 'white'); 
% 
% % d' by condition
% h_ax1 = subplot(1, 2, 1); 
% d_prime_data = [subj_stats.recog.d_prime_comp, subj_stats.recog.d_prime_iso];
% b1 = bar(d_prime_data, 'FaceColor', 'flat');
% b1.CData(1,:) = plot_color_comp;
% b1.CData(2,:) = plot_color_iso;
% 
% grid off;
% box off;
% title('d'' by condition', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
% ylabel('d-prime', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
% set(h_ax1, 'XTickLabel', {'Compared', 'Isolated'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
% 
% % % d' by item identity
% % h_ax2 = subplot(1, 2, 2);
% % d_prime_AB_data = [subj_stats.recog.d_prime_A, subj_stats.recog.d_prime_B];
% % b2 = bar(d_prime_AB_data, 'FaceColor', 'flat');
% % b2.CData(1,:) = plot_color_same;
% % b2.CData(2,:) = plot_color_similar;
% % 
% % grid off;
% % box off;
% % title('d'' by item identity', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
% % ylabel('d-prime', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
% % set(h_ax2, 'XTickLabel', {'A', 'B'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
% 
% % --- figure 2: hr vs far ---
% figure('color', 'white');
% hold on;
% b1 = bar(1, subj_stats.recog.HR_comp);
% set(b1, 'FaceColor', plot_color_comp);
% b2 = bar(2, subj_stats.recog.HR_iso);
% set(b2, 'FaceColor', plot_color_iso);
% baseline_val = subj_stats.recog.FAR_overall;
% l1 = yline(baseline_val, 'Color', plot_color_new, 'LineStyle', '--', 'LineWidth', 2);
% 
% grid off;
% box off;
% title('HR vs FAR', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
% ylabel('Proportion', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
% set(gca, 'XTick', [1 2], 'XTickLabel', {'Compared HR', 'Isolated HR'}, ...
%     'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
% legend([b1, b2, l1], {'Compared', 'Isolated', 'FAR'}, ...
%     'Location', 'best', 'FontName', plot_font_name);
% 
% ylim([0, 1.0]);
% % 
% % % --- rt by condition ---
% % h_ax3 = subplot(1, 2, 1);
% % rt_condition_data = [subj_stats.recog.rt_comp_hit, subj_stats.recog.rt_iso_hit, subj_stats.recog.rt_new_cr];
% % b4 = bar(rt_condition_data, 'FaceColor', 'flat');
% % b4.CData(1,:) = plot_color_comp;
% % b4.CData(2,:) = plot_color_iso;
% % b4.CData(3,:) = [0.5 0.5 0.5]; % medium grey
% % 
% % grid off;
% % box off;
% % title('RT by condition (correct)', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
% % ylabel('RT (s)', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
% % set(h_ax3, 'XTickLabel', {'Compared Hit', 'Isolated Hit', 'New CR'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
% % 
% % % --- rt by item identity ---
% % h_ax4 = subplot(1, 2, 2);
% % rt_identity_data = [subj_stats.recog.rt_A_hit, subj_stats.recog.rt_B_hit];
% % b5 = bar(rt_identity_data, 'FaceColor', 'flat');
% % b5.CData(1,:) = plot_color_same;
% % b5.CData(2,:) = plot_color_similar;
% 
% grid off;
% box off;
% title('RT by item identity (correct)', 'FontSize', plot_font_size_title, 'FontName', plot_font_name);
% ylabel('RT (s)', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
% set(h_ax4, 'XTickLabel', {'A', 'B'}, 'FontName', plot_font_name, 'FontSize', plot_font_size_axis);
% 
% 
% % save_file = fullfile(results_dir, sprintf('subj_%03d_stats.mat', subj_id));
% % save(save_file, 'subj_stats');
% 
% 
% 
% 
% 
% %% Analysis: Transfer of Learning (1-Back -> 2-Back)
% % Goal: Calculate P(2-Back Correct | 1-Back Correct) for Compared Lures
% % We strictly look at "Compared" "B" items (the lures).
% 
% % 1. Extract the relevant tables
% % Filter for Compared + Identity B (Lures)
% t1 = results_1_back(results_1_back.condition == "compared" & results_1_back.identity == "B", :);
% t2 = results_2_back(results_2_back.condition == "compared" & results_2_back.identity == "B", :);
% 
% % 2. Clean up columns for joining
% % We only need stim_id and the correctness info
% t1_lite = t1(:, {'stim_id', 'correct'});
% t1_lite.Properties.VariableNames = {'stim_id', 'correct_1back'};
% 
% t2_lite = t2(:, {'stim_id', 'correct'});
% t2_lite.Properties.VariableNames = {'stim_id', 'correct_2back'};
% 
% % 3. Join the tables on 'stim_id'
% % This matches the exact same image file presented in both tasks
% merged_stats = innerjoin(t1_lite, t2_lite, 'Keys', 'stim_id');
% 
% % 4. Calculate Conditional Probabilities
% % Group 1: Items they got RIGHT in 1-Back
% got_right_1back = merged_stats(merged_stats.correct_1back == 1, :);
% acc_2back_given_right = mean(got_right_1back.correct_2back);
% 
% % Group 2: Items they got WRONG in 1-Back
% got_wrong_1back = merged_stats(merged_stats.correct_1back == 0, :);
% acc_2back_given_wrong = mean(got_wrong_1back.correct_2back);
% 
% % Store in structure
% subj_stats.transfer.acc_given_1b_corr = acc_2back_given_right;
% subj_stats.transfer.acc_given_1b_err  = acc_2back_given_wrong;
% subj_stats.transfer.benefit = acc_2back_given_right - acc_2back_given_wrong;
% 
% fprintf('\n--- Transfer of Learning Analysis ---\n');
% fprintf('2-Back Accuracy given 1-Back Correct: %.3f (N=%d)\n', ...
%     acc_2back_given_right, height(got_right_1back));
% fprintf('2-Back Accuracy given 1-Back Error:   %.3f (N=%d)\n', ...
%     acc_2back_given_wrong, height(got_wrong_1back));
% fprintf('Learning Benefit: %+.3f\n', subj_stats.transfer.benefit);

%% Plot: Transfer of Learning
% figure('color', 'white', 'Name', 'Transfer Analysis');
% 
% % Prepare data
% y_data = [subj_stats.transfer.acc_given_1b_corr, subj_stats.transfer.acc_given_1b_err];
% x_labels = {'1-Back: Correct', '1-Back: Error'};
% 
% % Plot
% b = bar(y_data, 'FaceColor', 'flat');
% b.CData(1,:) = plot_color_comp;   % Use your 'compared' color (purple-ish)
% b.CData(2,:) = [0.6 0.6 0.6];     % Grey for the 'Error' baseline
% 
% % Aesthetics
% ylabel('Subsequent 2-Back Accuracy', 'FontSize', plot_font_size_axis, 'FontName', plot_font_name);
% xticklabels(x_labels);
% title({'Does 1-Back Performance Predict', '2-Back Success?'}, ...
%     'FontSize', plot_font_size_title, 'FontName', plot_font_name);
% ylim([0 1.05]);
% grid off; box off;
% 
% % Add text labels on top of bars
% text(1, y_data(1), sprintf('%.2f', y_data(1)), ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
% text(2, y_data(2), sprintf('%.2f', y_data(2)), ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

