clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%
% setup
%%%%%%%%%%%%%%%%%%%%%%%
subj_ids = [501, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610];
base_dir = '..'; 
min_rt = 0.150;
% colors
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
c_pts = [100 100 100]/255; c_reg = [180 174 211]/255; 
all_subjs = struct(); 

%%%%%%%%%%%%%%%%%%%%%%%
% math
%%%%%%%%%%%%%%%%%%%%%%%
for s = 1:length(subj_ids)
    curr_id = subj_ids(s);
    fprintf('processing subject %d...\n', curr_id);
    
    subj_dir = fullfile(base_dir, 'data', sprintf('sub%03d', curr_id));
    load(fullfile(subj_dir, sprintf('sub%03d_concat.mat', curr_id)), 'final_data_output');
    
    r1 = final_data_output.results_1_back_all;
    r2 = final_data_output.results_2_back_all;
    stats = [];

    %%%%%%%%%%%%%%%%%%%%%%%
    % 1-back stats
    %%%%%%%%%%%%%%%%%%%%%%%
    r1.resp_key = cellstr(r1.resp_key); r1.resp_key(strcmp(r1.resp_key,'NA')) = {'none'};
    r1.correct = strcmp(cellstr(r1.corr_resp), r1.resp_key);
    
    idx_sim = r1.condition=="compared" & r1.identity=="B";
    idx_sam = r1.condition=="repeat" & strcmp(r1.corr_resp,'j');
    idx_cr = (r1.condition=="compared"&r1.identity=="A")|r1.condition=="isolated"|(r1.condition=="repeat"&strcmp(r1.corr_resp,'none'));
    
    stats.one.acc_sim = mean(r1.correct(idx_sim));
    stats.one.acc_same = mean(r1.correct(idx_sam));
    stats.one.acc_new = mean(r1.correct(idx_cr));
    
    v_rt = r1.rt > min_rt;
    stats.one.rt_sim = median(r1.rt(idx_sim & r1.correct & v_rt), 'omitnan');
    stats.one.rt_same = median(r1.rt(idx_sam & r1.correct & v_rt), 'omitnan');
    
    n_sim = sum(idx_sim); n_sam = sum(idx_sam); n_cr = sum(idx_cr);
    stats.one.err_sim_as_same = sum(strcmp(r1.resp_key(idx_sim),'j'))/n_sim;
    stats.one.err_sim_as_new = sum(strcmp(r1.resp_key(idx_sim),'none'))/n_sim;
    stats.one.err_same_as_sim = sum(strcmp(r1.resp_key(idx_sam),'k'))/n_sam;
    stats.one.err_same_as_new = sum(strcmp(r1.resp_key(idx_sam),'none'))/n_sam;
    stats.one.err_new_as_same = sum(strcmp(r1.resp_key(idx_cr),'j'))/n_cr;
    stats.one.err_new_as_sim = sum(strcmp(r1.resp_key(idx_cr),'k'))/n_cr;

    %%%%%%%%%%%%%%%%%%%%%%%
    % 2-back stats
    %%%%%%%%%%%%%%%%%%%%%%%
    r2.resp_key = cellstr(r2.resp_key); r2.resp_key(strcmp(r2.resp_key,'NA')) = {'none'};
    r2.correct = strcmp(cellstr(r2.corr_resp), r2.resp_key);
    
    pan = zeros(height(r2),1); 
    for i=1:height(r2)-2, if strcmp(r2.goal(i),'A-N'), pan(i+2)=1; end; end % the corresponding item for A-N trials
    
    real = ~contains(r2.goal, "JUNK");
    v_rt = r2.rt > min_rt;
    aa_idx = real & strcmp(r2.goal,'A-A'); 
    ab_idx = real & strcmp(r2.goal,'A-B'); 
    an_idx = (pan==1);
    comp_idx = real & strcmp(r2.condition,'compared');
    iso_idx = real & strcmp(r2.condition,'isolated');
    nov_idx = real & strcmp(r2.condition,'novel');
    j_idx = real & strcmp(r2.corr_resp,'j'); 
    k_idx = real & strcmp(r2.corr_resp,'k');
    n_idx = real & strcmp(r2.corr_resp, 'none');

    calc_d = @(h,f,nh,nf) norminv(max(1/(2*nh), min(1-1/(2*nh), h))) - norminv(max(1/(2*nf), min(1-1/(2*nf), f)));

    stats.two.acc_AA_comp = mean(r2.correct(aa_idx & comp_idx & j_idx));
    stats.two.acc_AA_iso = mean(r2.correct(aa_idx & iso_idx & j_idx));
    stats.two.acc_AA_nov = mean(r2.correct(aa_idx & nov_idx & j_idx));
    stats.two.acc_AB_comp = mean(r2.correct(ab_idx & comp_idx & k_idx));
    stats.two.acc_AB_iso = mean(r2.correct(ab_idx & iso_idx & k_idx));
    stats.two.acc_AB_nov = mean(r2.correct(ab_idx & nov_idx & k_idx));
    stats.two.acc_AN_comp = mean(r2.correct(an_idx & comp_idx & n_idx));
    stats.two.acc_AN_iso = mean(r2.correct(an_idx & iso_idx & n_idx));
    stats.two.acc_AN_nov = mean(r2.correct(an_idx & nov_idx & n_idx));

    stats.two.rt_AA_comp = median(r2.rt(aa_idx & comp_idx & j_idx & r2.correct==1 & v_rt), 'omitnan');
    stats.two.rt_AA_iso = median(r2.rt(aa_idx & iso_idx & j_idx & r2.correct==1 & v_rt), 'omitnan');
    stats.two.rt_AA_nov = median(r2.rt(aa_idx & nov_idx & j_idx & r2.correct==1 & v_rt), 'omitnan');
    stats.two.rt_AB_comp = median(r2.rt(ab_idx & comp_idx & k_idx & r2.correct==1 & v_rt), 'omitnan');
    stats.two.rt_AB_iso = median(r2.rt(ab_idx & iso_idx & k_idx & r2.correct==1 & v_rt), 'omitnan');
    stats.two.rt_AB_nov = median(r2.rt(ab_idx & nov_idx & k_idx & r2.correct==1 & v_rt), 'omitnan');

    n_AA_comp = sum(aa_idx & comp_idx & j_idx);
    n_AA_iso = sum(aa_idx & iso_idx & j_idx);
    n_AA_nov = sum(aa_idx & nov_idx & j_idx);
    n_AB_comp = sum(ab_idx & comp_idx & k_idx);
    n_AB_iso = sum(ab_idx & iso_idx & k_idx);
    n_AB_nov = sum(ab_idx & nov_idx & k_idx);
    n_AN_comp = sum(an_idx & comp_idx & n_idx);
    n_AN_iso = sum(an_idx & iso_idx & n_idx);
    n_AN_nov = sum(an_idx & nov_idx & n_idx);

    stats.two.err_AA_comp_as_k = sum(strcmp(r2.resp_key(aa_idx & comp_idx & j_idx), 'k')) / n_AA_comp;
    stats.two.err_AA_comp_as_n = sum(strcmp(r2.resp_key(aa_idx & comp_idx & j_idx), 'none')) / n_AA_comp;
    stats.two.err_AB_comp_as_j = sum(strcmp(r2.resp_key(ab_idx & comp_idx & k_idx), 'j')) / n_AB_comp;
    stats.two.err_AB_comp_as_n = sum(strcmp(r2.resp_key(ab_idx & comp_idx & k_idx), 'none')) / n_AB_comp;
    stats.two.err_AN_comp_as_j = sum(strcmp(r2.resp_key(an_idx & comp_idx & n_idx), 'j')) / n_AN_comp;
    stats.two.err_AN_comp_as_k = sum(strcmp(r2.resp_key(an_idx & comp_idx & n_idx), 'k')) / n_AN_comp;

    stats.two.err_AA_iso_as_k = sum(strcmp(r2.resp_key(aa_idx & iso_idx & j_idx), 'k')) / n_AA_iso;
    stats.two.err_AA_iso_as_n = sum(strcmp(r2.resp_key(aa_idx & iso_idx & j_idx), 'none')) / n_AA_iso;
    stats.two.err_AB_iso_as_j = sum(strcmp(r2.resp_key(ab_idx & iso_idx & k_idx), 'j')) / n_AB_iso;
    stats.two.err_AB_iso_as_n = sum(strcmp(r2.resp_key(ab_idx & iso_idx & k_idx), 'none')) / n_AB_iso;
    stats.two.err_AN_iso_as_j = sum(strcmp(r2.resp_key(an_idx & iso_idx & n_idx), 'j')) / n_AN_iso;
    stats.two.err_AN_iso_as_k = sum(strcmp(r2.resp_key(an_idx & iso_idx & n_idx), 'k')) / n_AN_iso;
    
    stats.two.err_AA_nov_as_k = sum(strcmp(r2.resp_key(aa_idx & nov_idx & j_idx), 'k')) / n_AA_nov;
    stats.two.err_AA_nov_as_n = sum(strcmp(r2.resp_key(aa_idx & nov_idx & j_idx), 'none')) / n_AA_nov;
    stats.two.err_AB_nov_as_j = sum(strcmp(r2.resp_key(ab_idx & nov_idx & k_idx), 'j')) / n_AB_nov;
    stats.two.err_AB_nov_as_n = sum(strcmp(r2.resp_key(ab_idx & nov_idx & k_idx), 'none')) / n_AB_nov;
    stats.two.err_AN_nov_as_j = sum(strcmp(r2.resp_key(an_idx & nov_idx & n_idx), 'j')) / n_AN_nov;
    stats.two.err_AN_nov_as_k = sum(strcmp(r2.resp_key(an_idx & nov_idx & n_idx), 'k')) / n_AN_nov;
    
    stats.two.ldi_comp = stats.two.acc_AB_comp - stats.two.err_AN_comp_as_k;
    stats.two.ldi_iso = stats.two.acc_AB_iso - stats.two.err_AN_iso_as_k;
    stats.two.ldi_nov = stats.two.acc_AB_nov - stats.two.err_AN_nov_as_k;
    stats.two.dprime_comp = calc_d(stats.two.acc_AA_comp, stats.two.err_AN_comp_as_j, n_AA_comp, n_AN_comp);
    stats.two.dprime_iso = calc_d(stats.two.acc_AA_iso, stats.two.err_AN_iso_as_j, n_AA_iso, n_AN_iso);
    stats.two.dprime_nov = calc_d(stats.two.acc_AA_nov, stats.two.err_AN_nov_as_j, n_AA_nov, n_AN_nov);

    
    %%%%%%%%%%%%%%%%%%%%%%%
    % recognition stats
    %%%%%%%%%%%%%%%%%%%%%%%
    rec = final_data_output.results_recognition;
    rec.correct = strcmp(cellstr(rec.corr_resp), cellstr(rec.resp_key));
    old = rec(rec.trial_type=="old",:); new = rec(rec.trial_type~="old",:);
    n_new = height(new); n_fa = sum(strcmp(cellstr(new.resp_key),'j') & new.rt>min_rt);
    far = max(1/(2*n_new), min(1-1/(2*n_new), n_fa/n_new));

    tc = old(old.condition=="compared",:); nc = height(tc);
    hc = sum(tc.correct & tc.rt>min_rt)/nc;
    stats.rec.d_comp = calc_d(hc, far, nc, n_new);
    ti = old(old.condition=="isolated",:); ni = height(ti);
    hi = sum(ti.correct & ti.rt>min_rt)/ni;
    stats.rec.d_iso = calc_d(hi, far, ni, n_new);
    all_subjs(s).stats = stats;
        
end
get_v = @(f1, f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs);

%% visualization

%%%%%%%%%%%%%%%%%%%%%%%
% fig. 1-back
%%%%%%%%%%%%%%%%%%%%%%%
% figure('color','w','Position',[100 100 1200 400]); axis square;
% subplot(1,3,1); 
% data = [get_v('one','acc_same'); get_v('one','acc_sim'); get_v('one','acc_new')]';
% raincloud(data, {c_same, c_sim, c_new}, {'Hit(Same)','Hit(Sim)','CR(New)'}, 'Accuracy', '', [0,1]);
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(1,3,2); 
% data = [get_v('one','rt_same'); get_v('one','rt_sim')]';
% raincloud(data, {c_same, c_sim}, {'Hit(Same)','Hit(Sim)'}, 'RT (s)', '');
% add_sig(data, [1 2]);
% subplot(1,3,3); hold on;
% mat_1_back = [mean(get_v('one','acc_same')), mean(get_v('one','err_same_as_sim')), mean(get_v('one','err_same_as_new'));
%        mean(get_v('one','err_sim_as_same')), mean(get_v('one','acc_sim')), mean(get_v('one','err_sim_as_new'));
%        mean(get_v('one','err_new_as_same')), mean(get_v('one','err_new_as_sim')), mean(get_v('one','acc_new'))];
% draw_matrix(mat_1_back, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('Response matrix', 'FontSize', 12);
% sgtitle('1-Back Task Performance', 'FontSize', 16);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, '1Back_Figures.tiff', '-dtiff', '-r300'); 

%%%%%%%%%%%%%%%%%%%%%%%
% fig. 2-back
%%%%%%%%%%%%%%%%%%%%%%%
% figure('color','w','Position',[50 50 1200 900]);

%%%%%%%%%%%%%%%%%%%%%%%
% ldi & dprime & rt
%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,3,1);
% data = [get_v('two','ldi_comp'); get_v('two','ldi_iso'); get_v('two','ldi_nov')]';
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'LDI', 'Lure Discrimination');
% add_sig(data, [1 2; 2 3; 1 3]);
% 
% subplot(3,3,4);
% data = [get_v('two','dprime_comp'); get_v('two','dprime_iso'); get_v('two','dprime_nov')]';
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'd''', 'Recognition');
% add_sig(data, [1 2; 2 3; 1 3]);
% 
% subplot(3,3,2);
% data = [get_v('two','rt_AB_comp'); get_v('two','rt_AB_iso'); get_v('two','rt_AB_nov')]';
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Lure Discrimination)');
% add_sig(data, [1 2; 2 3; 1 3]);
% 
% subplot(3,3,5);
% data = [get_v('two','rt_AA_comp'); get_v('two','rt_AA_iso'); get_v('two','rt_AA_nov')]';
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Recognition)');
% add_sig(data, [1 2; 2 3; 1 3]);

%%%%%%%%%%%%%%%%%%%%%%%
% confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,3,7); hold on;
% mat_comp = [mean(get_v('two','acc_AA_comp')), mean(get_v('two','err_AA_comp_as_k')), mean(get_v('two','err_AA_comp_as_n'));
%        mean(get_v('two','err_AB_comp_as_j')), mean(get_v('two','acc_AB_comp')), mean(get_v('two','err_AB_comp_as_n'));
%        mean(get_v('two','err_AN_comp_as_j')), mean(get_v('two','err_AN_comp_as_k')), mean(get_v('two','acc_AN_comp'))];
% draw_matrix(mat_comp, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('compared', 'FontSize', 14);
% subplot(3,3,8); hold on;
% mat_iso = [mean(get_v('two','acc_AA_iso')), mean(get_v('two','err_AA_iso_as_k')), mean(get_v('two','err_AA_iso_as_n'));
%        mean(get_v('two','err_AB_iso_as_j')), mean(get_v('two','acc_AB_iso')), mean(get_v('two','err_AB_iso_as_n'));
%        mean(get_v('two','err_AN_iso_as_j')), mean(get_v('two','err_AN_iso_as_k')), mean(get_v('two','acc_AN_iso'))];
% draw_matrix(mat_iso, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('isolated', 'FontSize', 14);
% subplot(3,3,9); hold on;
% mat_nov = [mean(get_v('two','acc_AA_nov')), mean(get_v('two','err_AA_nov_as_k')), mean(get_v('two','err_AA_nov_as_n'));
%        mean(get_v('two','err_AB_nov_as_j')), mean(get_v('two','acc_AB_nov')), mean(get_v('two','err_AB_nov_as_n'));
%        mean(get_v('two','err_AN_nov_as_j')), mean(get_v('two','err_AN_nov_as_k')), mean(get_v('two','acc_AN_nov'))];
% draw_matrix(mat_nov, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('novel', 'FontSize', 14);
% sgtitle('2-Back Task Performance', 'FontSize', 16);
% print(gcf, '2Back_Figures.tiff', '-dtiff', '-r300'); 

%%%%%%%%%%%%%%%%%%%%%%%
% predicts 2-back performance from 1-back
%%%%%%%%%%%%%%%%%%%%%%%
x_rt_sim  = get_v('one', 'rt_sim')';
x_acc = get_v('one','acc_sim')';
y_d_c = get_v('two','dprime_comp')';
y_d_i = get_v('two','dprime_iso')';
y_d_n = get_v('two','dprime_nov')';
y_d_o = (y_d_c+y_d_i)/2;
y_l_c = get_v('two','ldi_comp')'; 
y_l_i = get_v('two','ldi_iso')'; 
y_l_n = get_v('two','ldi_nov')'; 
y_l_o = (y_l_c+y_l_i)/2;
y_r_c = get_v('two','rt_AB_comp')';
y_r_i = get_v('two','rt_AB_iso')';
y_r_n = get_v('two','rt_AB_nov')';

% acc ~ ldi
[r_lo, p_lo] = corr(x_acc, y_l_o); 
[r_lc, p_lc] = corr(x_acc, y_l_c); 
[r_li, p_li] = corr(x_acc, y_l_i);
[r_ln, p_ln] = corr(x_acc, y_l_n);
% rt ~ ldi/dprime
[r_rlc, p_rlc] = corr(x_rt_sim, y_l_c); 
[r_rli, p_rli] = corr(x_rt_sim, y_l_i);
[r_rln, p_rln] = corr(x_rt_sim, y_l_n);
[r_rlo, p_rlo] = corr(x_rt_sim, y_l_o);
[r_do, p_do] = corr(x_rt_sim, y_d_o);
[r_dc, p_dc] = corr(x_rt_sim, y_d_c);
[r_di, p_di] = corr(x_rt_sim, y_d_i);
% rt ~ rt
[r_rrc, p_rrc] = corr(x_rt_sim, y_r_c); 
[r_rri, p_rri] = corr(x_rt_sim, y_r_i); 
[r_rrn, p_rrn] = corr(x_rt_sim, y_r_n); 


figure('color','w','Position',[100 100 800 400]);
subplot(1,2,1); hold on; 
s_i2 = plot_layer(x_acc, y_l_i, c_iso, 60, 0.5, 2);
s_c2 = plot_layer(x_acc, y_l_c, c_comp, 60, 0.5, 2);
s_o2 = plot_layer(x_acc, y_l_o, [0.2 0.2 0.2], 60, 1, 2.5);
s_n2 = plot_layer(x_acc, y_l_n, c_nov, 60, 0.5, 2);
xlabel('1-Back Accuracy (Similar)','FontSize',14,'FontWeight','bold'); ylabel('2-Back LDI','FontSize',14,'FontWeight','bold');
title('1-Back Comparison Specifically Predicts Lure Discrimination','FontSize',12);
legend([s_o2, s_c2, s_i2, s_n2], {sprintf('compared+isolated (r=%.2f, p=%.3f)',r_lo,p_lo), ...
    sprintf('compared (r=%.2f, p=%.3f)',r_lc,p_lc), sprintf('isolated (r=%.2f, p=%.3f)',r_li,p_li), sprintf('novel (r=%.2f, p=%.3f)',r_ln,p_ln)}, ...
    'Location','southeast','FontSize',10);

subplot(1,2,2); hold on; 
s_i = plot_layer(x_rt_sim, y_l_i, c_iso, 60, 0.5, 2); 
s_c = plot_layer(x_rt_sim, y_l_c, c_comp, 60, 0.5, 2);
s_n = plot_layer(x_rt_sim, y_l_n, c_nov, 60, 0.5, 2);
s_o = plot_layer(x_rt_sim, y_l_o, [0.2 0.2 0.2], 60, 1, 2.5);
xlabel('1-Back RT (s)','FontSize',14,'FontWeight','bold'); ylabel('2-Back LDI','FontSize',14,'FontWeight','bold');
title('1-Back Faster Comparison Predicts Lure Discrimination','FontSize',12);
legend([s_o, s_c, s_i, s_n], {sprintf('compared+isolated (r=%.2f, p=%.3f)',r_rlo,p_rlo), ...
    sprintf('compared (r=%.2f, p=%.3f)',r_rlc,p_rlc), sprintf('isolated (r=%.2f, p=%.3f)',r_rli,p_rli), sprintf('novel (r=%.2f, p=%.3f)', r_rln,p_rln)}, ...
    'Location','northeast','FontSize',10);
grid off; set(gca,'GridAlpha',0.1); box off;


% subplot(2,3,3); hold on; 
% s_c = plot_layer(x_rt_sim, y_d_c, c_comp, 60, 0.5, 2);
% xlabel('1-Back RT (s)','FontSize',14,'FontWeight','bold'); ylabel('2-Back d''','FontSize',14,'FontWeight','bold');
% title('1-Back Comparison Speed Predicts 2-Back Recognition','FontSize',12);
% legend(s_c, {sprintf('compared (r=%.2f, p=%.3f)',r_dc,p_dc)}, ...
%     'Location','northeast','FontSize',10);
% grid off; set(gca,'GridAlpha',0.1); box off;
% 
% subplot(2,3,4); hold on; 
% s_c = plot_layer(x_rt_sim, y_r_c, c_comp, 60, 0.5, 2);
% s_i = plot_layer(x_rt_sim, y_r_i, c_iso, 60, 0.5, 2);
% s_n = plot_layer(x_rt_sim, y_r_n, c_nov, 60, 0.5, 2);
% xlabel('1-Back RT (s)','FontSize',14,'FontWeight','bold'); ylabel('2-Back RT (s)','FontSize',14,'FontWeight','bold');
% title('1-Back Comparison Speed Predicts 2-Back Recognition','FontSize',12);
% legend([s_c, s_i, s_n], {sprintf('compared (r=%.2f, p=%.3f)',r_rrc,p_rrc), ...
%     sprintf('isolated (r=%.2f, p=%.3f)',r_rri,p_rri), sprintf('novel (r=%.2f, p=%.3f)',r_rrn,p_rrn)}, ...
%     'Location','northeast','FontSize',10);
% grid off; set(gca,'GridAlpha',0.1); box off;

sgtitle('Predicting 2-back performance from 1-back','FontSize',16);


% subplot(2,3,3); hold on; 
% xline(0,'--','Color',[0.8 0.8 0.8]); yline(0,'--','Color',[0.8 0.8 0.8]);
% s_i2 = plot_layer(x_res, y_li_r, c_iso, 60, 0.5, 2);
% s_c2 = plot_layer(x_res, y_lc_r, c_comp, 60, 0.5, 2);
% s_o2 = plot_layer(x_res, y_lo_r, [0.2 0.2 0.2], 60, 1, 2.5);
% xlabel('1-Back Accuracy (Residuals)','FontSize',14,'FontWeight','bold'); ylabel('2-Back LDI (Residuals)','FontSize',14,'FontWeight','bold');
% title({'Comparison Memory Drives Lure Discrimination', '(controlling for speed)'}, 'FontSize',11);
% legend([s_o2, s_c2, s_i2], {sprintf('overall (r_{p}=%.2f, p=%.3f)',r_plo,p_plo), ...
%     sprintf('compared (r_{p}=%.2f, p=%.3f)',r_plc,p_plc), sprintf('isolated (r_{p}=%.2f, p=%.3f)',r_pli,p_pli)}, ...
%     'Location','southeast','FontSize',10);
% 

% partial regression: controlling for speed
% Z = x_rt; X = x_acc;
% mdl_x = fitlm(Z, X); x_res = mdl_x.Residuals.Raw;
% get_res = @(y) fitlm(Z, y).Residuals.Raw;
% y_dc_r = get_res(y_d_c); y_di_r = get_res(y_d_i); y_do_r = get_res(y_d_o);
% y_lc_r = get_res(y_l_c); y_li_r = get_res(y_l_i); y_lo_r = get_res(y_l_o);
% 
% [r_pdo, p_pdo] = partialcorr(X, y_d_o, Z); [r_pdc, p_pdc] = partialcorr(X, y_d_c, Z); [r_pdi, p_pdi] = partialcorr(X, y_d_i, Z);
% [r_plo, p_plo] = partialcorr(X, y_l_o, Z); [r_plc, p_plc] = partialcorr(X, y_l_c, Z); [r_pli, p_pli] = partialcorr(X, y_l_i, Z);
% 



% print(gcf, '2Back_1Back_Figures.tiff', '-dtiff', '-r300'); 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% % episodic memory
% %%%%%%%%%%%%%%%%%%%%%%%
% figure('color','w','Position',[100 100 600 500]);
% d_c = get_v('rec','d_comp'); d_i = get_v('rec','d_iso'); d_tot = (d_c+d_i)/2;
% data = [d_tot', d_c', d_i'];
% hold on;
% fill([-2, 5, 5, -2], [-2, -2, 0, 0], [0.92 0.92 0.92], 'EdgeColor', 'none'); 
% yline(0, 'r--', 'Chance', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
% raincloud(data, {[0.3 0.3 0.3], c_comp, c_iso}, {'overall','compared','isolated'}, 'd''', 'Validation of Episodic Encoding');
% add_sig(data, [1 0; 2 3]);
% [~,p_t,~,s_t] = ttest(d_tot); [~,p_d,~,s_d] = ttest(d_c, d_i);
% msg = {sprintf('\\bfoverall > 0:\\rm t(%d)=%.2f, p=%.3f', s_t.df, s_t.tstat, p_t), ...
%        sprintf('\\bfcompared vs isolated:\\rm t(%d)=%.2f, p=%.3f', s_d.df, s_d.tstat, p_d)};
% text(3.4, max(data(:))*1.15, msg, 'FontSize',10, 'BackgroundColor','w', 'EdgeColor','k', ...
%      'Margin',5, 'HorizontalAlignment','right', 'VerticalAlignment','top');
% ylim([-0.5, max(data(:))*1.3]);
% print(gcf, 'Recog_Figures.tiff', '-dtiff', '-r300'); 


%% functions
function raincloud(mat, cols, xlbls, ylbl, ttl, ylims)
    [n_rows, n_grps] = size(mat); hold on;
    jit = -0.15 - (rand(size(mat)) * 0.20); x_c = repmat(1:n_grps, n_rows, 1) + jit;
    plot(x_c', mat', '-', 'Color', [0.7 0.7 0.7, 0.4], 'LineWidth', 0.5);
    for i = 1:n_grps
        d = mat(:,i); d = d(~isnan(d)); [f, xi] = ksdensity(d); f = f/max(f)*0.4;
        patch([i+f, i*ones(1,length(f))], [xi, fliplr(xi)], cols{i}, 'EdgeColor','none','FaceAlpha',0.5);
        scatter(x_c(:,i), mat(:,i), 20, cols{i}, 'filled', 'MarkerFaceAlpha',0.6);
        q = quantile(d, [0.25, 0.5, 0.75]);
        rectangle('Position', [i-0.06, q(1), 0.12, q(3)-q(1)], 'FaceColor', cols{i}, 'EdgeColor','k', 'LineWidth',1.2);
        plot([i-0.06, i+0.06], [q(2), q(2)], 'k-', 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:n_grps, 'XTickLabel', xlbls, 'FontSize', 12);
    if ~isempty(ttl), title(ttl, 'FontSize', 14); end
    ylabel(ylbl,'FontSize',14,'FontWeight','bold'); xlim([0.2, n_grps+0.8]);
    if nargin>5 && ~isempty(ylims), ylim(ylims); end
    grid off; set(gca,'GridAlpha',0.1); box off; hold off;
end

function draw_matrix(mat, cols, ylbl, xlbl)
    hold on;
    for r=1:3
        for c=1:3
            v = mat(r,c); 
            t_c = (v>0.5)*[1 1 1] + (v<=0.5)*[0.2 0.2 0.2];           
            patch([c-0.48 c+0.48 c+0.48 c-0.48], [r-0.48 r-0.48 r+0.48 r+0.48], ...
                  cols{r}, 'EdgeColor','none','FaceAlpha',v);
            text(c, r, sprintf('%.2f',v), 'HorizontalAlignment','center', ...
                 'Color',t_c, 'FontWeight','bold', 'FontSize',10);
            if r==1 
                text(c, 0.4, xlbl{c}, 'HorizontalAlignment','center', ...
                     'Color',cols{c}, 'FontWeight','bold', 'FontSize',10); 
            end
        end
        text(0.4, r, ylbl{r}, 'HorizontalAlignment','right', ...
             'Color',cols{r}, 'FontWeight','bold', 'FontSize',10);
    end 
    axis ij equal; 
    xlim([0.2 3.8]); ylim([0.3 3.5]);
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
    hold off;
end

function s = plot_layer(x, y, c, sz, alph, lw)
    s = scatter(x, y, sz, c, 'filled', 'MarkerFaceAlpha', alph);
    p = polyfit(x, y, 1); x_r = linspace(min(x), max(x), 100);
    plot(x_r, polyval(p, x_r), 'Color', c, 'LineWidth', lw);
end

function add_sig(data, pairs)
    % pairs: [col_A col_B] -> paired t-test between A and B
    %        [col_A 0]     -> one-sample t-test of A against 0
    [~, n_grps] = size(data);
    y_max = max(data(:)); 
    rng = max(data(:)) - min(data(:));
    if rng == 0, rng = 1; end 
    step = rng * 0.15; 

    line_lvl = 1;
    hold on;
    for i = 1:size(pairs,1)
        c1 = pairs(i,1); c2 = pairs(i,2);
        if c2 == 0 
            [~, p] = ttest(data(:,c1)); 
            is_paired = false;
        else
            [~, p] = ttest(data(:,c1), data(:,c2));
            is_paired = true;
        end
        if p < 0.05
            if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
            y_pos = y_max + (line_lvl * step);
            if is_paired
                plot([c1, c1, c2, c2], [y_pos-step*0.2, y_pos, y_pos, y_pos-step*0.2], 'k-', 'LineWidth', 1.2);
                text(mean([c1 c2]), y_pos + step*0.05, txt, 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold');
            else
                text(c1, y_pos, txt, 'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
            end

            line_lvl = line_lvl + 1;
        end
    end
    if line_lvl > 1
        current_ylim = ylim;
        new_upper = y_max + (line_lvl * step) + step;
        if new_upper > current_ylim(2)
            ylim([current_ylim(1), new_upper]);
        end
    end
    hold off;
end
