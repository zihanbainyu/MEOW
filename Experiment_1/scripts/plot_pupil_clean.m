function plot_pupil_clean(d)
% plot pupil time series (baseline corrected: mean of first 5 pts)
n_r = length(d.pupilTS); n_tt = size(d.TEPR_TT{1}, 1);
lens = cellfun(@(x) size(x,2), d.TEPR_TT); min_len = min(lens);

% aggregate data & apply baseline correction
dat = nan(n_tt, min_len, n_r); sac = nan(n_tt, min_len, n_r);
for r = 1:n_r
    for t = 1:n_tt
        raw = d.TEPR_TT{r}(t, 1:min_len);
        bl = mean(raw(1:1)); % <--- BASELINE: MEAN OF FIRST 5 POINTS
        dat(t,:,r) = raw - bl; 
        sac(t,:,r) = d.sacRate{r}(t, 1:min_len);
    end
end
avg_pup = nanmean(dat, 3); avg_sac = nanmean(sac, 3);
t_ax = linspace(0, (min_len*d.downsampleRate)/d.sampleRate, min_len);

% saccade counts
n_sac = zeros(1, n_tt); n_tr = zeros(1, n_tt);
for r = 1:n_r
    st = d.nSaccsPerTrialType{r}; tt = d.trialTypes{r};
    for i = 1:n_tt
        n_tr(i) = n_tr(i) + sum(tt==i);
        if i<=length(st), n_sac(i) = n_sac(i) + st(i); end
    end
end
sac_rate = n_sac ./ n_tr;

% plotting setup
c_bl = [0 0.447 0.741]; c_pu = [0.494 0.184 0.556]; c_pi = [0.929 0.188 0.525];
cols = [c_bl; c_bl; c_pu; c_pu; c_pi; c_pi];
lsty = {'-'; '--'; '-'; '--'; '-'; '--'};
legs = {'Compared', 'Isolated', 'Novel'};

% fig 1: time series
figure('Name','Pupil TS (BL-Corr)','Position',[100 100 1500 800],'Color','w');
subplot(2,2,1); hold on;
plot(t_ax, avg_pup(1,:), 'LineWidth',2.5, 'Color',c_bl, 'LineStyle','-');
plot(t_ax, avg_pup(3,:), 'LineWidth',2.5, 'Color',c_pu, 'LineStyle','-');
plot(t_ax, avg_pup(5,:), 'LineWidth',2.5, 'Color',c_pi, 'LineStyle','-');
title('A-A (Target)','FontSize',14); xlabel('Time (s)'); ylabel('Pupil Change');
legend(legs,'Location','best'); grid off; box off;

subplot(2,2,2); hold on;
plot(t_ax, avg_pup(2,:), 'LineWidth',2.5, 'Color',c_bl, 'LineStyle','--');
plot(t_ax, avg_pup(4,:), 'LineWidth',2.5, 'Color',c_pu, 'LineStyle','--');
plot(t_ax, avg_pup(6,:), 'LineWidth',2.5, 'Color',c_pi, 'LineStyle','--');
title('A-B (Lure)','FontSize',14); xlabel('Time (s)'); 
legend(legs,'Location','best'); grid off; box off;

% % fig 2: bar chart (avg response)
% avg_scalar = mean(avg_pup, 2);
% dat_bar = [avg_scalar(1) avg_scalar(3) avg_scalar(5); avg_scalar(2) avg_scalar(4) avg_scalar(6)];
% figure('Name','Avg Pupil','Position',[200 200 800 500],'Color','w');
% b = bar(dat_bar); 
% b(1).FaceColor=c_bl; b(2).FaceColor=c_pu; b(3).FaceColor=c_pi;
% set(gca,'XTickLabel',{'A-A','A-B'},'FontSize',12);
% title('Avg Pupil Response (BL-Corr)','FontSize',14); ylabel('Pupil Change'); legend(legs); box off;
% 
% % fig 3: saccades
% figure('Name','Saccades','Position',[150 150 1400 450],'Color','w');
% subplot(1,2,1); 
% b1 = bar([sac_rate(1) sac_rate(3) sac_rate(5)], 'FaceColor','flat');
% b1.CData(1,:)=c_bl; b1.CData(2,:)=c_pu; b1.CData(3,:)=c_pi;
% set(gca,'XTickLabel',legs,'FontSize',12); title('A-A Saccades'); box off;
% 
% subplot(1,2,2);
% b2 = bar([sac_rate(2) sac_rate(4) sac_rate(6)], 'FaceColor','flat');
% b2.CData(1,:)=c_bl; b2.CData(2,:)=c_pu; b2.CData(3,:)=c_pi;
% set(gca,'XTickLabel',legs,'FontSize',12); title('A-B Saccades'); box off;
% sgtitle('Saccade Counts','FontSize',16);
end