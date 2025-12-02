clear; clc; close all;

%% setup
subj_ids = [601, 602, 603, 604, 605, 606, 607];
base_dir = '..'; 
data_dir = fullfile(base_dir, 'data');
out_dir = fullfile(data_dir, 'rec_mod');
min_rt = 0.150;

% colors
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
c_pts = [100 100 100]/255; c_reg = [180 174 211]/255; 
all_subjs = struct();
calc_d = @(h,f,nh,nf) norminv(max(1/(2*nh), min(1-1/(2*nh), h))) - norminv(max(1/(2*nf), min(1-1/(2*nf), f)));



%% math
for s = 1:length(subj_ids)
    curr_id = subj_ids(s);
    fprintf('processing subject %d...\n', curr_id);
   

    load(fullfile(out_dir, sprintf('sub%03d_concat_mod.mat', curr_id)), 'final_data_output');
    
    % --- recognition ---
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


% --- fig 3: recognition ---
figure('color','w','Position',[100 100 600 500]);
d_c = get_v('rec','d_comp'); d_i = get_v('rec','d_iso'); d_tot = (d_c + d_i)/2;
data = [d_tot', d_c', d_i'];
hold on;
fill([-2, 5, 5, -2], [-2, -2, 0, 0], [0.92 0.92 0.92], 'EdgeColor', 'none'); 
yline(0, 'r--', 'Chance', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
raincloud(data, {[0.3 0.3 0.3], c_comp, c_iso}, {'overall','compared','isolated'}, 'd''', 'Do they have memory?');
[~, p, ~, s_stat] = ttest(d_tot);
text(1.5, max(d_tot)*0.9, sprintf('Mean d''=%.2f\nt(%d)=%.2f, p=%.3f', mean(d_tot), s_stat.df, s_stat.tstat, p), ...
     'FontSize',11, 'BackgroundColor','w', 'EdgeColor','k', 'Margin',5);
ylim([-0.5, max(data(:))+0.5]);

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
    ylabel(ylbl,'FontSize',12,'FontWeight','bold'); xlim([0.2, n_grps+0.8]);
    if nargin>5 && ~isempty(ylims), ylim(ylims); end
    grid off; set(gca,'GridAlpha',0.1); box off; hold off;
end
