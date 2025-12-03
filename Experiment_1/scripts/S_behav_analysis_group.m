clear; clc; close all;

%% setup
subj_ids = [601, 602, 603, 604, 605, 606, 607, 608, 609, 610];
base_dir = '..'; 
min_rt = 0.150;

% colors
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
c_pts = [100 100 100]/255; c_reg = [180 174 211]/255; 
all_subjs = struct(); 

%% math
for s = 1:length(subj_ids)
    curr_id = subj_ids(s);
    fprintf('processing subject %d...\n', curr_id);
    
    f_dir = fullfile(base_dir, 'data', sprintf('sub%03d', curr_id));
    if ~exist(fullfile(f_dir, sprintf('sub%03d_concat.mat', curr_id)), 'file'), continue; end
    load(fullfile(f_dir, sprintf('sub%03d_concat.mat', curr_id)), 'final_data_output');
    
    r1 = final_data_output.results_1_back_all;
    r2 = final_data_output.results_2_back_all;
    stats = [];

    % --- 1-back stats ---
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

    % --- 2-back stats ---
    r2.resp_key = cellstr(r2.resp_key); r2.resp_key(strcmp(r2.resp_key,'NA')) = {'none'};
    r2.correct = strcmp(cellstr(r2.corr_resp), r2.resp_key);
    
    pan = zeros(height(r2),1); 
    for i=1:height(r2)-2, if strcmp(r2.goal(i),'A-N'), pan(i+2)=1; end; end
    
    real = ~contains(r2.goal, "JUNK"); v_rt = r2.rt > min_rt;
    aa = real & strcmp(r2.goal,'A-A'); ab = real & strcmp(r2.goal,'A-B'); an = (pan==1);
    comp = strcmp(r2.condition,'compared'); iso = strcmp(r2.condition,'isolated'); nov = strcmp(r2.condition,'novel');
    j_key = strcmp(r2.corr_resp,'j'); k_key = strcmp(r2.corr_resp,'k');

    calc_d = @(h,f,nh,nf) norminv(max(1/(2*nh), min(1-1/(2*nh), h))) - norminv(max(1/(2*nf), min(1-1/(2*nf), f)));

    % condition loop
    conds = {comp, iso, nov}; suffixes = {'comp', 'iso', 'nov'};
    for c=1:3
        msk = conds{c}; sx = suffixes{c};
        
        % metrics
        stats.two.(['acc_AA_' sx]) = mean(r2.correct(aa & msk & j_key));
        stats.two.(['acc_AB_' sx]) = mean(r2.correct(ab & msk & k_key));
        stats.two.(['rt_AA_' sx]) = median(r2.rt(aa & msk & j_key & r2.correct & v_rt), 'omitnan');
        stats.two.(['rt_AB_' sx]) = median(r2.rt(ab & msk & k_key & r2.correct & v_rt), 'omitnan');
        err_k = sum(strcmp(r2.resp_key(an & msk), 'k')) / sum(an & msk);
        stats.two.(['ldi_' sx]) = stats.two.(['acc_AB_' sx]) - err_k;
        nh = sum(aa & msk & j_key); nfa = sum(an & msk);
        hr = stats.two.(['acc_AA_' sx]); far = sum(strcmp(r2.resp_key(an & msk), 'j')) / nfa;
        stats.two.(['d_AA_' sx]) = calc_d(hr, far, nh, nfa);
        
        % Response proportions for visualization (matching your plot code)
        % A-A trials: correct = 'j' (same), error = 'k' (similar) or 'none' (new)
        stats.two.(['err_AA_' sx '_as_k']) = sum(strcmp(r2.resp_key(aa & msk), 'k')) / sum(aa & msk);
        stats.two.(['err_AA_' sx '_as_none']) = sum(strcmp(r2.resp_key(aa & msk), 'none')) / sum(aa & msk);
        
        % A-B trials: correct = 'k' (similar), error = 'j' (same) or 'none' (new)
        stats.two.(['err_AB_' sx '_as_j']) = sum(strcmp(r2.resp_key(ab & msk), 'j')) / sum(ab & msk);
        stats.two.(['err_AB_' sx '_as_none']) = sum(strcmp(r2.resp_key(ab & msk), 'none')) / sum(ab & msk);
        
        % A-N trials: correct = 'none' (new), error = 'j' (same) or 'k' (similar)
        stats.two.(['acc_AN_' sx]) = sum(strcmp(r2.resp_key(an & msk), 'none')) / sum(an & msk);
        stats.two.(['err_AN_' sx '_as_j']) = sum(strcmp(r2.resp_key(an & msk), 'j')) / sum(an & msk);
        stats.two.(['err_AN_' sx '_as_k']) = sum(strcmp(r2.resp_key(an & msk), 'k')) / sum(an & msk);
        
        % Optional: Keep confusion matrix with generic labels if needed elsewhere
        daa = sum(aa & msk); dab = sum(ab & msk); dan = sum(an & msk);
        stats.two.(['m_' sx '_aa_j']) = sum(strcmp(r2.resp_key(aa & msk),'j'))/daa;
        stats.two.(['m_' sx '_aa_k']) = sum(strcmp(r2.resp_key(aa & msk),'k'))/daa;
        stats.two.(['m_' sx '_aa_n']) = sum(strcmp(r2.resp_key(aa & msk),'none'))/daa;
        stats.two.(['m_' sx '_ab_j']) = sum(strcmp(r2.resp_key(ab & msk),'j'))/dab;
        stats.two.(['m_' sx '_ab_k']) = sum(strcmp(r2.resp_key(ab & msk),'k'))/dab;
        stats.two.(['m_' sx '_ab_n']) = sum(strcmp(r2.resp_key(ab & msk),'none'))/dab;
        stats.two.(['m_' sx '_an_j']) = sum(strcmp(r2.resp_key(an & msk),'j'))/dan;
        stats.two.(['m_' sx '_an_k']) = sum(strcmp(r2.resp_key(an & msk),'k'))/dan;
        stats.two.(['m_' sx '_an_n']) = sum(strcmp(r2.resp_key(an & msk),'none'))/dan;
    end
    
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

%% visualization

% --- fig 1: 1-back ---
figure('color','w','Position',[100 100 1200 400]);
subplot(1,3,1); 
data = [get_v('one','acc_same'); get_v('one','acc_sim'); get_v('one','acc_new')]';
raincloud(data, {c_same, c_sim, c_new}, {'Hit(Same)','Hit(Sim)','CR(New)'}, 'Accuracy', '', [0,1.05]);
subplot(1,3,2); 
data = [get_v('one','rt_same'); get_v('one','rt_sim')]';
raincloud(data, {c_same, c_sim}, {'Hit(Same)','Hit(Sim)'}, 'RT (s)', '');
subplot(1,3,3); hold on;
mat = [mean(get_v('one','acc_same')), mean(get_v('one','err_same_as_sim')), mean(get_v('one','err_same_as_new'));
       mean(get_v('one','err_sim_as_same')), mean(get_v('one','acc_sim')), mean(get_v('one','err_sim_as_new'));
       mean(get_v('one','err_new_as_same')), mean(get_v('one','err_new_as_sim')), mean(get_v('one','acc_new'))];
draw_matrix(mat, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
sgtitle('1-Back Task Performance', 'FontSize', 16);

% --- fig 2: 2-back
figure('color','w','Position',[50 50 1200 900]);

% row 1: metrics
subplot(3,3,1);
data = [get_v('two','ldi_comp'); get_v('two','ldi_iso'); get_v('two','ldi_nov')]';
raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'LDI', 'Lure Discrimination');
subplot(3,3,4);
data = [get_v('two','d_AA_comp'); get_v('two','d_AA_iso'); get_v('two','d_AA_nov')]';
raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'd''', 'Recognition');
subplot(3,3,2);
data = [get_v('two','rt_AB_comp'); get_v('two','rt_AB_iso'); get_v('two','rt_AB_nov')]';
raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Lure Discrimination)');

% row 2: rts
subplot(3,3,5);
data = [get_v('two','rt_AA_comp'); get_v('two','rt_AA_iso'); get_v('two','rt_AA_nov')]';
raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Recognition)');


% row 3: 3 confusion matrices
lbl_row = {'Exp Same(A-A)','Exp Sim(A-B)','Exp New(A-N)'}; lbl_col = {'Resp Same','Resp Sim','Resp New'};
mat_colors = {c_same, c_sim, c_new};

subplot(3,3,7); 
m_c = [mean(get_v('two','m_comp_aa_j')), mean(get_v('two','m_comp_aa_k')), mean(get_v('two','m_comp_aa_n'));
       mean(get_v('two','m_comp_ab_j')), mean(get_v('two','m_comp_ab_k')), mean(get_v('two','m_comp_ab_n'));
       mean(get_v('two','m_comp_an_j')), mean(get_v('two','m_comp_an_k')), mean(get_v('two','m_comp_an_n'))];
draw_matrix(m_c, mat_colors, lbl_row, lbl_col); title('Compared', 'FontSize',14);

subplot(3,3,8); 
m_i = [mean(get_v('two','m_iso_aa_j')), mean(get_v('two','m_iso_aa_k')), mean(get_v('two','m_iso_aa_n'));
       mean(get_v('two','m_iso_ab_j')), mean(get_v('two','m_iso_ab_k')), mean(get_v('two','m_iso_ab_n'));
       mean(get_v('two','m_iso_an_j')), mean(get_v('two','m_iso_an_k')), mean(get_v('two','m_iso_an_n'))];
draw_matrix(m_i, mat_colors, lbl_row, lbl_col); title('Isolated', 'FontSize',14);

subplot(3,3,9); 
m_n = [mean(get_v('two','m_nov_aa_j')), mean(get_v('two','m_nov_aa_k')), mean(get_v('two','m_nov_aa_n'));
       mean(get_v('two','m_nov_ab_j')), mean(get_v('two','m_nov_ab_k')), mean(get_v('two','m_nov_ab_n'));
       mean(get_v('two','m_nov_an_j')), mean(get_v('two','m_nov_an_k')), mean(get_v('two','m_nov_an_n'))];
draw_matrix(m_n, mat_colors, lbl_row, lbl_col); title('Novel', 'FontSize',14);
sgtitle('2-Back Task Performance', 'FontSize', 16);

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

function draw_matrix(mat, cols, ylbl, xlbl)
    hold on;
    for r=1:3
        for c=1:3
            v = mat(r,c); 
            t_c = (v>0.5)*[1 1 1] + (v<=0.5)*[0.2 0.2 0.2];           
            patch([c-0.48 c+0.48 c+0.48 c-0.48], [r-0.48 r-0.48 r+0.48 r+0.48], ...
                  cols{r}, 'EdgeColor','none','FaceAlpha',v);
            text(c, r, sprintf('%.2f',v), 'HorizontalAlignment','center', ...
                 'Color',t_c, 'FontWeight','bold', 'FontSize',12);
            if r==1 
                text(c, 0.4, xlbl{c}, 'HorizontalAlignment','center', ...
                     'Color',cols{c}, 'FontWeight','bold', 'FontSize',11); 
            end
        end
        text(0.4, r, ylbl{r}, 'HorizontalAlignment','right', ...
             'Color',cols{r}, 'FontWeight','bold', 'FontSize',11);
    end 
    axis ij equal; 
    xlim([0.2 3.8]); ylim([0.3 3.5]);
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
    hold off;
end