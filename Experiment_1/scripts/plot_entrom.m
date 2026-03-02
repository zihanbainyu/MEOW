%% Figure 1: Spatial entropy is NOT affected by condition (A and B items)
figure('Position', [100, 100, 900, 400]);

% A item entropy
subplot(1, 2, 1);
hold on;
data_a = [entropy_a_comp_subj, entropy_a_iso_subj];
raincloud(data_a, {[0.2, 0.6, 0.8], [0.8, 0.4, 0.4]}, {'Compared', 'Isolated'}, ...
    'Spatial Entropy (bits)', 'A Item Encoding', []);
title(sprintf('A Item: t(%d)=%.2f, p=%.3f', stats_entropy_a.df, stats_entropy_a.tstat, p_entropy_a));
ylim([5.8 6.8]);

% B item entropy
subplot(1, 2, 2);
hold on;
data_b = [entropy_b_comp_subj, entropy_b_iso_subj];
raincloud(data_b, {[0.2, 0.6, 0.8], [0.8, 0.4, 0.4]}, {'Compared', 'Isolated'}, ...
    'Spatial Entropy (bits)', 'B Item Encoding', []);
title(sprintf('B Item: t(%d)=%.2f, p=%.3f', stats_entropy_b.df, stats_entropy_b.tstat, p_entropy_b));
ylim([5.8 6.8]);

sgtitle('Individual Item Encoding: No Difference Between Conditions', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 2: A-B spatial overlap IS affected by condition (encoding only)
figure('Position', [100, 100, 500, 500]);
hold on;

% Compute subject-level means for encoding overlap
subj_overlap_comp = [];
subj_overlap_iso = [];
for sid = unique(ab_pairs_comp.subj_id)'
    idx_comp = ab_pairs_comp.subj_id == sid;
    idx_iso = ab_pairs_iso.subj_id == sid;
    
    overlap_comp_subj = overlap_comp(idx_comp);
    overlap_comp_subj = overlap_comp_subj(overlap_comp_subj ~= 0);
    
    overlap_iso_subj = overlap_iso(idx_iso);
    overlap_iso_subj = overlap_iso_subj(overlap_iso_subj ~= 0);
    
    if ~isempty(overlap_comp_subj), subj_overlap_comp = [subj_overlap_comp; mean(overlap_comp_subj)]; end
    if ~isempty(overlap_iso_subj), subj_overlap_iso = [subj_overlap_iso; mean(overlap_iso_subj)]; end
end

% Raincloud showing A-B overlap
data_overlap = [subj_overlap_comp, subj_overlap_iso];
raincloud(data_overlap, {[0.2, 0.6, 0.8], [0.8, 0.4, 0.4]}, {'Compared', 'Isolated'}, ...
    'A-B Spatial Overlap (r)', 'Relational Encoding (1-back)', []);
add_sig(data_overlap, [1 2]);

title(sprintf('Compared > Isolated: p < .001'), 'FontSize', 12);

sgtitle('Relational Encoding Differs Between Conditions', 'FontSize', 14, 'FontWeight', 'bold');

function raincloud(mat, cols, xlbls, ylbl, ttl, ylims)
    [n_rows, n_grps] = size(mat); hold on;
    d_v = mat(:); d_v = d_v(~isnan(d_v)); 
    mn = min(d_v); mx = max(d_v); rng = mx-mn; if rng==0, rng=1; end
    auto_lim = [mn-(rng*0.15), mx+(rng*0.15)];
    jit = -0.15 - (rand(size(mat)) * 0.20); x_c = repmat(1:n_grps, n_rows, 1) + jit;
    plot(x_c', mat', '-', 'Color', [0.7 0.7 0.7, 0.4], 'LineWidth', 0.5);
    for i = 1:n_grps
        d = mat(:,i); d = d(~isnan(d)); if isempty(d), continue; end
        [f, xi] = ksdensity(d); f = f/max(f)*0.4;
        patch([i+f, i*ones(1,length(f))], [xi, fliplr(xi)], cols{i}, 'EdgeColor','none','FaceAlpha',0.5);
        scatter(x_c(:,i), mat(:,i), 20, cols{i}, 'filled', 'MarkerFaceAlpha',0.6);
        q = quantile(d, [0.25, 0.5, 0.75]);
        rectangle('Position', [i-0.06, q(1), 0.12, q(3)-q(1)], 'FaceColor', cols{i}, 'EdgeColor','k', 'LineWidth',1.2);
        plot([i-0.06, i+0.06], [q(2), q(2)], 'k-', 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:n_grps, 'XTickLabel', xlbls, 'FontSize', 12);
    if ~isempty(ttl), title(ttl, 'FontSize', 14); end
    ylabel(ylbl,'FontSize',14,'FontWeight','bold'); xlim([0.2, n_grps+0.8]);
    if nargin>5 && ~isempty(ylims), ylim(ylims); else, ylim(auto_lim); end
    grid off; set(gca,'GridAlpha',0.1); box off; hold off;
end

function add_sig(data, pairs)
    [~, n_grps] = size(data);
    cl = ylim; y_top = cl(2); rng = cl(2)-cl(1); if rng==0, rng=1; end
    step = rng * 0.08; line_lvl = 0; hold on;
    real_mx = max(data(:)); if isnan(real_mx), real_mx=y_top; end
    base = max(y_top, real_mx + step*0.5);
    for i = 1:size(pairs,1)
        c1 = pairs(i,1); c2 = pairs(i,2);
        if c2==0, [~, p]=ttest(data(:,c1)); paired=false; else, [~, p]=ttest(data(:,c1),data(:,c2)); paired=true; end
        if p < 0.05
            line_lvl=line_lvl+1; y_p = base + (line_lvl*step);
            if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
            if paired, plot([c1,c1,c2,c2],[y_p-step*0.3,y_p,y_p,y_p-step*0.3],'k-','LineWidth',1.2);
               text(mean([c1 c2]), y_p+step*0.1, txt, 'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
            else, text(c1, y_p, txt, 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold'); end
        end
    end
    if line_lvl>0, ylim([cl(1), base+(line_lvl*step)+step]); end; hold off;
end
