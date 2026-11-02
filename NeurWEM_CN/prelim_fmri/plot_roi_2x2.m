T = readtable('/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/prelim_fmri/roi_2x2_persubj.csv');
FS = struct('tick',20,'lab',20,'ttl',20,'anno',20);
c_cor = [97 125 184]/255; c_inc = [196 80 80]/255;   % correct = house blue, incorrect = muted red
rois  = {'DLPFC','HC_whole'}; rlab = {'DLPFC','hippocampus'};
conds = {'compared','novel'};
subs  = unique(T.subj);
f = figure('Color','w','Position',[60 60 1120 860]);
p = 0;
for r = 1:numel(rois)
  for c = 1:numel(conds)
    p = p+1; subplot(2,2,p);
    M = nan(numel(subs),2);    % [incorrect correct]
    for i = 1:numel(subs)
      for a = [0 1]
        m = T.subj==subs(i) & strcmp(T.roi,rois{r}) & strcmp(T.condition,conds{c}) & T.acc==a;
        if any(m), M(i,2*(a==1)+ (a==0)) = T.z(m); end   % col1=incorrect(a0), col2=correct(a1)
      end
    end
    paired_plot(M, {'incorrect','correct'}, {c_inc,c_cor}, 'B2 activity (z)', ...
        sprintf('%s \x2014 %s', rlab{r}, conds{c}), 5, {[1 2]}, FS, subs);
    yline(0,'k:','LineWidth',0.6);
  end
end
exportgraphics(f, '/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/prelim_fmri/fig_dlpfc_vs_hc_2x2.pdf', 'ContentType','vector');
exportgraphics(f, '/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/prelim_fmri/fig_dlpfc_vs_hc_2x2.png', 'Resolution',200);
fprintf('wrote fig_dlpfc_vs_hc_2x2\n');

function paired_plot(M, lvllbl, cols, ylbl, ttl, min_n, pairs, FS, subj)
    [nS, nL] = size(M);
    hold on; jw = 0.07; bw = 0.34; dot_area = 90;
    X = (1:nL) + (rand(nS,nL)-0.5)*2*jw;
    s_col = lines(nS);
    mks = {'o','s','d','^','v','p','h','>','<'};
    s_mk = mks(mod(0:nS-1, numel(mks)) + 1);
    hs = gobjects(nS,1);
    for s = 1:nS
        hs(s) = plot(X(s,:), M(s,:), '-', 'Color', [s_col(s,:) 0.8], 'LineWidth', 1.4, ...
            'Marker', s_mk{s}, 'MarkerSize', 8, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor','none');
    end
    for L = 1:nL
        d = M(:,L); d = d(~isnan(d));
        if numel(d) >= 2
            q = quantile(d,[.25 .5 .75]); iqrv = q(3)-q(1);
            lo_w = max(min(d), q(1)-1.5*iqrv); hi_w = min(max(d), q(3)+1.5*iqrv);
            plot([L L],[lo_w q(1)],'k-','LineWidth',0.9);
            plot([L L],[q(3) hi_w],'k-','LineWidth',0.9);
            plot(L+[-.06 .06],[lo_w lo_w],'k-','LineWidth',0.9);
            plot(L+[-.06 .06],[hi_w hi_w],'k-','LineWidth',0.9);
            rectangle('Position',[L-bw/2, q(1), bw, max(iqrv,eps)], 'EdgeColor','k','LineWidth',1.1,'FaceColor','none');
            plot(L+[-bw/2 bw/2],[q(2) q(2)],'k-','LineWidth',1.8);
        end
        for s = 1:nS
            scatter(X(s,L), M(s,L), dot_area, cols{L}, s_mk{s}, 'filled', ...
                'MarkerFaceAlpha',0.85, 'MarkerEdgeColor',[.2 .2 .2], 'LineWidth',0.3);
        end
    end
    all_v = M(~isnan(M)); yr = range(all_v); if isempty(yr)||yr==0, yr = 1; end
    base = max(all_v) + 0.10*yr;
    if nS >= min_n
        for pp = 1:numel(pairs)
            ij = pairs{pp}; a = M(:,ij(1)); b = M(:,ij(2)); ok = ~isnan(a)&~isnan(b);
            pv = signrank(a(ok), b(ok));
            s = 'n.s.'; if pv<.001, s='***'; elseif pv<.01, s='**'; elseif pv<.05, s='*'; elseif pv<.1, s=sprintf('p=%.2f',pv); end
            y = base;
            plot([ij(1) ij(1) ij(2) ij(2)], [y-0.02*yr y y y-0.02*yr], 'k-','LineWidth',0.9);
            text(mean(ij), y, s, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FS.tick);
        end
    end
    set(gca,'XTick',1:nL,'XTickLabel',lvllbl,'FontSize',FS.tick);
    xlim([0.4 nL+0.6]); ylabel(ylbl,'FontSize',FS.lab);
    ylim([min(all_v)-0.10*yr, base+0.16*yr]);
    if ~isempty(ttl), title(ttl,'FontSize',FS.ttl); end
    box off;
end
