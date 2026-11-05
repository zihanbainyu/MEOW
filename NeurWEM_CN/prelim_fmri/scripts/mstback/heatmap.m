function heatmap()
% 7x7 network FC heatmaps (1-back, 2-back, diff-t), blue-grey-red diverging, Helvetica, vector PDF.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); S=load(fullfile(root,'data','mstback','heatmap_data.mat'));
NET=cellstr(S.netnames); cm=bgr(256);
f=figure('Color','w','Units','inches','Position',[1 1 13.5 4.8]);
tl=tiledlayout(f,1,3,'Padding','compact','TileSpacing','loose');
vm=max(abs([S.m1(:);S.m2(:)]));
tri(nexttile(tl),S.m1,'1-back',[-vm vm],cm,NET,[]);
a2=tri(nexttile(tl),S.m2,'2-back',[-vm vm],cm,NET,[]);
tv=max(abs(S.td(:)));
a3=tri(nexttile(tl),S.td,'2-back - 1-back (t)',[-tv tv],cm,NET,S.pd);
cb=colorbar(a2); cb.Label.String='connectivity (Fisher z)'; cb.Label.FontName='Helvetica'; cb.FontName='Helvetica';
cb2=colorbar(a3); cb2.Label.String='t'; cb2.Label.FontName='Helvetica'; cb2.FontName='Helvetica';
exportgraphics(f,fullfile(root,'figures','mstback','fig_fc_heatmap.pdf'),'ContentType','vector');
disp('wrote fig_fc_heatmap.pdf');
end
function ax=tri(ax,M,ttl,cl,cm,NET,P)
up=triu(true(7),1); M(up)=NaN;
imagesc(ax,M,'AlphaData',~isnan(M)); set(ax,'Color','w');
colormap(ax,cm); clim(ax,cl); axis(ax,'equal','tight'); set(ax,'YDir','reverse');
ax.XTick=1:7; ax.YTick=1:7; ax.XTickLabel=NET; ax.YTickLabel=NET; ax.XTickLabelRotation=35;
set(ax,'FontName','Helvetica','FontSize',9,'Box','off','TickLength',[0 0]);
ax.XColor=[0.15 0.15 0.18]; ax.YColor=[0.15 0.15 0.18];
title(ax,ttl,'FontName','Helvetica','FontWeight','normal','FontSize',11,'Interpreter','none');
if ~isempty(P)
    hold(ax,'on');
    for i=1:7
        for j=1:i
            if P(i,j)<0.05
                text(ax,j,i,'*','HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',13,'FontName','Helvetica','Color',[0.1 0.1 0.1]);
            end
        end
    end
end
end
function c=bgr(n)
h=round(n/2);
b2g=[linspace(0.13,0.5,h)' linspace(0.21,0.5,h)' linspace(0.85,0.5,h)'];
g2r=[linspace(0.5,0.82,n-h)' linspace(0.5,0.15,n-h)' linspace(0.5,0.15,n-h)'];
c=[b2g;g2r];
end
