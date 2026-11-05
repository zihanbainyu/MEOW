function matrix()
% 2-back - 1-back FC, equal-size network blocks (7x7 grid of parcel sub-matrices),
% blue-white-red (clim +/-5), black labels, significant network-pair blocks outlined bold.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
S=load(fullfile(root,'data','mstback','parcel_matrix.mat')); H=load(fullfile(root,'data','mstback','heatmap_data.mat'));
NET=cellstr(S.netnames); b=double(S.bounds(:)'); cm=bwr(256); P=H.pd;   % 7x7 network-pair p-values
f=figure('Color','w','Units','inches','Position',[1 1 7.6 6.8]);
tl=tiledlayout(f,7,7,'TileSpacing','none','Padding','compact');
for i=1:7
    for j=1:7
        ax=nexttile(tl); sub=S.TD(b(i)+1:b(i+1), b(j)+1:b(j+1));
        imagesc(ax,sub); colormap(ax,cm); clim(ax,[-5 5]);
        set(ax,'XTick',[],'YTick',[]); axis(ax,'tight'); hold(ax,'on');
        n1=size(sub,1); n2=size(sub,2);
        lw=0.5; if P(i,j)<0.05, lw=2.6; end           % bold outline = significant pair
        plot(ax,[0.5 n2+0.5 n2+0.5 0.5 0.5],[0.5 0.5 n1+0.5 n1+0.5 0.5],'k-','LineWidth',lw);
        if i==1, title(ax,NET{j},'FontName','Helvetica','FontWeight','normal','FontSize',9,'Color','k'); end
        if j==1, ylabel(ax,NET{i},'FontName','Helvetica','FontSize',9,'Color','k','Rotation',0, ...
                'HorizontalAlignment','right','VerticalAlignment','middle'); end
    end
end
cb=colorbar('peer',gca); cb.Layout.Tile='east'; cb.Label.String='t'; cb.FontName='Helvetica'; cb.Ticks=-5:2.5:5;
exportgraphics(f,fullfile(root,'figures','mstback','fig_fc_matrix.pdf'),'ContentType','vector');
disp('wrote fig_fc_matrix.pdf');
end
function c=bwr(n)
h=round(n/2);
b2w=[linspace(0.13,1,h)' linspace(0.21,1,h)' linspace(0.85,1,h)'];
w2r=[linspace(1,0.82,n-h)' linspace(1,0.15,n-h)' linspace(1,0.15,n-h)'];
c=[b2w;w2r];
end
