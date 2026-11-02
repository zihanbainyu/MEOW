function make_mst_behav()
% MST recognition behaviour, house style: three measures (recognition, lure
% discrimination, novel rejection), per-subject spaghetti + group mean/SEM.
here=fileparts(mfilename('fullpath')); S=load(fullfile(here,'mst_behav.mat'));
Y=[S.hit(:) S.lure(:) S.newCR(:)];            % N x 3
lab={'Old recognition','Lure discrimination','Novel rejection'};
N=size(Y,1); x=1:3;
INK=[0.149 0.149 0.184]; GREY=[0.72 0.72 0.76];
f=figure('Color','w','Units','inches','Position',[1 1 4.4 4.2]);
ax=axes(f); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',11);
% per-subject connected lines
jit=(linspace(-1,1,N))*0.06;
for i=1:N
    plot(ax,x+jit(i),Y(i,:),'-','Color',[GREY 0.9],'LineWidth',0.6);
    scatter(ax,x+jit(i),Y(i,:),18,GREY,'filled','MarkerEdgeColor','none');
end
% group mean + SEM
m=mean(Y,1); se=std(Y,0,1)/sqrt(N);
for k=1:3
    plot(ax,[x(k)-0.22 x(k)+0.22],[m(k) m(k)],'-','Color',INK,'LineWidth',2.4);
    plot(ax,[x(k) x(k)],[m(k)-se(k) m(k)+se(k)],'-','Color',INK,'LineWidth',1.2);
end
ax.XTick=x; ax.XTickLabel=lab; ax.XTickLabelRotation=15;
xlim(ax,[0.4 3.6]); ylim(ax,[0 1]); ax.YTick=0:0.25:1;
ylabel(ax,'proportion correct','FontName','Helvetica','FontSize',12);
title(ax,'MST recognition behaviour (N = 7)','FontName','Helvetica', ...
    'FontWeight','normal','FontSize',13,'HorizontalAlignment','left');
set(get(ax,'Title'),'Units','normalized','Position',[0 1.02 0]);
box(ax,'off'); ax.TickLength=[0 0]; ax.XColor=INK; ax.YColor=INK; ax.LineWidth=0.8;
exportgraphics(f,fullfile(here,'fig_mst_behav.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(here,'fig_mst_behav.png'),'Resolution',300);
disp('wrote fig_mst_behav.pdf / .png');
end
