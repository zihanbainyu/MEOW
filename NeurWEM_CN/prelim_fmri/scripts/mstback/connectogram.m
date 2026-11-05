function connectogram()
% 2-colour FC connectogram (red = increase, blue = decrease)
here = fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
S = load(fullfile(root,'data','mstback','connectogram_data.mat'));   % edges(Mx3 i,j,t), net(1x200), ang(1x200), yeo(7x3), netnames
pos = [cos(S.ang(:)), sin(S.ang(:))];
RED=[0.76 0.21 0.18]; BLUE=[0.17 0.43 0.70];
f=figure('Color','w','Units','inches','Position',[1 1 8 8]); ax=axes(f); hold(ax,'on');
axis(ax,'equal','off'); set(ax,'FontName','Helvetica');
E=S.edges; [~,si]=sort(abs(E(:,3))); E=E(si,:); tv=max(abs(E(:,3)));
for e=1:size(E,1)
    i=E(e,1); j=E(e,2); t=E(e,3); p0=pos(i,:); p2=pos(j,:); c=0.12*(p0+p2);
    u=linspace(0,1,40)'; B=(1-u).^2*p0 + 2*(1-u).*u*c + u.^2*p2;
    col=RED; if t<0, col=BLUE; end
    lw=0.3+1.8*min(abs(t)/tv,1); al=0.30+0.5*min(abs(t)/tv,1);
    plot(ax,B(:,1),B(:,2),'Color',[col al],'LineWidth',lw);
end
for k=1:7
    idx=find(S.net==k); scatter(ax,pos(idx,1),pos(idx,2),16,S.yeo(k,:),'filled','MarkerEdgeColor','none');
end
for k=1:7
    idx=find(S.net==k); a=S.ang(idx); a0=min(a); a1=max(a); aa=linspace(a0,a1,40);
    plot(ax,1.08*cos(aa),1.08*sin(aa),'Color',S.yeo(k,:),'LineWidth',6);
    am=(a0+a1)/2; r=1.2; rot=rad2deg(am); if cos(am)<0, rot=rot+180; end
    text(ax,r*cos(am),r*sin(am),S.netnames{k},'Color',S.yeo(k,:),'FontName','Helvetica', ...
        'FontSize',12,'HorizontalAlignment','center','Rotation',rot-90);
end

plot(ax,[-0.28 -0.12],[0.06 0.06],'Color',RED,'LineWidth',3);
text(ax,-0.08,0.06,'increase (2 > 1)','FontName','Helvetica','FontSize',11,'VerticalAlignment','middle');
plot(ax,[-0.28 -0.12],[-0.06 -0.06],'Color',BLUE,'LineWidth',3);
text(ax,-0.08,-0.06,'decrease (2 < 1)','FontName','Helvetica','FontSize',11,'VerticalAlignment','middle');
xlim(ax,[-1.35 1.35]); ylim(ax,[-1.35 1.35]);
exportgraphics(f,fullfile(root,'figures','mstback','fig_fc_connectogram.pdf'),'ContentType','vector');
disp('wrote fig_fc_connectogram.pdf');
end
