function em_coupling_brain()
% Anterior HC -> cortical reinstatement coupling on a real brain (BrainNet style).
% ROIs whose coupling is above chance: Angular, Perirhinal, Ventral temporal, LOC.
% Edge thickness proportional to the group coupling (Fisher z, vs 0). No title.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
S=load(fullfile(root,'data','mstback','em_network_meshes.mat'));
INK=[0.149 0.149 0.184]; SHELL=[0.80 0.80 0.83];
HC =[0.77 0.31 0.32]; ANG=[0.30 0.52 0.74]; PERI=[0.88 0.54 0.12]; VT=[0.24 0.60 0.40]; LOC=[0.49 0.36 0.65];
% targets: key, color, coupling r (Fisher z, vs 0)
T={'ang',ANG,0.200; 'peri',PERI,0.208; 'vt',VT,0.204; 'loc',LOC,0.217};
f=figure('Color','w','Units','inches','Position',[1 1 8.4 8.4]);
ax=axes(f,'Position',[0.0 0.10 1.0 0.88]); hold(ax,'on');
% painted ROIs
roipatch(ax,S.hc_v,S.hc_f,HC);
for k=1:size(T,1), roipatch(ax,S.([T{k,1} '_v']),S.([T{k,1} '_f']),T{k,2}); end
% edges from anterior HC to each target, width proportional to coupling
for k=1:size(T,1)
    c=S.(['cen_' T{k,1}]); h=S.cen_hc;
    plot3(ax,[h(1) c(1)],[h(2) c(2)],[h(3) c(3)],'-','Color',[0.17 0.17 0.21],'LineWidth',1.8+T{k,3}*15);
end
% translucent cortical shell
patch(ax,'Vertices',S.shell_v,'Faces',S.shell_f,'FaceColor',SHELL,'FaceAlpha',0.15,...
    'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.65,'SpecularStrength',0.06,'BackFaceLighting','unlit');
axis(ax,'equal','off','vis3d'); view(ax,[-90 10]); camzoom(ax,1.05);
camlight(ax,'headlight'); camlight(ax,-28,32); lighting(ax,'gouraud'); material(ax,'dull');
% legend (bottom)
lg={'ant HC',HC; 'Angular',ANG; 'Perirhinal',PERI; 'vent temp',VT; 'LOC',LOC};
x0=0.035; w=0.152;
for k=1:size(lg,1)
    xx=x0+(k-1)*w;
    annotation(f,'rectangle',[xx 0.046 0.022 0.028],'FaceColor',lg{k,2},'EdgeColor','none');
    annotation(f,'textbox',[xx+0.028 0.033 w-0.03 0.052],'String',lg{k,1},'EdgeColor','none',...
        'FontName','Helvetica','FontSize',13,'Color',INK,'VerticalAlignment','middle','FitBoxToText','off');
end
annotation(f,'line',[0.800 0.840],[0.060 0.060],'Color',[0.17 0.17 0.21],'LineWidth',5);
annotation(f,'textbox',[0.848 0.033 0.15 0.052],'String','edge \propto r','EdgeColor','none',...
    'FontName','Helvetica','FontSize',13,'Color',INK,'VerticalAlignment','middle','HorizontalAlignment','left','FitBoxToText','off');
exportgraphics(f,fullfile(root,'figures','mstback','fig_coupling_brain.pdf'),'ContentType','image','Resolution',400);
close(f);
end

function roipatch(ax,v,fc,c)
patch(ax,'Vertices',v,'Faces',fc,'FaceColor',c,'FaceAlpha',0.98,'EdgeColor','none',...
    'FaceLighting','gouraud','AmbientStrength',0.5,'DiffuseStrength',0.75,'SpecularStrength',0.12,'BackFaceLighting','reverselit');
end
