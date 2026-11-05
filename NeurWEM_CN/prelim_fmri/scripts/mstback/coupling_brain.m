function coupling_brain()
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
S=load(fullfile(root,'data','mstback','coupling_meshes.mat'));
INK=[0.149 0.149 0.184]; SHELL=[0.80 0.80 0.83];
HC=[0.768 0.306 0.321]; ANG=[0.192 0.431 0.702]; DLP=[0.882 0.541 0.122]; VT=[0.251 0.624 0.420]; VIS=[0.490 0.361 0.651];
C=com(S);
E={S.cen_ang,S.r_ang; S.cen_dlpfc,S.r_dlpfc; S.cen_vt,S.r_vt; S.cen_vis,S.r_vis};
V={-90 8; -115 10; -75 6};
for q=1:size(V,1)
    f=figure('Color','w','Units','inches','Position',[1 1 8.2 7.8]);
    ax=axes(f,'Position',[0.0 0.155 1.0 0.70]); hold(ax,'on');
    roipatch(ax,S.hc_v,S.hc_f,HC); roipatch(ax,S.ang_v,S.ang_f,ANG);
    roipatch(ax,S.dlpfc_v,S.dlpfc_f,DLP); roipatch(ax,S.vt_v,S.vt_f,VT); roipatch(ax,S.vis_v,S.vis_f,VIS);
    for k=1:size(E,1), edge(ax,C.hc,E{k,1},E{k,2},INK); end
    patch(ax,'Vertices',S.shell_v,'Faces',S.shell_f,'FaceColor',SHELL,'FaceAlpha',0.20,'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.6,'SpecularStrength',0.08,'BackFaceLighting','unlit');
    axis(ax,'equal','off','vis3d'); view(ax,[V{q,1} V{q,2}]); camzoom(ax,1.52);
    camlight(ax,'headlight'); camlight(ax,-28,32); lighting(ax,'gouraud'); material(ax,'dull');
    annotation(f,'textbox',[0.04 0.915 0.92 0.07],'String','Anterior hippocampus \rightarrow cortical reinstatement coupling','EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',18,'Color',INK);
    lg={'ant HC',HC,0.040;'Angular',ANG,0.165;'DLPFC',DLP,0.300;'vent temporal',VT,0.435;'visual',VIS,0.655};
    for k=1:size(lg,1)
        xx=lg{k,3};
        annotation(f,'rectangle',[xx 0.055 0.026 0.032],'FaceColor',lg{k,2},'EdgeColor','none');
        annotation(f,'textbox',[xx+0.033 0.040 0.18 0.052],'String',lg{k,1},'EdgeColor','none','FontName','Helvetica','FontSize',13,'Color',INK,'VerticalAlignment','middle','FitBoxToText','off');
    end
    annotation(f,'line',[0.800 0.835],[0.071 0.071],'Color',INK,'LineWidth',5);
    annotation(f,'textbox',[0.843 0.040 0.15 0.052],'String','edge \propto r','EdgeColor','none','FontName','Helvetica','FontSize',13,'Color',INK,'VerticalAlignment','middle','FitBoxToText','off');
    exportgraphics(f,fullfile(root,'figures','mstback',sprintf('fig_coupling_brain_v%d.pdf',q)),'ContentType','image','Resolution',320);
    close(f);
end
end

function roipatch(ax,v,fc,c)
patch(ax,'Vertices',v,'Faces',fc,'FaceColor',c,'FaceAlpha',1.0,'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.5,'DiffuseStrength',0.75,'SpecularStrength',0.15);
end

function edge(ax,a,b,r,INK)
plot3(ax,[a(1) b(1)],[a(2) b(2)],[a(3) b(3)],'-','Color',INK,'LineWidth',1.8+r*15);
end

function C=com(S)
C.hc=S.cen_hc; C.ang=S.cen_ang; C.dlpfc=S.cen_dlpfc; C.vt=S.cen_vt; C.vis=S.cen_vis;
end
