function em_network_brain()
% A priori EM network on a real (left-hemisphere) brain: hippocampal source +
% cortical reinstatement regions, painted inside a translucent shell, with
% leader-line callouts and citations. BrainNet-style render.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
S=load(fullfile(root,'data','mstback','em_network_meshes.mat'));
INK=[0.149 0.149 0.184]; SHELL=[0.80 0.80 0.83];
HC  =[0.77 0.31 0.32]; PERI=[0.90 0.57 0.53]; ANG=[0.30 0.52 0.74];
VT  =[0.24 0.60 0.40]; LOC =[0.55 0.74 0.45];
R={ 'hc',HC; 'peri',PERI; 'ang',ANG; 'vt',VT; 'loc',LOC };
% callouts: key, display name, citation, screen offset [right up] (mm), alignment
L={ 'hc',  'Hippocampus',       'Marr 1971; Nakazawa et al. 2002',[-88  16],'right';
    'peri','Perirhinal cortex', 'Brown & Aggleton 2001',          [-78 -46],'right';
    'ang', 'Angular gyrus',     'Vilberg & Rugg 2008',            [-34  74],'right';
    'vt',  'Ventral temporal',  'Polyn et al. 2005',              [-16 -66],'center';
    'loc', 'Lateral occipital', 'Wheeler et al. 2000',            [ 58  40],'left' };
f=figure('Color','w','Units','inches','Position',[1 1 11.0 6.6]);
ax=axes(f,'Position',[0.0 0.0 1.0 1.0]); hold(ax,'on');
for k=1:size(R,1), roipatch(ax,S.([R{k,1} '_v']),S.([R{k,1} '_f']),R{k,2}); end
patch(ax,'Vertices',S.shell_v,'Faces',S.shell_f,'FaceColor',SHELL,'FaceAlpha',0.15,...
    'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.65,'SpecularStrength',0.06,'BackFaceLighting','unlit');
axis(ax,'equal','off','vis3d'); view(ax,[-90 10]); camzoom(ax,0.82);
camlight(ax,'headlight'); camlight(ax,-28,32); lighting(ax,'gouraud'); material(ax,'dull');
% screen (image-plane) basis from the camera
vd=camtarget(ax)-campos(ax); vd=vd/norm(vd); u0=camup(ax);
rt=cross(vd,u0); rt=rt/norm(rt); up=cross(rt,vd); up=up/norm(up);
col=containers.Map(R(:,1),R(:,2));
for k=1:size(L,1)
    c=S.(['cen_' L{k,1}]); c=c(:)';
    off=L{k,4}; Lp=c+off(1)*rt+off(2)*up;
    plot3(ax,[c(1) Lp(1)],[c(2) Lp(2)],[c(3) Lp(3)],'-','Color',INK,'LineWidth',0.75);
    cc=col(L{k,1});
    str=sprintf('\\color[rgb]{%.3f,%.3f,%.3f}\\bf%s\\rm\n\\color[rgb]{0.42,0.42,0.46}\\it%s',cc(1),cc(2),cc(3),L{k,2},L{k,3});
    text(ax,Lp(1),Lp(2),Lp(3),str,'Interpreter','tex','HorizontalAlignment',L{k,5},...
        'VerticalAlignment','middle','FontName','Helvetica','FontSize',13,'Clipping','off');
end
exportgraphics(f,fullfile(root,'figures','mstback','fig_em_network_lat.pdf'),'ContentType','image','Resolution',400);
close(f);
end

function roipatch(ax,v,f,c)
patch(ax,'Vertices',v,'Faces',f,'FaceColor',c,'FaceAlpha',0.98,'EdgeColor','none',...
    'FaceLighting','gouraud','AmbientStrength',0.5,'DiffuseStrength',0.75,'SpecularStrength',0.12,'BackFaceLighting','reverselit');
end
