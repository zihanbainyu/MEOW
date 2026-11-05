function compose_brain()
% Compose BrainNet 2-colour 3-view into a PNG: equal-size frames, even spacing, Helvetica legends.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); V=fullfile(root,'figures','mstback','brainnet_views');
imgs={'c2_side.png','c2_bottom.png','c2_back.png'};
A=cellfun(@(x) trimwhite(imread(fullfile(V,x))), imgs, 'UniformOutput',false);
H=max(cellfun(@(x) size(x,1),A)); W=max(cellfun(@(x) size(x,2),A));
A=cellfun(@(x) padto(x,H,W), A, 'UniformOutput',false);        % all same size -> even tiles
f=figure('Color','w','Units','inches','Position',[1 1 13.5 5.6]);
tl=tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
for i=1:3, ax=nexttile(tl); imshow(A{i},'Parent',ax); end
lax=axes(f,'Position',[0 0 1 1],'Visible','off'); hold(lax,'on'); xlim(lax,[0 1]); ylim(lax,[0 1]);
NET={'Visual','SomMot','DorsAttn','VentAttn','Limbic','Control','Default'}; C=jet(7);
x0=0.17; dx=0.097; y=0.085;
for k=1:7
    scatter(lax,x0+(k-1)*dx,y,70,C(k,:),'filled','MarkerEdgeColor','none');
    text(lax,x0+(k-1)*dx+0.011,y,NET{k},'FontName','Helvetica','FontSize',11,'VerticalAlignment','middle');
end
RED=[0.76 0.21 0.18]; BLUE=[0.17 0.43 0.70]; y2=0.03;
plot(lax,[0.34 0.37],[y2 y2],'Color',RED,'LineWidth',3);
text(lax,0.38,y2,'increase (2 > 1)','FontName','Helvetica','FontSize',11,'VerticalAlignment','middle');
plot(lax,[0.53 0.56],[y2 y2],'Color',BLUE,'LineWidth',3);
text(lax,0.57,y2,'decrease (2 < 1)','FontName','Helvetica','FontSize',11,'VerticalAlignment','middle');
exportgraphics(f,fullfile(root,'figures','mstback','fig_fc_brain_2color.png'),'Resolution',300);
disp('wrote fig_fc_brain_2color.png');
end
function B=trimwhite(A)
g=min(A,[],3); m=g<250; r=find(any(m,2)); c=find(any(m,1));
if isempty(r), B=A; else B=A(r(1):r(end),c(1):c(end),:); end
end
function B=padto(A,H,W)
[h,w,~]=size(A); t=floor((H-h)/2); l=floor((W-w)/2);
B=255*ones(H,W,size(A,3),'uint8'); B(t+1:t+h, l+1:l+w, :)=A;
end
