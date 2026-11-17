function sl_housefig()
% House-style searchlight figure: anatomy + hot overlay, single colorbar, Helvetica, minimal.
base = '/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/viewer/searchlight';
A = ld(fullfile(base,'tpl_on_grid.nii.gz'));
T = ld(fullfile(base,'group_t_fwe05.nii.gz')); T(~isfinite(T))=0;
INK=[0.149 0.149 0.184]; VMAX=6; THR=1e-6; cmap=hot(256);
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');

% axial Z=-8 (k=34) and sagittal X=-14 (i=31)
axial = overlay(rot90(squeeze(A(:,:,34))), rot90(squeeze(T(:,:,34))), THR,VMAX,cmap);
sag   = overlay(rot90(squeeze(A(31,:,:))), rot90(squeeze(T(31,:,:))), THR,VMAX,cmap);

figW=7.2; figH=3.3; bot=0.08; hH=0.84; inchH=hH*figH;
aw = inchH*(size(axial,2)/size(axial,1))/figW;
sw = inchH*(size(sag,2)/size(sag,1))/figW;
f=figure('Color','w','Units','inches','Position',[1 1 figW figH]);

ax1=axes(f,'Position',[0.015 bot aw hH]); image(ax1,axial); axis(ax1,'image','off');
lab(ax1,size(axial,2),'Z = −8',true,INK);
ax2=axes(f,'Position',[0.02+aw+0.02 bot sw hH]); image(ax2,sag); axis(ax2,'image','off');
lab(ax2,size(sag,2),'X = −14',false,INK);

cx=0.015+aw+0.02+sw+0.05;
cax=axes(f,'Position',[cx 0.20 0.014 0.56]);
imagesc(cax,[0 1],[0 VMAX],reshape(cmap,[256 1 3])); set(cax,'YDir','normal');
cax.XTick=[]; cax.YAxisLocation='right'; cax.YTick=[0 VMAX]; cax.YColor=INK; cax.FontSize=12; cax.LineWidth=0.75; cax.TickLength=[0 0];
text(cax,0.5,1.10,'t','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12,'Color',INK);

exportgraphics(f,fullfile(base,'fig_searchlight_house.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(base,'fig_searchlight_house.png'),'Resolution',300);
close(f); disp('wrote fig_searchlight_house.pdf/png');
end

function lab(ax,ncol,txt,lr,INK)
text(ax,0.5,1.015,sprintf(txt),'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',13,'Color',INK);
if lr
    text(ax,0.05,1.015,'L','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',13,'Color',INK);
    text(ax,0.95,1.015,'R','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',13,'Color',INK);
end
end

function rgb = overlay(an,st,thr,vmax,cmap)
an=double(an); an=(an-min(an(:)))/(max(an(:))-min(an(:))+eps);
g=repmat(an,[1 1 3]);
m=st>thr; v=min(max(st,0),vmax)/vmax; ci=max(1,round(v*255)+1);
R=reshape(cmap(ci,1),size(st)); G=reshape(cmap(ci,2),size(st)); Bl=reshape(cmap(ci,3),size(st));
rr=g(:,:,1); gg=g(:,:,2); bb=g(:,:,3);
rr(m)=R(m); gg(m)=G(m); bb(m)=Bl(m); rgb=cat(3,rr,gg,bb);
end

function V = ld(f)
if endsWith(f,'.gz'), t=tempname; g=gunzip(f,t); f=g{1}; end
V=double(niftiread(f));
end
