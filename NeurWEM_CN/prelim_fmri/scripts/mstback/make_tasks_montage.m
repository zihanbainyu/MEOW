function make_tasks_montage()
% Task GLM slice montage, MATLAB house style. 6 contrasts x (3 sagittal, 3 coronal,
% 3 axial). High-contrast 1mm MNI anatomy (local-contrast enhanced), cold-hot t overlay
% (|t|>3, clipped +/-10). Helvetica, vector-friendly PNG/PDF.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); g=load(fullfile(root,'maps','hires','geom.mat'));
xi=double(g.xidx)+1; yi=double(g.yidx)+1; zi=double(g.zidx)+1;
tpl=double(niftiread(fullfile(root,'maps','tpl.nii.gz')));
cons={'group_act_1back_t','group_act_2back_t','group_act_mst_t','group_c_2v1_t','group_c_2vmst_t','group_c_1vmst_t'};
ttl={'1-back > baseline','2-back > baseline','MST > baseline','2-back > 1-back','2-back > MST','1-back > MST'};
thr=3; vmax=10; INK=[0.149 0.149 0.184];
collab={'x = -30','x = 0','x = 30','y = -70','y = -20','y = 20','z = -12','z = 10','z = 40'};
tw=0.085; g1=0.006; G=0.020; x0=0.125;
cx=zeros(1,9); cx(1)=x0;
for k=2:9
    step=tw+g1; if k==4||k==7, step=tw+G; end
    cx(k)=cx(k-1)+step;
end
rh=0.132; rg=0.013; ytop=0.915; ry=ytop-rh-(0:5)*(rh+rg);
f=figure('Color','w','Units','inches','Position',[1 1 15 11]);
for r=1:6
    T=double(niftiread(fullfile(root,'maps','hires',[cons{r} '.nii.gz'])));
    for c=1:9
        if c<=3, v='x'; i=xi(c); elseif c<=6, v='y'; i=yi(c-3); else, v='z'; i=zi(c-6); end
        [an,ts]=getslice(tpl,T,v,i);
        rgb=overlay(an,ts,thr,vmax);
        ax=axes(f,'Position',[cx(c) ry(r) tw rh]); imshow(rgb,'Parent',ax); axis(ax,'image','off');
        if r==1 && c>=4   % L / R on coronal + axial
            text(ax,0.06,0.07,'L','Units','normalized','Color','w','FontName','Helvetica','FontSize',9,'FontWeight','bold');
            text(ax,0.94,0.07,'R','Units','normalized','Color','w','FontName','Helvetica','FontSize',9,'FontWeight','bold','HorizontalAlignment','right');
        end
    end
    annotation(f,'textbox',[0.005 ry(r) 0.115 rh],'String',ttl{r},'EdgeColor','none', ...
        'HorizontalAlignment','left','VerticalAlignment','middle','FontName','Helvetica','FontSize',11,'FontWeight','bold','Color',INK);
end
for c=1:9
    annotation(f,'textbox',[cx(c) ry(6)-0.028 tw 0.025],'String',collab{c},'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',9.5,'Color',INK);
end
% colorbar (custom cold-hot with grey sub-threshold gap)
cbx=0.955; cba=axes(f,'Position',[cbx 0.30 0.012 0.40]);
vv=linspace(vmax,-vmax,256)'; strip=zeros(256,1,3);
cc=tcolor(vv,thr,vmax); strip(:,1,:)=reshape(cc,256,1,3);
image(cba,strip); set(cba,'XTick',[],'YDir','reverse');
set(cba,'YTick',linspace(1,256,5),'YTickLabel',{'10','5','0','-5','-10'});
cba.YAxisLocation='right'; set(cba,'FontName','Helvetica','FontSize',9,'YColor',INK,'XColor','none','TickLength',[0 0]);
ylabel(cba,'t','FontName','Helvetica','FontSize',10,'Color',INK,'Rotation',0);
exportgraphics(f,fullfile(root,'figures','mstback','fig_tasks_montage.pdf'),'ContentType','vector');
disp('wrote fig_tasks_montage.pdf');
end
function [an,ts]=getslice(tpl,T,v,i)
switch v
    case 'x', a=flipud(squeeze(tpl(i,:,:))'); b=flipud(squeeze(T(i,:,:))');
    case 'y', a=flipud(squeeze(tpl(:,i,:))'); b=flipud(squeeze(T(:,i,:))');
    case 'z', a=flipud(squeeze(tpl(:,:,i))'); b=flipud(squeeze(T(:,:,i))');
end
an=enhance(a); ts=b;
end
function an=enhance(sl)
sl=double(sl); m=sl>0; an=zeros(size(sl));
if ~any(m(:)), return; end
v=sort(sl(m)); lo=v(max(1,round(0.02*numel(v)))); hi=v(round(0.995*numel(v)));
an=(sl-lo)/max(hi-lo,eps); an=min(max(an,0),1); an=an.^0.82;
try, an=adapthisteq(an,'ClipLimit',0.006,'NumTiles',[8 8]); catch, end
an(~m)=0;
end
function rgb=overlay(an,ts,thr,vmax)
rgb=repmat(an,[1 1 3]); brain=an>0;
idx=find((abs(ts)>=thr)&brain); if isempty(idx), return; end
cc=tcolor(ts(idx),thr,vmax); n=numel(an);
for ch=1:3
    layer=rgb(:,:,ch); layer(idx)=cc(:,ch); rgb(:,:,ch)=layer;
end
end
function cc=tcolor(tv,thr,vmax)
tv=tv(:); cc=repmat([0.80 0.80 0.82],numel(tv),1);  % sub-threshold grey
red0=[0.79 0.12 0.14]; yel=[1.00 0.92 0.23]; blu0=[0.13 0.30 0.73]; cya=[0.27 0.79 1.00];
p=tv>=thr; fp=min((tv(p)-thr)/(vmax-thr),1); cc(p,:)=red0+fp.*(yel-red0);
n=tv<=-thr; fn=min((-tv(n)-thr)/(vmax-thr),1); cc(n,:)=blu0+fn.*(cya-blu0);
end
