function make_contrasts()
% One PDF per univariate map (no-fix GLM). 3 views (sag/cor/ax) x 3 cuts, every
% slice padded to a common square cell so panel sizes are consistent. MATLAB house
% style, high-contrast 1mm anatomy, cold-hot t overlay (or solid for binary). PDF only.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
maps={ 'v_1back','1-back > rest','t';
       'v_2back','2-back > rest','t';
       'v_mst','MST > rest','t';
       'c_2v1','2-back > 1-back','t';
       'c_2vmst','2-back > MST','t';
       'c_1vmst','1-back > MST','t';
       'union','Task-active (union, t>3 any)','bin';
       'conjunction','Common to all three (t>3)','bin' };
for i=1:size(maps,1)
    render(here,maps{i,1},maps{i,2},maps{i,3});
end
disp('wrote all contrast PDFs');
end

function render(here,mapname,ttl,kind)
root=fullfile(here,'..','..');
g=load(fullfile(root,'maps','hires','geom.mat')); xi=double(g.xidx)+1; yi=double(g.yidx)+1; zi=double(g.zidx)+1;
tpl=double(niftiread(fullfile(root,'maps','tpl.nii.gz')));
T=double(niftiread(fullfile(root,'maps','hires',[mapname '.nii.gz'])));
thr=3; vmax=10; INK=[0.149 0.149 0.184];
views={'x',xi,{'x = -30','x = 0','x = 30'},'sagittal';
       'y',yi,{'y = -70','y = -20','y = 20'},'coronal';
       'z',zi,{'z = -12','z = 10','z = 40'},'axial'};
CELL=233; GP=10; Wc=3*CELL+2*GP; Hc=3*CELL+2*GP;      % one square canvas
canvas=zeros(Hc,Wc,3);
for r=1:3
    v=views{r,1}; idx=views{r,2};
    for c=1:3
        [an,ts]=getslice(tpl,T,v,idx(c));
        if strcmp(kind,'bin'), rgb=overlay_bin(an,ts); else, rgb=overlay(an,ts,thr,vmax); end
        rgb=padcell(rgb,CELL);
        y0=(r-1)*(CELL+GP); x0=(c-1)*(CELL+GP); canvas(y0+1:y0+CELL, x0+1:x0+CELL, :)=rgb;
    end
end
f=figure('Color','w','Units','inches','Position',[1 1 7.6 7.8]);
axL=0.11; axB=0.075; axW=0.74; axH=axW*7.6/7.8;        % square axes (inches)
ax=axes(f,'Position',[axL axB axW axH]); imshow(canvas,'Parent',ax); axis(ax,'image','off'); hold(ax,'on');
fx=@(px) axL+(px/Wc)*axW; fy=@(py) axB+(1-py/Hc)*axH;
for r=1:3
    cy=(r-1)*(CELL+GP)+CELL/2;
    annotation(f,'textbox',[0.005 fy(cy)-0.03 0.095 0.06],'String',views{r,4},'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',10,'Color',INK);
    for c=1:3
        clab=views{r,3};
        text(ax,(c-1)*(CELL+GP)+10,(r-1)*(CELL+GP)+18,clab{c},'Color','w','FontName','Helvetica','FontSize',8.5);
        if r>=2
            text(ax,(c-1)*(CELL+GP)+14,(r-1)*(CELL+GP)+CELL-12,'L','Color','w','FontName','Helvetica','FontSize',8,'FontWeight','bold');
            text(ax,(c-1)*(CELL+GP)+CELL-14,(r-1)*(CELL+GP)+CELL-12,'R','Color','w','FontName','Helvetica','FontSize',8,'FontWeight','bold','HorizontalAlignment','right');
        end
    end
end
annotation(f,'textbox',[0 0.93 1 0.05],'String',ttl,'EdgeColor','none','HorizontalAlignment','center', ...
    'FontName','Helvetica','FontSize',14,'Color',INK);
if strcmp(kind,'t')
    cba=axes(f,'Position',[0.895 0.33 0.016 0.34]); vv=linspace(vmax,-vmax,256)'; strip=zeros(256,1,3);
    strip(:,1,:)=reshape(tcolor(vv,thr,vmax),256,1,3); image(cba,strip);
    set(cba,'XTick',[],'YDir','reverse','YTick',linspace(1,256,5),'YTickLabel',{'10','5','0','-5','-10'}, ...
        'YAxisLocation','right','FontName','Helvetica','FontSize',8.5,'YColor',INK,'TickLength',[0 0]); cba.XColor='none';
    ylabel(cba,'t','FontName','Helvetica','FontSize',10,'Color',INK,'Rotation',0);
end
exportgraphics(f,fullfile(root,'figures','mstback',['fig_' mapname '.pdf']),'ContentType','vector'); close(f);
fprintf('wrote fig_%s.pdf\n',mapname);
end

function B=padcell(A,N)
[h,w,~]=size(A); B=zeros(N,N,3); t=floor((N-h)/2); l=floor((N-w)/2);
B(t+1:t+h, l+1:l+w, :)=A;
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
rgb=repmat(an,[1 1 3]); brain=an>0; idx=find((abs(ts)>=thr)&brain); if isempty(idx), return; end
cc=tcolor(ts(idx),thr,vmax);
for ch=1:3, L=rgb(:,:,ch); L(idx)=cc(:,ch); rgb(:,:,ch)=L; end
end
function rgb=overlay_bin(an,ts)
rgb=repmat(an,[1 1 3]); brain=an>0; idx=find((ts>0.5)&brain); if isempty(idx), return; end
col=[0.79 0.12 0.14];
for ch=1:3, L=rgb(:,:,ch); L(idx)=col(ch); rgb(:,:,ch)=L; end
end
function cc=tcolor(tv,thr,vmax)
tv=tv(:); cc=repmat([0.80 0.80 0.82],numel(tv),1);
red0=[0.79 0.12 0.14]; yel=[1.00 0.92 0.23]; blu0=[0.13 0.30 0.73]; cya=[0.27 0.79 1.00];
p=tv>=thr; fp=min((tv(p)-thr)/(vmax-thr),1); cc(p,:)=red0+fp.*(yel-red0);
n=tv<=-thr; fn=min((-tv(n)-thr)/(vmax-thr),1); cc(n,:)=blu0+fn.*(cya-blu0);
end
