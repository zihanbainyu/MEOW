%% setup
base = '/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/viewer/searchlight';
c_same=[97 125 184]/255; c_sim=[255 191 205]/255; c_null=[150 150 150]/255;

%% searchlight map  (group t, TFCE sign-flip FWE<.05, N=9)
V = nv_loadnii(fullfile(base,'group_t_fwe05.nii.gz'));
nv_montage(V, [0 6]);

%% graded same / similar / null  (visual cortex)
T = readtable(fullfile(base,'visual_graded.csv'));
M = [T.same T.similar T.null];
nv_graded(M, {c_same,c_sim,c_null}, {'same','similar','null'});
fprintf('same %.4f | similar %.4f | null %.4f  (N=%d)\n', mean(M,1,'omitnan'), size(M,1));
fprintf('same>similar p=%.4f  similar>null p=%.4f\n', signrank(M(:,1),M(:,2)), signrank(M(:,2),M(:,3)));

%% functions
function V = nv_loadnii(f)
    if endsWith(f,'.gz'), t=tempname; g=gunzip(f,t); f=g{1}; end
    V = double(niftiread(f));
end

function nv_montage(V, clim)
    V(~isfinite(V))=0; z=squeeze(sum(sum(V~=0,1),2)); zk=find(z>0);
    sl=round(linspace(zk(1),zk(end),12));
    f=figure('Color','k','Position',[60 60 980 720]); t=tiledlayout(f,3,4,'Padding','none','TileSpacing','none');
    cmap=hot(256);
    for i=1:numel(sl)
        nexttile; im=rot90(V(:,:,sl(i))); imagesc(im,clim); axis image off; colormap(cmap);
    end
    cb=colorbar; cb.Color='w'; cb.Layout.Tile='east';
    title(t,sprintf('searchlight model-RSA  (group t, TFCE FWE<.05, N=9)'),'Color','w','FontSize',14);
end

function nv_graded(M, cols, labs)
    [n,k]=size(M); set(0,'DefaultAxesFontName','Helvetica');
    f=figure('Color','w','Position',[80 80 460 460]); ax=axes(f); hold(ax,'on');
    jit=-0.15-rand(n,k)*0.2; X=repmat(1:k,n,1)+jit;
    for s=1:n, plot(X(s,:),M(s,:),'-','Color',[0.7 0.7 0.7 0.4],'LineWidth',0.5); end
    for i=1:k
        d=M(:,i); d=d(~isnan(d)); [ff,xi]=ksdensity(d); ff=ff/max(ff)*0.4;
        patch([i+ff,i*ones(1,numel(ff))],[xi,fliplr(xi)],cols{i},'EdgeColor','none','FaceAlpha',0.5);
        scatter(X(:,i),M(:,i),26,cols{i},'filled','MarkerFaceAlpha',0.7);
        q=quantile(d,[.25 .5 .75]); rectangle('Position',[i-0.06 q(1) 0.12 q(3)-q(1)],'EdgeColor','k','LineWidth',1.2,'FaceColor',cols{i});
        plot([i-0.06 i+0.06],[q(2) q(2)],'k-','LineWidth',2);
    end
    yline(0,'k:'); set(ax,'XTick',1:k,'XTickLabel',labs,'FontSize',15,'TickDir','out','LineWidth',1);
    xlim([0.4 k+0.6]); ylabel('neural similarity (r)','FontSize',15); box off;
    base=max(M(:))+0.05*range(M(:));
    for p=1:k-1
        pv=signrank(M(:,p),M(:,p+1)); st=repmat('*',1,(pv<.05)+(pv<.01)+(pv<.001)); if isempty(st), st='n.s.'; end
        y=base+(p-1)*0.06*range(M(:));
        plot([p p p+1 p+1],[y-0.01*range(M(:)) y y y-0.01*range(M(:))],'k-','LineWidth',1);
        text(p+0.5,y,st,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);
    end
end
