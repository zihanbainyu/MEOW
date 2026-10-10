%% NeurWEM viewer  —  open in MATLAB, run sections one at a time (Run Section / Ctrl+Enter)
% Self-serve viewing of the behavior, analysis results, figures, and group brain maps.
% Nothing here talks to the cluster: it reads local files. Refresh the results with the
% sync command printed by Section 0. Big single-trial betas stay on the cluster; Section 6
% loads one only if you have copied its betas_wb.mat locally.

%% 0. Setup  (edit paths if your repo lives elsewhere)
clear; clc;
REPO   = '/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN';
VIEWER = fullfile(REPO,'viewer');              % synced results live here
BEHAV  = fullfile(REPO,'data');                % sub<ID>/sub<ID>_concat.mat
SUBJS  = [101 102 103 104 105 106 107];        % add new subjects here
addpath(VIEWER);
fprintf('Results in: %s\n', VIEWER);
fprintf('To refresh results from the cluster, run in a terminal:\n');
fprintf(['  scp ''zb2164@login.torch.hpc.nyu.edu:/scratch/zb2164/CBI/NeurWEM/analysis/group/prelim/*.csv'' %s/data/\n'], VIEWER);
fprintf(['  scp ''zb2164@...:/scratch/zb2164/CBI/NeurWEM/analysis/group/mst/*.p*'' %s/figures/\n'], VIEWER);

%% 1. Behavior  —  MST-back discrimination + final MST accuracy  (tables)
fprintf('\n===== MST-back: B2 similar-discrimination accuracy (no-response = incorrect) =====\n');
fprintf('subj  compared(acc  n)   novel(acc  n)\n');
for s = SUBJS
    t = nv_concat(BEHAV,s,'results_2_back_all');
    if isempty(t), fprintf('%4d   (no 2-back table)\n',s); continue; end
    [cc,cn] = nv_b2acc(t,'compared'); [nc,nn] = nv_b2acc(t,'novel');
    fprintf('%4d   %.2f  (%2d)       %.2f  (%2d)\n', s, cc, cn, nc, nn);
end
fprintf('\n===== Final in-scanner MST: old / lure / new accuracy =====\n');
fprintf('subj   old    lure   new\n');
for s = SUBJS
    t = nv_concat(BEHAV,s,'results_mst');
    if isempty(t), fprintf('%4d   (no MST)\n',s); continue; end
    a = nv_mstacc(t);
    fprintf('%4d   %.2f   %.2f   %.2f\n', s, a.old, a.lure, a.new);
end

%% 2. Results  —  what analysis CSVs are available to plot
d = dir(fullfile(VIEWER,'data','*.csv'));
fprintf('\nAvailable result CSVs (use in Section 3):\n');
for i=1:numel(d), fprintf('  %s\n', d(i).name); end

%% 3. Plot a result   (set WHICH, then Run Section)
WHICH = 'b2_separation';     % one of: b2_separation | mst_roi | lss_vif_n6 | fig2_glm | fig3_rsa
nv_plot_result(VIEWER, WHICH);

%% 4. Open a saved figure (PDF/PNG)   (set WHICH, then Run Section)
d = dir(fullfile(VIEWER,'figures','*.*'));
fprintf('Figures available:\n'); for i=1:numel(d), if ~d(i).isdir, fprintf('  %s\n',d(i).name); end, end
WHICH = 'fig_pairmate_sep_3net.pdf';
f = fullfile(VIEWER,'figures',WHICH);
if ~exist(f,'file'), fprintf('not found: %s\n',f);
elseif usejava('desktop'), open(f);              % opens in MATLAB/Preview when interactive
else, fprintf('(would open %s)\n', f); end

%% 5. View a group brain map as axial slice montage   (set WHICH, then Run Section)
% maps/ names (drop .nii.gz): group_perm_t (searchlight), group_glm_interaction_t,
%   grp_old_net / grp_sim_net / grp_new_net / grp_recognition / grp_discrimination ...
WHICH = 'grp_old_net';
nv_showmap(fullfile(VIEWER,'maps',[WHICH '.nii.gz']), WHICH);

%% 6. Load one subject's single-trial betas  (needs betas_wb.mat copied locally)
% scp zb2164@...:/scratch/.../analysis/data/sub-105/betas_wb.mat  <VIEWER>/
SUBJ = 105;
f = fullfile(VIEWER,sprintf('betas_wb_sub-%d.mat',SUBJ));
if ~exist(f,'file'), f = fullfile(VIEWER,'betas_wb.mat'); end
if exist(f,'file')
    S = load(f);
    fprintf('betas: %d trials x %d voxels\n', size(S.betas,1), size(S.betas,2));
    tt = string(S.task); gl = string(S.goal); id = string(S.identity); cn = string(S.condition);
    fprintf('tasks: %s\n', strjoin(unique(tt),', '));
    b2 = tt=="2_back" & gl=="A-B" & id=="B";
    fprintf('B2 discrimination trials: %d (compared %d, novel %d)\n', sum(b2), sum(b2&cn=="compared"), sum(b2&cn=="novel"));
else
    fprintf('No local betas_wb for sub-%d. Copy it with the scp line above, then re-run.\n',SUBJ);
end

%% 7. Name the regions in a map  (load a *_regions.csv -> atlas labels + MNI)
% Glass brains can't localize; this table names each cluster (Harvard-Oxford) with its MNI peak.
d = dir(fullfile(VIEWER,'data','*_regions.csv'));
fprintf('Region tables available:\n'); for i=1:numel(d), fprintf('  %s\n', d(i).name); end
WHICH = 'grp_em_correct_vs_incorrect_regions.csv';
f = fullfile(VIEWER,'data',WHICH);
if exist(f,'file'), disp(readtable(f)); else, fprintf('not found: %s\n',f); end
% To label ANY map yourself, on the cluster run:
%   python3 analysis/code/label_map.py <map.nii.gz> <threshold> [minvox]
% e.g.  label_map.py analysis/group/group_perm_t.nii.gz 3.365 40

%% ================= helper functions =================
function t = nv_concat(base, s, field)
    t = [];
    f = fullfile(base, sprintf('sub%d',s), sprintf('sub%d_concat.mat',s));
    if ~exist(f,'file'), return; end
    S = load(f); if ~isfield(S,'final_data_output'), return; end
    fdo = S.final_data_output;
    if isfield(fdo,field) && ~isempty(fdo.(field)) && istable(fdo.(field)), t = fdo.(field); end
end
function [acc,n] = nv_b2acc(t, cond)
    iden=string(t.identity); goal=string(t.goal); c=string(t.condition);
    cr=string(t.corr_resp);  rk=string(t.resp_key);
    m = iden=="B" & goal=="A-B" & c==cond;     % similar-discrimination probe
    n = sum(m); acc = sum(rk(m)==cr(m)) / max(n,1);   % no-response (rk~=cr) counts wrong
end
function a = nv_mstacc(t)
    tt=string(t.trial_type); cr=string(t.corr_resp); rk=string(t.resp_key);
    miss = rk=="NA" | rk=="" | ismissing(rk);
    a.old  = mean(rk(tt=="old") ==cr(tt=="old"));
    a.lure = mean(rk(tt=="lure")==cr(tt=="lure"));
    a.new  = mean(miss(tt=="new"));            % new is correct when withheld
end
function nv_plot_result(VIEWER, which)
    T = readtable(fullfile(VIEWER,'data',[which '.csv']));
    cB=[0.16 0.44 0.74]; cR=[0.80 0.25 0.22];
    figure('Color','w');
    switch which
      case 'b2_separation'   % A2-B2 distance by condition x correctness, per ROI
        rois = unique(T.roi,'stable'); hold on; xt=[]; xl={};
        for r=1:numel(rois)
            x0=(r-1)*3;
            for c=0:1
                d=T.dist(strcmp(T.roi,rois{r})&T.correct==c);
                nv_box(x0+1+(c==0), d, ternary(c==1,cB,cR));
            end
            xt=[xt x0+1.5]; xl{end+1}=rois{r}; %#ok<AGROW>
        end
        set(gca,'XTick',xt,'XTickLabel',xl); ylabel('A2\leftrightarrowB2 distance (1-r)');
        title('Pairmate separation  (blue=correct, red=incorrect)');
      case 'mst_roi'
        cons=unique(T.contrast,'stable'); rois=unique(T.roi,'stable');
        for k=1:numel(cons)
            subplot(ceil(numel(cons)/2),2,k); hold on;
            for r=1:numel(rois)
                v=T.value(strcmp(T.contrast,cons{k})&strcmp(T.roi,rois{r}));
                bar(r,mean(v),0.6,'FaceColor',[.3 .45 .7],'FaceAlpha',.3);
                scatter(r+0*v,v,20,'k','filled','MarkerFaceAlpha',.5);
            end
            yline(0,'k-'); set(gca,'XTick',1:numel(rois),'XTickLabel',rois); title(cons{k},'Interpreter','none');
        end
      case 'lss_vif_n6'
        subs=unique(T.subj); hold on; xt=[]; xl={};
        for i=1:numel(subs)
            x0=(i-1)*3;
            nv_box(x0+1, log10(min(max(T.vif_lss(T.subj==subs(i)),1),1000)), cB);
            nv_box(x0+2, log10(min(max(T.vif_lsa(T.subj==subs(i)),1),1000)), cR);
            xt=[xt x0+1.5]; xl{end+1}=sprintf('s%d',subs(i)); %#ok<AGROW>
        end
        yline(log10(5),'k--'); set(gca,'XTick',xt,'XTickLabel',xl,'YTick',log10([1 5 30 1000]),'YTickLabel',{'1','5','30','1000'});
        ylabel('VIF'); title('LS-S (blue) vs LS-A (red) design collinearity');
      case 'fig2_glm'
        rois=unique(T.roi,'stable'); hold on; xt=[]; xl={};
        for r=1:numel(rois)
            x0=(r-1)*3;
            nv_box(x0+1,T.effect_pscb(strcmp(T.roi,rois{r})&strcmp(T.condition,'compared')),cB);
            nv_box(x0+2,T.effect_pscb(strcmp(T.roi,rois{r})&strcmp(T.condition,'novel')),cR);
            xt=[xt x0+1.5]; xl{end+1}=rois{r}; %#ok<AGROW>
        end
        yline(0,'k--'); set(gca,'XTick',xt,'XTickLabel',xl); ylabel('correct-incorrect (% signal)');
        title('GLM discrimination effect  (blue=compared, red=novel)');
      case 'fig3_rsa'
        mech=unique(T.mechanism,'stable');
        for k=1:numel(mech)
            subplot(1,numel(mech),k); hold on;
            subs=unique(T.subj); C=[];I=[];
            for s=subs'
                C(end+1)=mean(T.pair_r(T.subj==s&strcmp(T.mechanism,mech{k})&T.correct==1),'omitnan'); %#ok<AGROW>
                I(end+1)=mean(T.pair_r(T.subj==s&strcmp(T.mechanism,mech{k})&T.correct==0),'omitnan'); %#ok<AGROW>
            end
            nv_box(1,C',cB); nv_box(2,I',cR); yline(0,'k--');
            set(gca,'XTick',[1 2],'XTickLabel',{'correct','incorrect'}); title(mech{k},'Interpreter','none');
        end
      otherwise
        disp(head(T)); title('raw table head printed to console');
    end
end
function nv_showmap(f, name)
    V = nv_loadnii(f); if isempty(V), return; end
    V(~isfinite(V))=0;
    zc = round(linspace(0.30,0.80,12)*size(V,3));      % 12 axial slices
    lim = prctile(abs(V(V~=0)),99); if isempty(lim)||lim==0, lim=1; end
    figure('Color','w','Name',name);
    for i=1:numel(zc)
        subplot(3,4,i); imagesc(rot90(squeeze(V(:,:,zc(i)))),[-lim lim]); axis image off;
    end
    colormap(nv_coldhot); c=colorbar('Position',[0.93 0.3 0.015 0.4]); c.Label.String='value';
    sgtitle(sprintf('%s   (axial, clim \\pm%.2g)', name, lim),'Interpreter','none');
end
function V = nv_loadnii(f)
    V=[];
    if ~exist(f,'file'), fprintf('not found: %s\n',f); return; end
    try
        if endsWith(f,'.gz'), tmp=gunzip(f,tempdir); V=niftiread(tmp{1}); delete(tmp{1});
        else, V=niftiread(f); end
        V=double(V);
    catch e
        fprintf('Could not read NIfTI (need Image Processing Toolbox niftiread): %s\n', e.message);
    end
end
function nv_box(x,d,c)
    d=d(~isnan(d)); if isempty(d), return; end
    q=quantile(d,[.25 .5 .75]); w=0.3;
    patch([x-w x+w x+w x-w],[q(1) q(1) q(3) q(3)],c,'FaceAlpha',.15,'EdgeColor',c,'LineWidth',1.1); hold on;
    plot([x-w x+w],[q(2) q(2)],'-','Color',c,'LineWidth',2);
    scatter(x+(rand(size(d))-.5)*.2, d, 18, c,'filled','MarkerFaceAlpha',.6);
end
function cm = nv_coldhot()
    n=128;
    b=[linspace(0,1,n)' linspace(0,1,n)' ones(n,1)];     % blue -> white  (low)
    r=[ones(n,1) linspace(1,0,n)' linspace(1,0,n)'];     % white -> red   (high); 0 = white
    cm=[b; r];
end
function o = ternary(c,a,b), if c, o=a; else, o=b; end, end
