function gaze_reinstat()
% NeurWEM gaze-reinstatement replication (paper-1 metric ported to the scanner data).
% Parses EyeLink .asc -> per-trial duration-weighted fixation-density maps (20x20 over
% the stimulus ROI, sigma=2, normalised), computes item-specific gaze similarity
% reinst_index = corr(map_X,map_Y) - mean corr to other items' Y maps, and tests:
%   (1) Fig-1 foundational: same-item vs similar-item reinstatement (>0, same>similar)
%   (2) mechanism indices A1-B1 (separation), B2-B1 (completion), A2-B1 (predictive)
here=fileparts(mfilename('fullpath')); root=fileparts(here);
ascdir=fullfile(root,'asc'); metadir=fullfile(root,'meta'); resdir=fullfile(root,'results');
subs=[101 102 103 105 107 108 109];
P=struct('roi_x',[760 1160],'roi_y',[340 740],'grid_x',20,'grid_y',20,'sigma',2);

% ---------- parse + build per-trial fixation maps, tagged with metadata ----------
T=table();                                   % one row per trial: subj,task,item,identity,condition,b2correct,map
for s=subs
    M=readtable(fullfile(metadir,sprintf('sub-%d_meta.csv',s)),'TextType','string');
    for tsk=1:2
        for blk=1:4
            f=fullfile(ascdir,sprintf('%d_%d_%d.asc',s,tsk,blk));
            if ~isfile(f), continue; end
            runstr=metarun(s,tsk,blk);
            trials=parse_asc(f);             % struct array: .trialid .stim .x .y .dur
            for k=1:numel(trials)
                tr=trials(k); row=M(M.run==runstr & M.order==(tr.trialid-1),:);
                if height(row)~=1, continue; end
                map=create_fixation_map(tr.x,tr.y,tr.dur,P);
                if ~any(map(:)), continue; end
                T=[T; cell2table({s,string(row.task),double(row.item),string(row.identity), ...
                    string(row.condition),string(row.stim_id),string(row.resp),string(row.resp_key),{map}}, ...
                    'VariableNames',{'subj','task','item','identity','condition','stim','resp','respkey','map'})]; %#ok<AGROW>
            end
        end
    end
    fprintf('parsed sub-%d (%d trials so far)\n',s,height(T));
end
save(fullfile(resdir,'gaze_maps.mat'),'T','P');

% ---------- B2 accuracy per item (from the 2-back B probe) ----------
isB2=T.task=="2_back" & T.identity=="B";
acc=containers.Map('KeyType','char','ValueType','double');
for i=find(isB2)'
    r=T.resp(i); k=T.respkey(i);
    if ismember(r,["1","2"]); acc(key(T.subj(i),T.item(i)))=double(strrep(k,'.0','')==r); end
end

% ---------- compute reinst_index for the key comparisons ----------
% roles: A1=1back/A, B1=1back/B, A2=2back/A, B2=2back/B (per subj,item)
rows=[];
for s=subs
    Ts=T(T.subj==s,:);
    getmap=@(tk,id,it) firstmap(Ts,tk,id,it);
    items=unique(Ts.item(Ts.task=="2_back" & Ts.identity=="B"));
    % baseline pools (other items' maps) per role
    for it=items'
        A1=getmap("1_back","A",it); B1=getmap("1_back","B",it);
        A2=getmap("2_back","A",it); B2=getmap("2_back","B",it);
        cond=firstval(Ts,"2_back","B",it,'condition');
        ak=key(s,it); a=NaN; if isKey(acc,ak), a=acc(ak); end
        % same-item
        r_a1a2=ri(A1,A2,poolmaps(Ts,"2_back","A",it));
        r_b1b2=ri(B1,B2,poolmaps(Ts,"2_back","B",it));
        % similar-item (pairmates)
        r_a1b1=ri(A1,B1,poolmaps(Ts,"1_back","B",it));   % separation
        r_b2b1=ri(B2,B1,poolmaps(Ts,"1_back","B",it));   % completion
        r_a2b1=ri(A2,B1,poolmaps(Ts,"1_back","B",it));   % predictive
        rows=[rows; s,it,double(cond=="compared"),a,r_a1a2,r_b1b2,r_a1b1,r_b2b1,r_a2b1]; %#ok<AGROW>
    end
end
R=array2table(rows,'VariableNames',{'subj','item','compared','b2correct','a1a2','b1b2','a1b1','b2b1','a2b1'});
writetable(R,fullfile(resdir,'gaze_reinstat_trials.csv'));

% ---------- (1) Fig-1: same-item vs similar-item ----------
sm=[]; for s=subs
    Rs=R(R.subj==s,:);
    same=nanmean([Rs.a1a2;Rs.b1b2]); sim=nanmean([Rs.a1b1;Rs.b2b1;Rs.a2b1]);
    sm=[sm; s,same,sim]; end %#ok<AGROW>
[~,ps]=ttest(sm(:,2)); [~,psi]=ttest(sm(:,3)); [~,pss]=ttest(sm(:,2),sm(:,3));
fprintf('\n=== Fig-1 reinstatement (N=%d) ===\n',size(sm,1));
fprintf('same-item   M=%.4f vs0 p=%.4f\n',nanmean(sm(:,2)),ps);
fprintf('similar-item M=%.4f vs0 p=%.4f\n',nanmean(sm(:,3)),psi);
fprintf('same vs similar p=%.4f\n',pss);

% ---------- (2) mechanism indices (overall, vs 0) ----------
for nm=["a1b1" "b2b1" "a2b1"]
    v=[]; for s=subs, v=[v; nanmean(R.(char(nm))(R.subj==s))]; end %#ok<AGROW>
    [~,p]=ttest(v); fprintf('%-5s M=%+.4f vs0 p=%.4f\n',nm,nanmean(v),p);
end
save(fullfile(resdir,'gaze_reinstat_summary.mat'),'R','sm');
disp('done gaze_reinstat');
end

% ================= helpers =================
function s=metarun(subj,tsk,blk)
if subj==102, s=string(sprintf('run-%02d',(blk-1)*2+tsk));
else, s=string(sprintf('task-%dback_run-%02d',(tsk==2)+1,blk)); end
end
function k=key(s,it), k=sprintf('%d_%d',s,it); end
function trials=parse_asc(f)
fid=fopen(f); trials=struct('trialid',{},'stim',{},'x',{},'y',{},'dur',{});
cur=0; st=-1; stim="";
while ~feof(fid)
    ln=fgetl(fid); if ~ischar(ln), break; end
    if startsWith(ln,'MSG')
        if contains(ln,'TRIALID')
            t=sscanf(ln,'MSG %*d TRIALID %d');
            if ~isempty(t)&&t>0, cur=t; st=-1; stim=""; trials(end+1)=struct('trialid',t,'stim','','x',[],'y',[],'dur',[]); end %#ok<AGROW>
        elseif contains(ln,'STIM_ONSET')&&cur>0
            c=textscan(ln,'%*s %f %*s %s'); st=c{1}; stim=string(c{2}{1}); trials(end).stim=stim;
        end
    elseif startsWith(ln,'EFIX')&&cur>0&&st>0
        d=sscanf(ln,'EFIX %*c %d %d %d %f %f %*d');  % start end dur x y
        if numel(d)>=5 && d(1)>=st
            trials(end).dur(end+1)=d(3); trials(end).x(end+1)=d(4); trials(end).y(end+1)=d(5);
        end
    end
end
fclose(fid);
end
function map=create_fixation_map(x,y,dur,p)
in=x>=p.roi_x(1)&x<=p.roi_x(2)&y>=p.roi_y(1)&y<=p.roi_y(2); x=x(in);y=y(in);dur=dur(in);
map=zeros(p.grid_y,p.grid_x); if isempty(x), return; end
xb=linspace(p.roi_x(1),p.roi_x(2),p.grid_x+1); yb=linspace(p.roi_y(1),p.roi_y(2),p.grid_y+1);
for i=1:numel(x)
    xi=find(x(i)>=xb(1:end-1)&x(i)<xb(2:end),1); yi=find(y(i)>=yb(1:end-1)&y(i)<yb(2:end),1);
    if ~isempty(xi)&&~isempty(yi), map(yi,xi)=map(yi,xi)+dur(i); end
end
if p.sigma>0, map=imgaussfilt(map,p.sigma); end
if sum(map(:))>0, map=map/sum(map(:)); end
end
function m=firstmap(Ts,tk,id,it)
r=find(Ts.task==tk&Ts.identity==id&Ts.item==it,1); if isempty(r), m=[]; else, m=Ts.map{r}; end
end
function v=firstval(Ts,tk,id,it,col)
r=find(Ts.task==tk&Ts.identity==id&Ts.item==it,1); if isempty(r), v=""; else, v=Ts.(col)(r); end
end
function pool=poolmaps(Ts,tk,id,it)
idx=find(Ts.task==tk&Ts.identity==id&Ts.item~=it); pool=Ts.map(idx);
end
function r=ri(mx,my,pool)
if isempty(mx)||isempty(my), r=NaN; return; end
match=corr(mx(:),my(:));
b=nan(numel(pool),1); for i=1:numel(pool), if ~isempty(pool{i}), b(i)=corr(mx(:),pool{i}(:)); end, end
r=match-nanmean(b);
end
