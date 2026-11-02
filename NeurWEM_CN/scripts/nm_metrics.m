function M = nm_metrics(Dc)
% Per-subject SDT + bias + NA/RT metrics for the 2-back similar-discrimination
% task, split by condition (compared/novel). Mirrors S_group_mstback scoring:
%   AA (goal A-A, corr '1') = same/target ; AB (goal A-B, corr '2') = lure/similar ;
%   AN (goal A-N, corr none) = foil. "similar" response = key '2'.
two=Dc.two;
g=string(two.goal); cond=string(two.condition); cr=string(two.corr_resp); rk=string(two.resp_key);
% paralleling the A-N selection used in the group script (A-N two trials after, same block)
pan=false(height(two),1); blk=two.block;
for i=1:height(two)-2
    if g(i)=="A-N" && blk(i+2)==blk(i), pan(i+2)=true; end
end
aa = g=="A-A" & cr=="1";
ab = g=="A-B" & cr=="2";
an = pan & cr=="none";
M=struct();
conds={'compared','novel'};
for c=1:2
    cm = cond==conds{c};
    % similar-discrimination: DI = p(sim|lure) - p(sim|foil)
    pL = pk(rk,ab&cm,'2'); pF = pk(rk,an&cm,'2');
    M.di(c)  = pL - pF;
    M.hitSim(c)=pL; M.faSim(c)=pF;
    % SDT on similar judgement
    nL=sum(ab&cm); nF=sum(an&cm);
    zH=zc(pL,nL); zF=zc(pF,nF);
    M.dSim(c) = zH - zF;                 % sensitivity for "similar"
    M.cSim(c) = -0.5*(zH + zF);          % criterion (bias toward "similar")
    % same-detection d' (as in group script): p(old|target) vs p(old|foil)
    pHo=pk(rk,aa&cm,'1'); pFo=pk(rk,an&cm,'1');
    M.dSame(c)= zc(pHo,sum(aa&cm)) - zc(pFo,sum(an&cm));
    % NA / miss by condition
    req = cm & (cr=="1"|cr=="2");
    M.miss(c) = mean(rk(req)=="none");
    M.na(c)   = mean(rk(cm)=="none");
    % median correct RT (press trials)
    rt=two.rt; good = cm & two.correct & (cr=="1"|cr=="2") & rt>0.15;
    M.rt(c)   = median(rt(good),'omitnan');
end
end
function p=pk(rk,mask,key), if sum(mask)==0, p=NaN; else, p=mean(rk(mask)==key); end, end
function z=zc(p,n), if isnan(p)||n==0, z=NaN; return; end, z=norminv(max(1/(2*n),min(1-1/(2*n),p))); end
