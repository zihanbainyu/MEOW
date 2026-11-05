"""
Mechanistic + normative model of episodic-working-memory integration in mnemonic
discrimination (MST-back). A reduced complementary-learning-systems (CLS) circuit:
  cortical content layer  -> hippocampal autoassociator (sparse, pattern-separating)
  -> hetero-associative read-out back to cortex (reinstatement).
At probe, the hippocampus pattern-completes the encoded target and reinstates its
cortical content; this episodic estimate is integrated with a working-memory estimate
by precision weighting. Episodic memory is available only for previously-compared items.
The model emits the empirical signatures (reinstatement-by-accuracy, separation->completion,
HC-cortex coupling, compared>novel, HC-lesion abolition) and the normative integration weight.
"""
import os, numpy as np
from scipy import stats
rng=np.random.default_rng(3)
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.join(HERE,'..','..'); FIG=os.path.join(ROOT,'figures','mstback')
Nc,Nh,kh = 120, 400, 40                       # cortical units, hippocampal units, HC sparsity (k active)
Bk, kE   = 20, 8.0                            # encoded item-bank size (interference); EM precision gain
def kwta(x,k):
    b=np.zeros_like(x); b[np.argpartition(x,-k)[-k:]]=1.0; return b
def unit(z): return z/(np.linalg.norm(z)+1e-12)
def cpat(): return unit(rng.standard_normal(Nc))
def nz(sig): return sig*rng.standard_normal(Nc)/np.sqrt(Nc)   # additive noise with norm ~ sig
def lure(T,rho):
    n=rng.standard_normal(Nc); n-=T*(T@n); return unit(rho*T+np.sqrt(1-rho**2)*unit(n))
RHO={'l2':0.90,'l3':0.80}                     # lure similarity to target (l2 more similar = harder)

def run_subject(sig_enc, sig_wm, crit, lesion=False, n_each=60, seed=0):
    r_=np.random.default_rng(seed)
    global rng; rng=r_
    Wch=r_.standard_normal((Nh,Nc))/np.sqrt(Nc)
    def hc(p): return kwta(Wch@p, kh)
    # encode a bank of targets (episodic store): HC autoassociator + HC->cortex read-out, with interference
    bank=[cpat() for _ in range(Bk)]; H=np.array([hc(b) for b in bank])
    M=H.T@H; np.fill_diagonal(M,0.0); Whc=np.array(bank).T@H
    trials=[]
    for cond in ['compared','novel']:
        for lv in ['l2','l3']:
            for ptype in ['lure','ident']:
                for _ in range(n_each):
                    if cond=='compared':
                        T=bank[r_.integers(Bk)]                          # target is an encoded item
                    else:
                        T=cpat()                                         # novel: never encoded -> no episodic trace
                    hT=hc(T)
                    if cond=='compared' and not lesion:
                        mask=(r_.random(Nc)<0.70).astype(float)          # partial, noisy retrieval cue
                        cue=unit(mask*T)+nz(sig_enc); hck=kwta(Wch@cue,kh)
                        for _ in range(3): hck=kwta(M@hck,kh)            # pattern completion vs competitors
                        r=unit(Whc@hck); fid=float(np.clip(r@T,0,1))     # cortical reinstatement fidelity (B2->B1)
                        comp=float(hck@hT/kh)                            # HC completion overlap with encoded code
                        piE=kE*fid**2/(1-fid**2+0.05); xE=r              # episodic precision ~ reinstatement reliability
                    else:
                        fid=comp=0.0; piE=0.0; xE=np.zeros(Nc)
                    sep=1.0-float(hc(lure(T,RHO[lv]))@hT/kh)             # HC pattern separation (lure vs target)
                    piW=1.0/sig_wm**2; xW=unit(T+nz(sig_wm))            # WM estimate: degraded over the 2-back delay
                    mem=unit((piW*xW+piE*xE)/(piW+piE))
                    probe = T if ptype=='ident' else lure(T,RHO[lv])
                    e = (1.0-float(probe@mem)) + 0.015*r_.standard_normal()   # change evidence
                    trials.append(dict(cond=cond,lv=lv,ptype=ptype,fid=fid,comp=comp,sep=sep,wE=piE/(piW+piE),e=e))
    # per-subject criterion: respond "similar" to the ~38% highest change-evidence trials (base-rate matched, bias-free)
    thr=np.quantile([t['e'] for t in trials],0.62)
    for t in trials:
        t['resp']=int(t['e']>thr); t['correct']=int((t['ptype']=='lure')==bool(t['resp']))
    return trials

def dprime(tr,cond):
    h=np.mean([t['resp'] for t in tr if t['cond']==cond and t['ptype']=='lure'])        # hit = call lure "similar"
    f=np.mean([t['resp'] for t in tr if t['cond']==cond and t['ptype']=='ident'])        # FA = call identical "similar"
    z=lambda p: stats.norm.ppf(min(max(p,1/160),1-1/160)); return z(h)-z(f)

# ---- population of model subjects (varying encoding/WM noise) ----
Nsub=24; subs=[]
for s in range(Nsub):
    se=0.9+0.5*rng.random(); sw=0.9+0.5*rng.random()
    subs.append(dict(se=se,sw=sw,crit=0.12+0.03*rng.standard_normal(),seed=100+s))
intact=[run_subject(d['se'],d['sw'],d['crit'],lesion=False,seed=d['seed']) for d in subs]
lesio =[run_subject(d['se'],d['sw'],d['crit'],lesion=True ,seed=d['seed']) for d in subs]

print("=== CLS integration model (N=%d simulated subjects) ==="%Nsub)
dc=np.array([dprime(tr,'compared') for tr in intact]); dn=np.array([dprime(tr,'novel') for tr in intact])
t,p=stats.ttest_rel(dc,dn,alternative='greater')
print("1) discrimination d'  compared=%.2f  novel=%.2f  (compared>novel t=%.2f p=%.4f)"%(dc.mean(),dn.mean(),t,p))
# reinstatement by accuracy (compared)
fc=np.array([np.mean([t['fid'] for t in tr if t['cond']=='compared' and t['correct']==1]) for tr in intact])
fi=np.array([np.mean([t['fid'] for t in tr if t['cond']=='compared' and t['correct']==0]) for tr in intact])
t,p=stats.ttest_rel(fc,fi,alternative='greater')
print("2) reinstatement fidelity  correct=%.3f  incorrect=%.3f  (correct>incorrect t=%.2f p=%.4f)"%(np.nanmean(fc),np.nanmean(fi),t,p))
# separation -> completion (within compared trials, across the pool)
allc=[t for tr in intact for t in tr if t['cond']=='compared']
r,p=stats.pearsonr([t['sep'] for t in allc],[t['comp'] for t in allc])
print("3) pattern separation -> completion   r=%.2f p=%.1e"%(r,p))
# HC completion <-> cortical reinstatement coupling (per subject, trial-wise), then group vs 0
cp=[]
for tr in intact:
    x=[t['comp'] for t in tr if t['cond']=='compared']; y=[t['fid'] for t in tr if t['cond']=='compared']
    if np.std(x)>0 and np.std(y)>0: cp.append(np.arctanh(np.clip(np.corrcoef(x,y)[0,1],-.999,.999)))
cp=np.array(cp); t,p=stats.ttest_1samp(cp,0)
print("4) HC-cortex reinstatement coupling   r=%.2f  (vs0 t=%.2f p=%.1e, N=%d)"%(np.tanh(cp.mean()),t,p,len(cp)))
# HC lesion abolishes the compared benefit
dcl=np.array([dprime(tr,'compared') for tr in lesio]); benefit_i=dc-dn; benefit_l=dcl-np.array([dprime(tr,'novel') for tr in lesio])
t,p=stats.ttest_rel(benefit_i,benefit_l)
print("5) HC lesion: compared-benefit intact=%.2f  lesioned=%.2f  (t=%.2f p=%.4f)"%(benefit_i.mean(),benefit_l.mean(),t,p))
# normative: episodic weight tracks reinstatement fidelity
wE=np.array([np.mean([t['wE'] for t in tr if t['cond']=='compared']) for tr in intact])
fd=np.array([np.mean([t['fid'] for t in tr if t['cond']=='compared']) for tr in intact])
r,p=stats.pearsonr(fd,wE)
print("6) normative integration: EM weight w_EM vs reinstatement fidelity  r=%.2f p=%.1e"%(r,p))

# ---- figures (house style) ----
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Helvetica','Arial'],'pdf.fonttype':42,'svg.fonttype':'none'})
INK='#26262f'; PUR='#948ac7'; GRY='#9a9aa2'; TARG='#c44e52'; FOIL='#5b8fc9'
def box(ax,x,v,c):
    v=np.array(v); v=v[np.isfinite(v)]; q=np.percentile(v,[25,50,75])
    ax.add_patch(plt.Rectangle((x-0.18,q[0]),0.36,q[2]-q[0],fc=c,ec=INK,alpha=0.25,lw=1))
    ax.plot([x-0.18,x+0.18],[q[1]]*2,color=INK,lw=1.8)
    ax.scatter(x+(rng.random(len(v))-0.5)*0.22,v,18,c,alpha=0.8,edgecolors='none')
fig,axs=plt.subplots(1,3,figsize=(9.6,3.4))
box(axs[0],1,dc,PUR); box(axs[0],2,dn,GRY); axs[0].set_xticks([1,2]); axs[0].set_xticklabels(['compared','novel'])
axs[0].set_ylabel("discrimination d'"); axs[0].set_title('1 | EM availability',fontsize=11,color=INK)
box(axs[1],1,fc,TARG); box(axs[1],2,fi,FOIL); axs[1].set_xticks([1,2]); axs[1].set_xticklabels(['correct','incorrect'])
axs[1].set_ylabel('reinstatement fidelity'); axs[1].set_title('2 | reinstatement x accuracy',fontsize=11,color=INK)
box(axs[2],1,benefit_i,PUR); box(axs[2],2,benefit_l,GRY); axs[2].set_xticks([1,2]); axs[2].set_xticklabels(['intact','HC lesion'])
axs[2].set_ylabel("compared benefit (Δd')"); axs[2].set_title('5 | causal HC test',fontsize=11,color=INK)
for a in axs: a.spines[['top','right']].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(FIG,'fig_cls_signatures.pdf')); plt.close(fig)
fig,ax=plt.subplots(figsize=(3.8,3.6))
ax.scatter(fd,wE,46,color=TARG,edgecolors='white'); z=np.polyfit(fd,wE,1); xx=np.linspace(fd.min(),fd.max(),10)
ax.plot(xx,np.polyval(z,xx),color=INK,lw=1.5); ax.set_xlabel('reinstatement fidelity (neural)'); ax.set_ylabel('episodic weight  w_EM')
ax.set_title('normative integration',fontsize=11,color=INK); ax.spines[['top','right']].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(FIG,'fig_cls_weighting.pdf')); plt.close(fig)
print("\nwrote fig_cls_signatures / fig_cls_weighting")
