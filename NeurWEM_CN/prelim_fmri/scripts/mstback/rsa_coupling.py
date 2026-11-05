import csv, os, numpy as np
from collections import defaultdict
from scipy import stats

rows=[r for r in csv.DictReader(open(os.path.join(os.path.dirname(__file__),'..','..','data','mstback','rsa_reinstate.csv'))) if r['cell']=='B2-B1_PC']
D=defaultdict(dict); CF=defaultdict(list)
for r in rows:
    try: v=float(r['sim'])
    except: v=np.nan
    D[r['subj']].setdefault(r['roi'],[]).append(v)
    if r['roi']=='HC': CF[r['subj']].append(int(r['correct']))
fz=lambda r: np.arctanh(max(min(r,0.999),-0.999))

for cx in ['Angular','vmPFC','DLPFC']:
    zc=[]; zi=[]
    for s in sorted(D):
        if 'HC' not in D[s] or cx not in D[s]: continue
        hc=np.array(D[s]['HC']); co=np.array(D[s][cx]); cf=np.array(CF[s])
        if len(hc)!=len(cf): continue
        m=~np.isnan(hc)&~np.isnan(co)
        for acc,st in [(1,zc),(0,zi)]:
            sel=m&(cf==acc)
            st.append(fz(np.corrcoef(hc[sel],co[sel])[0,1]) if sel.sum()>=5 and np.std(hc[sel])>0 and np.std(co[sel])>0 else np.nan)
    zc=np.array(zc); zi=np.array(zi); ok=~np.isnan(zc)&~np.isnan(zi)
    rc=np.tanh(np.nanmean(zc[ok])); ri=np.tanh(np.nanmean(zi[ok]))
    t,p=stats.ttest_rel(zc[ok],zi[ok]); tc,pc=stats.ttest_1samp(zc[ok],0)
    print(f"HC-{cx:8s} r_cor={rc:+.3f} r_inc={ri:+.3f} | cor>inc t({ok.sum()-1})={t:+.2f} p={p:.3f} | cor>0 t={tc:+.2f} p={pc:.3f} (N={ok.sum()})")
