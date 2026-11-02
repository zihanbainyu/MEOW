#!/usr/bin/env python3
import sys, os, numpy as np, nibabel as nib
from scipy.ndimage import label
from multiprocessing import Pool
G='/scratch/zb2164/CBI/NeurWEM/analysis/group'
H=2.0; E=0.5; DH=0.1; B=5000; SEED=0; CONN=np.ones((3,3,3),int)
subs=sys.argv[1:]; n=len(subs)
imgs=[nib.load(f'{G}/sl_rho_sub-{s}.nii.gz') for s in subs]
aff=imgs[0].affine; grid=imgs[0].shape
A=np.stack([im.get_fdata().astype(np.float64) for im in imgs],0)
full=np.isfinite(A).sum(0)==n
idx=np.where(full)
X=A[:,idx[0],idx[1],idx[2]]
V=X.shape[1]
rt=np.sqrt(n)
def tstat(M):
    mean=M.mean(0); sd=M.std(0,ddof=1); t=np.zeros(V); ok=sd>0; t[ok]=mean[ok]*rt/sd[ok]; return t
def tfce_pos(s):
    out=np.zeros(V); smax=s.max()
    if smax<=0: return out
    vol=np.zeros(grid,bool)
    for h in np.arange(DH,smax+DH,DH):
        vol[:]=False; vol[idx]=s>=h
        lab,nl=label(vol,structure=CONN)
        if nl==0: continue
        sizes=np.bincount(lab.ravel()); sizes[0]=0
        lv=lab[idx]; ext=np.where(lv>0,sizes[lv],0)
        out+=(ext.astype(np.float64)**E)*(h**H)*DH
    return out
def perm(sg):
    t=tstat(sg[:,None]*X)
    p=tfce_pos(np.clip(t,0,None)); q=tfce_pos(np.clip(-t,0,None))
    return max(p.max(),q.max())
if __name__=='__main__':
    rng=np.random.default_rng(SEED)
    signs=np.ones((B,n)); signs[1:]=rng.choice([1.0,-1.0],size=(B-1,n))
    t0=tstat(X); pos=tfce_pos(np.clip(t0,0,None)); neg=tfce_pos(np.clip(-t0,0,None))
    mag=pos+neg; signed=pos-neg
    nw=int(os.environ.get('SLURM_CPUS_PER_TASK','8'))
    with Pool(nw) as pl:
        nullmax=np.array(pl.map(perm,[signs[i] for i in range(B)],chunksize=8))
    fwe=np.ones(V)
    order=np.argsort(nullmax)
    srt=nullmax[order]
    rank=np.searchsorted(srt,mag,side='left')
    fwe=(B-rank+1)/(B+1)
    def save(vec,nm,fill=0.0):
        vol=np.full(grid,fill,np.float32); vol[idx]=vec; nib.save(nib.Nifti1Image(vol,aff),f'{G}/{nm}.nii.gz')
    save(signed,'group_tfce'); save(fwe,'group_tfce_fwe_p',fill=1.0)
    tvol=t0.copy(); sigmask=fwe<0.05
    tfwe=np.where(sigmask,t0,0.0); save(tfwe,'group_t_fwe05')
    save(sigmask.astype(np.float64),'group_tfce_fwe05_mask')
    print('N=%d  V=%d  B=%d  TFCE H=%.1f E=%.1f dh=%.2f'%(n,V,B,H,E,DH),flush=True)
    print('null max-TFCE: p95=%.1f max=%.1f'%(np.percentile(nullmax,95),nullmax.max()),flush=True)
    print('voxels FWE<.05 = %d (two-sided)'%int(sigmask.sum()),flush=True)
    print('of those t>0 = %d ; t<0 = %d'%(int((sigmask&(t0>0)).sum()),int((sigmask&(t0<0)).sum())),flush=True)
