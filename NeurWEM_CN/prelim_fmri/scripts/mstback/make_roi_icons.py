"""Render ROI thumbnail icons (sagittal T1 slice, ROI in red) matching the house
style used in prelim_fmri/icons/. Run in the fMRIPrep singularity container on Torch.
Each ROI -> group/rois/roi_<Name>.png (736x736, white bg). Masks mirror extract_2x2.py."""
import numpy as np, nibabel as nib
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from nilearn import datasets, plotting
OUT='/scratch/zb2164/CBI/NeurWEM/analysis/group/rois'
tmpl=datasets.load_mni152_template(resolution=1)
cort=datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm'); cd=cort.maps.get_fdata()
ju  =datasets.fetch_atlas_juelich('maxprob-thr25-2mm')
jd  =ju.maps.get_fdata() if hasattr(ju.maps,'get_fdata') else nib.load(ju.maps).get_fdata()
jaff=ju.maps.affine if hasattr(ju.maps,'affine') else nib.load(ju.maps).affine
def cmask(names):
    m=np.zeros(cd.shape,bool)
    for i,l in enumerate(cort.labels):
        if any(n in l for n in names): m|=(cd==i)
    return nib.Nifti1Image(m.astype(np.uint8),cort.maps.affine)
def jmask(names):
    m=np.zeros(jd.shape,bool)
    for i,l in enumerate([str(x) for x in ju.labels]):
        if any(n in l for n in names): m|=(jd==i)
    return nib.Nifti1Image(m.astype(np.uint8),jaff)
# {output name : (atlas, labels)}   atlas: 'cort' = Harvard-Oxford cortical, 'ju' = Juelich
ROIS={
 'LOC':        ('cort',['Lateral Occipital Cortex']),
 'Perirhinal': ('cort',['Temporal Fusiform Cortex, anterior division','Parahippocampal Gyrus, anterior division']),
 'EarlyVisual':('ju',  ['V1 BA17','V2 BA18']),
}
red=ListedColormap(['#e8392f'])
for name,(atlas,labels) in ROIS.items():
    img=(jmask if atlas=='ju' else cmask)(labels); d=img.get_fdata()>0
    world=nib.affines.apply_affine(img.affine,np.argwhere(d))
    rx=float(np.median(world[world[:,0]>0,0])) if (world[:,0]>0).any() else float(world[:,0].mean())
    fig=plt.figure(figsize=(4.6,4.6),facecolor='white'); ax=fig.add_axes([0,0,1,1])
    plotting.plot_roi(img,bg_img=tmpl,display_mode='x',cut_coords=[rx],cmap=red,
                      black_bg=False,annotate=False,draw_cross=False,alpha=0.85,axes=ax,figure=fig)
    fig.savefig(f'{OUT}/roi_{name}.png',dpi=160,facecolor='white'); plt.close(fig)
    print('roi_%s.png  x=%.0f  vox=%d'%(name,rx,int(d.sum())),flush=True)
