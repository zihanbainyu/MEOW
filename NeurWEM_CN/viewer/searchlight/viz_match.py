import numpy as np, nibabel as nib
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Helvetica','Arial','DejaVu Sans'],'font.weight':'normal'})
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
from nilearn import plotting, datasets
G='/scratch/zb2164/CBI/NeurWEM/analysis/group'
im=nib.load(f'{G}/group_t_fwe05.nii.gz'); t=im.get_fdata(); t[~np.isfinite(t)]=0; aff=im.affine
THR=1e-6; VMAX=6.0; LF=16; CF=18
BG=datasets.load_mni152_template(resolution=1)
img=nib.Nifti1Image(t.astype('float32'),aff)
ch=plt.get_cmap('cold_hot'); pos_cmap=ListedColormap(ch(np.linspace(0.55,1.0,128))); neg_cmap=ListedColormap(ch(np.linspace(0.45,0.0,128)))
fig=plt.figure(figsize=(10.5,5.8)); fig.patch.set_facecolor('white')
specs=[('z',-8,[0.02,0.21,0.47,0.72]),('x',-14,[0.49,0.21,0.35,0.72])]
disps={}
for md,c,pos in specs:
    ax=fig.add_axes(pos)
    d=plotting.plot_stat_map(img,display_mode=md,cut_coords=[c],threshold=THR,cmap='cold_hot',vmax=VMAX,bg_img=BG,dim=-0.2,
        colorbar=False,axes=ax,black_bg=False,annotate=False,draw_cross=False)
    disps[(md,c)]=(d,ax)
    lbl=('Z = %d' if md=='z' else 'X = %d')%c
    ax.text(0.5,1.0,lbl,transform=ax.transAxes,ha='center',va='bottom',fontsize=CF)
    if md=='z':
        ax.text(0.11,1.0,'L',transform=ax.transAxes,ha='center',va='bottom',fontsize=CF)
        ax.text(0.89,1.0,'R',transform=ax.transAxes,ha='center',va='bottom',fontsize=CF)

def _ann(axx,items):
    for name,xy,xytext in items:
        axx.annotate(name,xy=xy,xytext=xytext,fontsize=LF,ha="center",va="center",
            annotation_clip=False,arrowprops=dict(arrowstyle="-",color="k",lw=0.9,shrinkA=0,shrinkB=3))
axZ1=list(disps[("z",-8)][0].axes.values())[0].ax
axX =list(disps[("x",-14)][0].axes.values())[0].ax
_ann(axZ1,[("early visual",(-8,-92),(-8,-156)),("fusiform",(-36,-68),(-62,-132)),("lateral occipital",(48,-80),(54,-132))])
cpos=fig.add_axes([0.880,0.52,0.013,0.26]); cb1=fig.colorbar(ScalarMappable(Normalize(0,VMAX),pos_cmap),cax=cpos); cb1.set_ticks([0,VMAX]); cb1.ax.tick_params(labelsize=14)
cneg=fig.add_axes([0.880,0.23,0.013,0.26]); cb2=fig.colorbar(ScalarMappable(Normalize(0,VMAX),neg_cmap),cax=cneg); cb2.set_ticks([0,VMAX]); cb2.ax.tick_params(labelsize=14)
fig.text(0.862,0.80,'pos',ha='right',va='center',fontsize=16,color='#b11')
fig.text(0.862,0.51,'neg',ha='right',va='center',fontsize=16,color='#15c')
fig.savefig(f'{G}/fig_searchlight_t_match.png',dpi=240,facecolor='white')
fig.savefig(f'{G}/fig_searchlight_t_match.pdf',facecolor='white')
print('wrote fig_searchlight_t_match.png+pdf',flush=True)
