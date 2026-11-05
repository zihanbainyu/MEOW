# NeurWEM neural pipeline (N=9, sub-101..109)

Input: LS-S single-trial betas at analysis/data/sub-NN/betas_wb.mat (MNI152NLin2009cAsym 2mm).
Behavior fix: patch_102.py (sub-102 resp_key realignment).

## 1. Searchlight model-RSA  -> visual cortex (reportable)
sl_rho.py SUBJ            (array: sl_rho.sbatch)   -> group/sl_rho_sub-NN.nii.gz
sl_tfce.py               (sl_tfce.sbatch, B=5000 sign-flip TFCE FWE) -> group/group_t_fwe05.nii.gz, group_tfce*
viz_match.py             -> group/fig_searchlight_t_match.pdf
peaks.py / cluster_table.py / region_count.py -> peak + region tables
Result: FWE<.05 in early visual + lateral occipital + fusiform only.

## 2. ROI condition x accuracy 2x2  (B2 probe, correct vs incorrect x compared vs novel)
roi_2x2.py               (roi_2x2.sbatch)          -> group/prelim/roi_2x2_persubj.csv
local: roi_2x2_anova.m   (rmANOVA interaction + paired-t post-hoc, two-sided)
local: python plotting   -> fig_{dlpfc,angular,hippocampus}_2x2.pdf
Result (trend, N=9): DLPFC F(1,8)=4.37 p=.070; lateral-temporal p=.072; angular p=.101; hippocampus n.s.

## 3. Whole-brain condition x accuracy GLM
glm_uni.py SUBJ          (array: glm_uni_arr.sbatch) -> group/glm_{interaction,compeff,noveleff}_sub-NN.nii.gz
glm_group_tfce.py        (glm_group_tfce.sbatch)     -> group/group_glm_interaction_{t,fwe05}.nii.gz
viz_interaction.py / viz_inter_match.py / interaction_loc.py -> maps + localization
Result: 0 voxels survive FWE; uncorrected map is noise -> not reportable.

archive/ holds superseded scripts (old searchlight group, surface renders, dlpfc/hc_glm/b2 iterations).

---

# Local layout & reproduction

Directory (type -> phase; analysis type is in the filename: fig_2x2_*, fig_rsa_*, fig_behav_*, fig_mem_*):

    scripts/mstback/   scripts/mst/     analysis + figure scripts
    data/mstback/      data/mst/        inputs (.csv/.mat/.npz); + mstback/fc, mstback/brainnet
    figures/mstback/   figures/mst/     outputs (.pdf/.png); + backup/
    maps/              maps/hires/      whole-brain group t-maps + bg_MNI2009cAsym_T1 + tpl (hires = source maps + geom)
    icons/                              ROI icons
    scripts/archive/                    legacy scripts (external /Users/bai/Downloads data; not repo-reproducible)

Reproduce (run from prelim_fmri/):

    bash sync_from_cluster.sh                 # pull group summaries + icons into data/ and icons/
    matlab -batch "run('scripts/mstback/make_prelim_panels.m')"   # 2-back 2x2 + RSA + MST-mem + behaviour panels
    matlab -batch "run('scripts/mstback/make_contrasts.m')"       # whole-brain task-contrast montages (reads maps/hires)
    matlab -batch "run('scripts/mstback/coupling_brain.m')"       # 3D HC->cortex coupling (reads data/mstback/coupling_meshes.mat)
    matlab -batch "run('scripts/mstback/{connectogram,heatmap,matrix,compose_brain}.m')"  # FC figures
    matlab -batch "run('scripts/mst/make_mst_{idx,m3,behav,sep,combined}.m')"             # MST figures

All script paths are relative to the script location, so any of the above runs from prelim_fmri/ without edits.
Whole-brain interaction/MST t-maps for eyeballing: maps/grp_interaction_discrim_t, grp_interaction_samedetect_t, grp_mst_emlure_acc_t (overlay on maps/bg_MNI2009cAsym_T1).
