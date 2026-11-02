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
