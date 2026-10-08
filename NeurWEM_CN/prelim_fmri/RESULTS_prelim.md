# MST-back fMRI — preliminary results (N = 6)

Six participants (sub-101/102/103/104/105/106) contribute to the model-RSA searchlight; five
(sub-101/103/104/105/106) carry the trial-wise responses needed for the accuracy-dependent
analyses (sub-102 was acquired without response logging). Single-trial responses were estimated
with LS-S over the eight n-back runs in MNI152NLin2009cAsym (2 mm). Analyses are reported in
standard space with a-priori anatomical ROIs (Harvard-Oxford); hippocampus is treated as whole
and along the anterior/posterior long axis, not at the subfield scale, which the 2 mm functional
resolution does not support.

## Behavior
Similar-pair discrimination (B2 probe) is more accurate when the pair was previously compared than
when it is novel (compared 57%, novel 50%, n = 5), reproducing the direction of the behavioral
EM-benefit the pilot was designed to localize.

## Figure 1 — The brain tracks image similarity (localizer)
A whole-brain model-RSA searchlight (neural 1-back × 2-back pattern similarity against a graded
identity model: same image = 1, similar pair-mate = 0.5, otherwise = 0) was tested at the group
level with a sign-flip permutation (all 2^6 flips, cluster-forming p < .01, cluster-extent FWE).
A single cluster survives correction: a large bilateral ventral-occipital visual region
(7,494 voxels, peak t = 18.5 at MNI [18, −76, −22], p_FWE = .031; **fig1a**). The effect is
concentrated in early and ventral visual cortex and tapers through posterior parietal and temporal
cortex. In an a-priori early-visual ROI the graded structure is present but shallow — same = 0.15,
similar = 0.14, null = 0.11 (**fig1b**): item-level similarity (same and similar > null) is tracked,
while exact repeats and similar exemplars are not separated, consistent with shared low-level
features across the perceptually similar objects. No hippocampal or parietal cluster survives, as
expected for weaker effects at this sample size.

## Figure 2 — Univariate discrimination network (primary univariate hypothesis)
Estimated with a proper condition-level first-level GLM (the four B2 cells modeled with a canonical
HRF in percent-signal-change units, same confounds as the LS-S), the discrimination effect
(correct − incorrect) is larger in the EM-available (compared) than the no-EM (novel) condition in
all four a-priori ROIs — hippocampus (+0.08 vs +0.03), angular (+0.14 vs +0.02), DLPFC (+0.11 vs
−0.02), IPS (+0.06 vs −0.02) % signal (**fig2b**). This is the predicted "accuracy modulation is
stronger when episodic memory is available" direction, consistent across regions, though not
significant at this N: the whole-brain sign-flip interaction yields no FWE-surviving cluster, with
the largest positive focus in left inferior parietal cortex (angular/supramarginal, MNI
[−50, −34, 34], p_FWE = .21; **fig2a**). Method note: this GLM supersedes an initial attempt that
averaged LS-S single-trial betas, which was unusable because raw single-trial amplitudes carry
arbitrary per-participant scaling and heavy outliers (see fig_b2_trialbetas_HC: within-subject
correct/incorrect distributions overlap on very different scales). The GLM's percent-signal-change
units are comparable across participants by construction.

## Figure 3 — Multivariate mechanisms of EM support (hypotheses 1–3)
Pattern similarity between presentations of a compared pair, split by subsequent B2 discrimination
success (**fig3**):
- **Pattern separation (A1–B1, posterior hippocampus).** Lower one-back similarity on correct
  trials, in the predicted direction — correct r = 0.45 < incorrect r = 0.52, in all five
  participants (d = −1.06; signed-rank p = .13). The temporal-proximity of consecutive A1/B1 trials
  inflates the baseline similarity, so the effect should be read as the within-pair correct-vs-
  incorrect contrast rather than the absolute level.
- **Pattern completion (B2–B1, anterior hippocampus).** No difference (d = −0.11).
- **Predictive recall (A2–B1, angular gyrus).** No reliable difference (d = −0.59, p = .13, and in
  the non-predicted direction).

## Reading of the prelim
The pilot reproduces the behavioral EM-benefit and cleanly localizes a visual similarity-tracking
network (the method works and is group-significant at N = 6). Of the three mechanistic predictions,
only pattern separation in posterior hippocampus shows the predicted direction, consistently across
participants though not significant at this N. The univariate shared-subnetwork prediction is not
yet supported. These are small-sample, standard-resolution results: the sign-flip searchlight has a
permutation floor of p = 1/64 = .016 (1/32 for the five-participant interaction), and subfield-level
and native-space tests await the full protocol (1.5 mm functional + coronal T2 + ASHS).
