# How to analyze the MST-back fMRI data: an evidence-based method review

Scope: this maps each of your hypotheses to the way the field actually tests that claim, cites the
method literature, and — per your instruction — judges whether each hypothesis is appropriately
testable with *this* data (2 mm isotropic, MNI space, event-related, small but growing N, no
subfield-resolution T2 in the pilot). Where a hypothesis is confounded or non-standard as stated, I
say so and give the validated alternative. Nothing here assumes your hypotheses are correct.

## Bottom line up front
1. The conceptual framework (separation / completion / predictive recall supporting WM discrimination)
   is coherent and largely well-precedented — **but the three mechanistic tests are not equally
   sound.** Predictive recall and completion are strong and directly testable; pattern separation as
   you've framed it is the weakest.
2. The single most important method issue is **temporal autocorrelation**, and it undercuts the
   prelim numbers. Your key comparisons (A1–B1 consecutive; A2–B1, B2–B1 within session) are exactly
   the comparisons BOLD autocorrelation biases. The fix is the one your proposal already names —
   **cross-validated Mahalanobis (crossnobis) distance with cross-run folds** — which the prelim did
   not use (it used within-session Pearson correlation). Treat the prelim similarity values as biased.
3. "Pattern separation in DG" is **not resolvable at 2 mm** and is not what a raw A1–B1 correlation
   measures. At this resolution the honest granularity is whole-hippocampus and an anterior/posterior
   long-axis split; true subfields wait for the 1.5 mm + coronal-T2 protocol.

---

## Part 1 — Foundational choices (apply to every analysis)

### 1.1 Single-trial estimation
LS-S (iterative one-GLM-per-trial; Mumford et al. 2012; Turner et al. 2012) is standard and defensible.
But **GLMsingle** (Prince et al. 2022, *eLife*) is the current best practice for event-related
single-trial estimates: it fits a voxelwise HRF from a library, derives data-driven noise regressors,
and — critically for you — applies ridge regularization that **explicitly decorrelates temporally
adjacent trials**. Because every one of your key contrasts is within-session, the adjacent-trial
decorrelation matters more here than in a typical decoding study. Recommendation: re-estimate with
GLMsingle, or at minimum report that LS-S betas carry residual adjacent-trial correlation.

### 1.2 Distance / similarity metric
Pearson correlation (what the prelim used) is biased by measurement noise and by run membership.
The field standard for RSA on fMRI is **crossnobis = cross-validated Mahalanobis distance**
(Walther et al. 2016, *NeuroImage*; Diedrichsen & Kriegeskorte 2017, *PLoS Comp Biol*; Nili et al.
2014, *PLoS Comp Biol*). It is (a) multivariate-noise-normalized, (b) **unbiased by noise level**
because of cross-validation, and (c) when the cross-validation folds are *runs*, it removes the
within-run temporal autocorrelation that inflates similarity. Your proposal already specifies
cross-validated Mahalanobis — use it, and structure the folds across the four runs. For model-RDM
comparison use whitened RDM (cosine) similarity (Diedrichsen & Kriegeskorte 2017).

### 1.3 The temporal-autocorrelation confound — the main threat to your design
Empirically, single-trial BOLD patterns from the same run can correlate up to r ≈ 0.5 for consecutive
trials, decaying to ~0 over seconds (reviewed in the GLMsingle and MVPA-design literature; Mumford
et al. 2014). Consequences for your hypotheses:
- **A1–B1 (separation)**: A1 and B1 are *consecutive* in the compared condition → their similarity is
  inflated by autocorrelation, not representation. And the novel condition has no A1/B1, so there is
  no matched baseline.
- **A2–B1, B2–B1 (predictive recall, completion)**: within-session, variable temporal lag across pairs.
Mitigations, in order of strength: (i) crossnobis with cross-run folds (removes same-run bias by
construction); (ii) GLMsingle betas; (iii) **match or regress out temporal lag** between the two
presentations as a nuisance covariate; (iv) never contrast a near-in-time pair against a far-apart
pair without lag control. Any similarity result that is not lag-controlled is not interpretable.

### 1.4 Inference at small N
Do not use parametric cluster statistics at tiny N — they are unstable (your N=3 group t blew up for
exactly this reason). Use **sign-flip permutation with TFCE** via PALM/randomise (Nichols & Holmes
2002, *HBM*; Winkler et al. 2014, *NeuroImage*; Smith & Nichols 2009, *NeuroImage*). The permutation
floor is 1/2^N, so report it honestly (e.g., N=6 → min p ≈ .016). Favor effect sizes + the permutation
p over thresholded maps alone.

### 1.5 Circularity / non-independence
Kriegeskorte et al. (2009, *Nat Neurosci*) and Vul et al. (2009, *Perspect Psychol Sci*): never select
voxels/ROIs with the same data or contrast you then test. Your "localizer → ROI → mechanism test"
plan is circular if it reuses the same data. Legitimate options: define ROIs from **independent
anatomy** (Harvard-Oxford / FreeSurfer), from an **independent localizer run**, or use leave-one-run-out
cross-validation; or skip ROI selection and report whole-brain corrected maps.

---

## Part 2 — Hypothesis by hypothesis

### Localizer: regions that track image similarity (same > similar > null)
Sound and standard — this is model-RDM RSA against a graded identity model. Keep it, but its main
legitimate role is (a) a positive control that your patterns carry stimulus information, and (b) a
*data-independent* ROI definer only if estimated on held-out data. **Appropriate.**

### Primary univariate: shared subnetwork = [correct−incorrect | EM] − [correct−incorrect | novel]
A subsequent-memory × condition interaction. Precedented in principle, but three problems:
- **Difficulty/effort confound**: incorrect trials are typically harder, slower, lower-attention, so an
  activation difference can reflect task difficulty rather than EM. Control RT (trial-wise covariate),
  and check the effect is not explained by RT.
- **Power**: ~12–28 trials per cell per subject; the interaction is low-powered and noisy at single-subject level.
- **Directness**: BOLD amplitude is an indirect index of "EM availability." The multivariate reinstatement
  measures test the actual mechanism.
- **Estimation**: estimate with a proper condition-level GLM (canonical HRF, the four B2 cells modeled
  explicitly), not by averaging single-trial betas — raw LS-S amplitudes carry arbitrary per-subject
  scaling and heavy outliers (this exact failure showed up in the pilot).
**Verdict: keep as a secondary / confirmatory analysis, not the headline.**

### H1 — Pattern separation: A1–B1 lower similarity predicts B2 discrimination
Weakest of the three, for three reasons:
1. **Non-standard operationalization.** The field's validated signature of hippocampal pattern
   separation is *activation-based*: DG/CA3 treat a similar lure like a novel item (reduced
   repetition suppression / "lure" response), measured against targets and foils (Bakker et al. 2008,
   *Science*; Kirwan & Stark 2007; Yassa & Stark 2011, *Trends Neurosci*). A raw A1–B1 pattern
   correlation is an indirect proxy for "the pair was encoded distinctly," not the canonical measure.
2. **Confounds.** A1–B1 are consecutive (autocorrelation ↑ similarity), the novel condition lacks any
   A1/B1 baseline, and lower correlation can simply mean noisier (less reliable) patterns — a
   reliability artifact that crossnobis, not Pearson, controls.
3. **Resolution.** "DG" is unresolvable at 2 mm (subfields need ≤1 mm / high-res coronal T2; Carr,
   Rissman & Wagner 2010, *Neuron*; partial-volume mixing dominates ≥2 mm). At 2 mm the honest unit is
   whole-HC or posterior HC (posterior = detail/separation on the long axis; Poppenk et al. 2013,
   *Trends Cogn Sci*; Strange et al. 2014, *Nat Rev Neurosci*).
**Verdict: as stated, not cleanly testable.** Either (a) test separation the validated way — the
mnemonic-discrimination activation contrast (lure vs. target/foil) from the final in-scanner MST task,
in HC; or (b) if you keep the pattern-similarity version, use crossnobis with cross-run folds, match
temporal lag, restrict to anterior/posterior HC, and call it "neural distinctiveness of the pair,"
not DG pattern separation. Defer true subfield separation to the full protocol.

### H2 — Pattern completion: B2–B1 higher similarity predicts success
This is **encoding–retrieval similarity (ERS) / reinstatement**, which is well-precedented and sound
(Ritchey et al. 2013, *Cereb Cortex* — hippocampus mediates cortical ERS–memory link; Gordon et al.
2013; Wing et al. 2015 — item-level ERS in occipitotemporal cortex predicts memory). Requirements:
- **Item-specific baseline**: contrast same-item ERS (B2–B1) against different-item ERS, so the effect
  is reinstatement of *this* item, not generic similarity. Crossnobis handles the noise/autocorrelation.
- CA3 is unresolvable at 2 mm, but completion is also expected in cortex — temporal/parietal retrieval
  network (Rugg & Vilberg 2013) and angular gyrus — which *is* resolvable at 2 mm.
**Verdict: appropriate and well-grounded.** Run as item-level ERS with crossnobis, success-split, in
a-priori HC + angular/lateral-parietal ROIs and a whole-brain searchlight.

### H3 — Predictive recall: A2–B1 higher similarity predicts success
**The strongest and most novel of the three.** Viewing A2 reinstating the *pairmate's* earlier pattern
(B1) *before* B2 appears is exactly anticipatory/associative reinstatement, which has direct fMRI
precedent: hippocampal pattern completion driving anticipatory reinstatement in sensory cortex (Hindy,
Ng & Turk-Browne 2016, *Nat Neurosci*; Kok & Turk-Browne 2018; Kok et al. 2012, *Neuron*). Notes:
- **Measure**: similarity(A2 pattern, B1 pattern), item-baselined against similarity(A2, other B1s),
  crossnobis. Expected in lateral parietal / precuneus (and HC) — all resolvable at 2 mm.
- **Key confound**: A2 and B1 are *different* items but perceptually similar, so shared low-level
  features could inflate similarity independent of recall. Two defenses: the similar-but-unstudied
  baseline, and — decisively — the **success-dependence** test (correct > incorrect), because
  low-level perceptual similarity cannot explain why the effect tracks subsequent discrimination.
**Verdict: appropriate, best-grounded, and the one I would prioritize.**

---

## Part 3 — The analysis I would actually run

1. **Betas**: re-estimate single-trial responses with GLMsingle (ridge-regularized, adjacent-trial
   decorrelated). Keep LS-S as a robustness check.
2. **Localizer RSA** (identity model, crossnobis, searchlight, PALM + TFCE): positive control + optional
   independent ROIs from held-out runs.
3. **Primary multivariate** — all crossnobis, cross-run folds, item-baselined, split by B2 success:
   - **Predictive recall A2–B1** (priority): lateral parietal / precuneus + HC; success-dependence is
     the critical test.
   - **Completion B2–B1** (ERS): HC + angular / lateral parietal; item-specific.
   - **Separation A1–B1**: reframed as pair distinctiveness in anterior/posterior HC, lag-controlled;
     and/or the activation-based lure contrast from the final MST task.
4. **Secondary univariate interaction**: condition-level GLM (not averaged single-trial betas),
   RT-controlled, permutation inference. Confirmatory only.
5. **ROIs**: a-priori anatomical — HC (whole + anterior/posterior long-axis split), angular, precuneus,
   IPS, DLPFC — to avoid circularity. Subfields only with the full-protocol high-res T2 + ASHS.
6. **Inference**: PALM sign-flip + TFCE; report effect sizes and the permutation floor; state small-N
   limits plainly.

## Part 4 — Honest verdict on the hypotheses
- **Coherent and mostly well-precedented**, but not uniform in quality.
- **Predictive recall (H3)** and **completion (H2)** are the appropriate, directly-testable, novel
  tests with this data — prioritize them, as ERS/reinstatement with crossnobis.
- **Pattern separation (H1)** as a raw A1–B1 similarity is confounded (autocorrelation, no matched
  baseline, reliability), non-standard (the validated measure is activation-based lure discrimination),
  and mis-scaled (not "DG" at 2 mm). Reframe it or test it the validated way.
- **Univariate interaction** is legitimate but secondary and confounded by difficulty.
- **Biggest single fix**: crossnobis + cross-run cross-validation + lag control. Your proposal already
  commits to cross-validated Mahalanobis; the prelim's within-session Pearson numbers are biased and
  should not be over-interpreted.

## Why this matters for framing (WM–EM)
The broader claim — that EM contributes to WM — is itself contested and worth citing directly:
Beukers et al. (2021, *Trends Cogn Sci*, "Is activity-silent working memory simply episodic memory?")
and Lewis-Peacock & Postle (2008, *J Neurosci*) argue that long-term/episodic traces support
ongoing WM. Anchoring your design in that debate strengthens the significance and pre-empts the
reviewer question of whether you are measuring EM or just WM.

---

## Selected references
- Bakker A, Kirwan CB, Miller M, Stark CEL (2008). Pattern separation in the human hippocampal CA3 and dentate gyrus. *Science*.
- Beukers AO, Buschman TJ, Cohen JD, Norman KA (2021). Is activity-silent working memory simply episodic memory? *Trends Cogn Sci*.
- Carr VA, Rissman J, Wagner AD (2010). Imaging the human medial temporal lobe with high-resolution fMRI. *Neuron*.
- Diedrichsen J, Kriegeskorte N (2017). Representational models: a common framework... *PLoS Comput Biol*.
- Gordon AM, Rissman J, Kiani R, Wagner AD (2013). Cortical reinstatement and the confidence/accuracy of memory. *Cereb Cortex*.
- Hindy NC, Ng FY, Turk-Browne NB (2016). Linking pattern completion in the hippocampus to predictive coding in visual cortex. *Nat Neurosci*.
- Kok P, Turk-Browne NB (2018). Associative prediction of visual shape in the hippocampus. *J Neurosci*.
- Kriegeskorte N, Mur M, Bandettini P (2008). Representational similarity analysis. *Front Syst Neurosci*.
- Kriegeskorte N, Simmons WK, Bellgowan PSF, Baker CI (2009). Circular analysis in systems neuroscience: the dangers of double dipping. *Nat Neurosci*.
- Lewis-Peacock JA, Postle BR (2008). Temporary activation of long-term memory supports working memory. *J Neurosci*.
- Mumford JA, Turner BO, Ashby FG, Poldrack RA (2012). Deconvolving BOLD activation in event-related designs for multivoxel pattern analysis. *NeuroImage*.
- Nichols TE, Holmes AP (2002). Nonparametric permutation tests for functional neuroimaging. *Hum Brain Mapp*.
- Nili H, Wingfield C, Walther A, Su L, Marslen-Wilson W, Kriegeskorte N (2014). A toolbox for representational similarity analysis. *PLoS Comput Biol*.
- Poppenk J, Evensmoen HR, Moscovitch M, Nadel L (2013). Long-axis specialization of the human hippocampus. *Trends Cogn Sci*.
- Prince JS, Charest I, Kurzawski JW, Pyles JA, Tarr MJ, Kay KN (2022). Improving single-trial fMRI response estimates using GLMsingle. *eLife*.
- Ritchey M, Wing EA, LaBar KS, Cabeza R (2013). Neural similarity between encoding and retrieval is related to memory via hippocampal interactions. *Cereb Cortex*.
- Rugg MD, Vilberg KL (2013). Brain networks underlying episodic memory retrieval. *Curr Opin Neurobiol*.
- Smith SM, Nichols TE (2009). Threshold-free cluster enhancement. *NeuroImage*.
- Strange BA, Witter MP, Lein ES, Moser EI (2014). Functional organization of the hippocampal longitudinal axis. *Nat Rev Neurosci*.
- Vul E, Harris C, Winkielman P, Pashler H (2009). Puzzlingly high correlations in fMRI studies of emotion, personality, and social cognition. *Perspect Psychol Sci*.
- Walther A, Nili H, Ejaz N, Alink A, Kriegeskorte N, Diedrichsen J (2016). Reliability of dissimilarity measures for multi-voxel pattern analysis. *NeuroImage*.
- Wing EA, Ritchey M, Cabeza R (2015). Reinstatement of individual past events revealed by the similarity of distributed activation patterns during encoding and retrieval. *J Cogn Neurosci*.
- Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE (2014). Permutation inference for the general linear model. *NeuroImage*.
- Yassa MA, Stark CEL (2011). Pattern separation in the hippocampus. *Trends Neurosci*.
