# LS-S single-trial estimation: exact pipeline, mapped to the code

This is the full chain from raw data to single-trial betas, step by step, each step tied to the exact
lines in `lss_new.py`. The method is Least-Squares-Separate (LS-S; Mumford et al. 2012, *NeuroImage*).

---

## Stage 0 — What fMRIPrep produced (the inputs; before any of this code runs)
The betas are estimated on fMRIPrep 25.2.1 outputs, not on scanner raw. For each run fMRIPrep did:
motion correction, susceptibility-distortion correction, slice-timing correction (reference = middle
of the TR), coregistration, and normalization to `MNI152NLin2009cAsym` at 2 mm. It wrote:
- `*_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz` — the preprocessed 4-D series.
- `*_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz` — the brain mask.
- `*_desc-confounds_timeseries.tsv` (+ `.json`) — nuisance regressors.
No spatial smoothing is applied at any point (required for pattern analysis).

The event timing comes from `events.tsv` (built from the behavioral logs): `onset` in seconds from the
first volume, `duration = 1.5 s`.

---

## Stage 1 — Load timing, confounds, and the BOLD matrix

**Confound selection** — `pick_confounds()`:
- `motion = [trans_x, trans_y, trans_z, rot_x, rot_y, rot_z] + their _derivative1` → **12 motion
  regressors** (6 rigid-body + 6 temporal derivatives).
- `acc = a_comp_cor_* with Mask=='combined', sorted by SingularValue, top 6` → **6 aCompCor** components
  (combined WM+CSF anatomical noise).
- `cosine = columns starting 'cosine'` → fMRIPrep's **discrete-cosine high-pass regressors**. Including
  these *is* the high-pass filter, which is why `drift_model=None` later (we don't add a second filter).
- `outl = 'motion_outlier*' / 'non_steady*'` → one-hot **spike/censor regressors** for high-motion and
  dummy volumes.
- Returns the matrix `conf` (time × nuisance) and the column names `conf_cols`.

**Load BOLD into a 2-D matrix:**
```
nvols = nib.load(bold).shape[3]; n = len(ev)                 # #volumes, #trials
frame_times = np.arange(nvols)*TR + TR/2.0                   # TR=1.5; sampling at mid-TR (0.5 ref)
masker = NiftiMasker(mask_img=mask, standardize=False, t_r=TR)
Y = masker.fit_transform(bold)                               # Y is (time × voxels), raw BOLD units
```
- `frame_times` places each volume at the **middle of its TR** (`+TR/2`), matching fMRIPrep's
  slice-timing reference.
- `standardize=False` → betas stay in native BOLD units (no z-scoring). No `smoothing_fwhm` → no smoothing.
- `Y` = every in-mask voxel's time course, one column per voxel.

---

## Stage 2 — The LS-S loop: one GLM per trial

`betas = np.zeros((n, nvoxels))` then `for i in range(n):` — a **separate GLM is fit for each trial i**.

### 2a. Build the trial's two-condition event table
```
ev_i = DataFrame(onset=all onsets, duration=all durations, trial_type=['other']*n)
ev_i.loc[i,'trial_type'] = 'trial'
```
Every trial is labeled `'other'` **except trial i**, which is `'trial'`. This is the defining LS-S move:
the design has exactly **two task conditions** — trial i by itself, and all remaining trials collapsed
into a single regressor.

### 2b. HRF convolution → the design matrix (`make_first_level_design_matrix`)
```
dm = make_first_level_design_matrix(frame_times, ev_i, hrf_model='glover',
                                     drift_model=None, add_regs=conf.values, add_reg_names=conf_cols)
```
This is where the HRF convolution happens. For each condition, nilearn:
1. builds a boxcar = 1 from `onset` to `onset+duration` (1.5 s), 0 elsewhere;
2. convolves that boxcar with the **Glover canonical HRF** `h(t)` (no temporal/dispersion derivatives);
3. samples the convolved signal at `frame_times`.

So the two task columns are, in continuous form:
- **`trial`** = (boxcar of trial *i*) ⊛ h, sampled at frame_times.
- **`other`** = (Σ over all trials *j≠i* of boxcar_j) ⊛ h, i.e. the summed predicted response of every
  other trial, as one regressor.

The full design `dm` (its columns, in order) is:
`[other, trial, <12 motion>, <6 aCompCor>, <cosine…>, <outliers…>, constant]`
(nilearn alphabetizes the two conditions → `other` before `trial`, and adds the `constant` intercept).
`drift_model=None` because the cosine high-pass is already supplied through `add_regs`.

### 2c. Fit the GLM (ordinary least squares)
```
labels, est = run_glm(Y, dm.values, noise_model='ols')
```
Solves `Y = dm · β + ε` per voxel by OLS — one β per design column per voxel.

### 2d. Pull out this trial's amplitude
```
con = zeros(ncols); con[index of 'trial'] = 1
betas[i,:] = compute_contrast(labels, est, con, contrast_type='t').effect_size()
```
The contrast selects only the `trial` column, so `effect_size()` returns **β̂_trial** — the estimated
response amplitude of trial *i* at every voxel. (`contrast_type='t'` only affects the test statistic,
not the effect size we store.) That voxel vector is row *i* of `betas`.

The loop repeats for all `n` trials, each with its own refit GLM.

---

## Stage 3 — Outputs
```
np.save('sub-{S}_task-{T}_run-{R}_betas.npy', betas)   # (n_trials × n_voxels), float32
ev.to_csv('..._trials.tsv')                             # the trial metadata, aligned row-for-row
```
Downstream, `build_betas_wb_new.py` intersects the per-run brain masks, stacks the eight n-back runs
into `betas_wb.mat` (trials × common-mask voxels) with per-trial labels.

---

## Why this is textbook LS-S (the point to make to your colleagues)
- **LS-A** (Least-Squares-All) would put *n* separate trial regressors in **one** GLM. In a rapid design
  (median neighbor SOA here ≈ 2.6 s) those regressors are massively collinear → variance inflation
  (we measured LS-A VIF ≈ 31 median, >98% of trials above 5).
- **LS-S** refits *n* GLMs, each with **one** trial regressor plus **one** combined "all-other-trials"
  regressor. Collapsing the neighbors into a single column removes the trial-by-trial collinearity, at
  the cost of a small, well-characterized bias. That is exactly why the `'trial'` VIF is ≈ 1.3 here —
  the designed behavior of LS-S, not an artifact (Mumford et al. 2012; Abdulrahman & Henson 2016).
- Everything else is standard: Glover HRF convolution, fMRIPrep confounds (24-param-style motion +
  aCompCor + cosine high-pass + spike censors), OLS, no smoothing, native units.

### Exact parameter list (for the methods text)
TR = 1.5 s; Glover canonical HRF, no derivatives; no drift model (cosine high-pass via fMRIPrep);
frame-time reference = mid-TR; mask = fMRIPrep MNI brain mask; no spatial smoothing; `standardize=False`;
nuisance = 6 motion + 6 motion derivatives + 6 aCompCor (combined) + cosine + motion-outlier/non-steady
spikes; OLS fit; single-trial amplitude = `trial` contrast effect size; one GLM per trial.
