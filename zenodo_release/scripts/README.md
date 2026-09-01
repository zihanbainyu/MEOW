# Reproduction code

*Episodic memory rescues working memory via pattern separation, pattern completion, and predictive recall.* Bai, Fougnie & Michelmann.

These scripts reproduce every statistic in the manuscript from the consolidated
group data. Zenodo doi: 10.5281/zenodo.22017336 (code); 10.5281/zenodo.22017191 (data).

## layout

Put the two archives side by side:

```
scripts/   this archive; run run_all.m from here
data/      the data archive: behavior.mat, gaze.mat, pupil.mat
```

If your folders differ, edit `MEOW_DATA` at the top of `run_all.m`.

## how to run

In MATLAB, from `scripts/`: `run_all`. It prints all statistics to the console;
compare against the manuscript. Modules are independent and use only the three
group files.

| module          | reproduces                                         |
|-----------------|----------------------------------------------------|
| `repro_behavior`| figure 2, supplementary figure 1                   |
| `repro_pupil`   | figure 3, supplementary figure 3                   |
| `repro_gaze`    | figure 4, figure 5, supplementary entropy, controls|

## reporting convention (`lib/`)

- repeated-measures anova (`fitrm`/`ranova`), integer degrees of freedom
- greenhouse-geisser on the p when numerator df > 1 (3-level condition effects), epsilon reported; 1-df effects need no correction
- effect size: partial eta^2 = ss_effect / (ss_effect + ss_error)
- a-priori directional predictions one-tailed; control / no-difference tests two-tailed
- multiple comparisons within a family: benjamini-hochberg fdr

Helpers: `rm_oneway`, `rm_2x2`, `rm_3x2`, `gg_eps`, `pstr`, `bh_fdr`, `cohend`.

## how the group data was produced

The three data files were generated once from the raw individual-level recordings
(per-fixation eye table and per-trial pupil series, which are not deposited),
including the precomputed pupil cluster-permutation results. The generation scripts
are kept with the authors and are not needed to reproduce anything here.

## dependencies

- MATLAB with the Statistics and Machine Learning Toolbox (`fitrm`, `ranova`)
- Image Processing Toolbox for `make_group_data` only (reads the stimulus images)
- the bayesFactor toolbox (Krekelberg) is bundled in `lib/bayesFactor-master` and
  added to the path automatically; it provides the reported bayes factors

## notes

- pupil cluster p-values are permutation-based; the seed is fixed (`rng(1)`) so
  they are reproducible.
- bayes factors here are computed by deterministic numerical integration
  (quadrature) and are therefore reproducible. they may differ slightly from the
  values reported in the manuscript, which used the toolbox's default monte-carlo
  integration; this reflects monte-carlo estimation error (typically about 1
  percent, at most a few percent) and does not change any conclusion.
- this archive reproduces the statistics; the published figures are drawn by the
  figure scripts in the original analysis repository.
