# Reproduction code — *Episodic memory rescues working memory via pattern separation, pattern completion, and predictive recall*

Bai, Fougnie & Michelmann.

These scripts reproduce every statistic reported in the manuscript from the
**group-level data** (Zenodo doi: 10.5281/zenodo.22017191). They are a clean,
consolidated re-implementation of the analysis pipeline; all statistics use one
consistent convention (see below). The heavier per-trial preprocessing (raw
parsing, gaze-density and pupil preprocessing) is not repeated here — the scripts
start from the deposited group-level `.mat` files.

## Layout

```
zenodo/
  run_all.m              % master script — runs every module
  repro_behavioral.m     % Figure 2, Supplementary Figure 1
  repro_pupil.m          % Figure 3, Supplementary Figure 3
  repro_gaze.m           % Figure 4d–f
  repro_cumulative.m     % Figure 5
  repro_entropy.m        % Supplementary Information (spatial entropy)
  repro_controls.m       % Self-selection control (image similarity + lure bin)
  make_image_similarity.m% builds the objective A–B similarity lookup from the images
  lib/                   % shared statistical helpers (see "Reporting convention")
```

## How to run

1. Put this `zenodo/` folder inside `Experiment_1/` (beside `results/`, `data/`,
   `stimulus/`), **or** edit the three paths at the top of `run_all.m`.
2. In MATLAB: `run_all`.
3. Each module prints its statistics to the console. Compare against the manuscript.

Modules are independent; you can run any one on its own (each adds `lib/` to the path via `run_all`, so run `run_all` once or `addpath lib` first).

## Reporting convention (lib/)

All ANOVAs are **repeated-measures ANOVAs** (`fitrm`/`ranova`):

- **Integer degrees of freedom** for every *F*.
- **Greenhouse–Geisser** correction applied to the *p*-value of any effect with
  numerator df > 1 (the 3-level condition effects), with ε reported; 1-df
  (2-level) effects need no correction.
- Effect size = **partial η²ₚ** = SS_effect / (SS_effect + SS_error).
- Directional, a-priori predictions are tested **one-tailed** (graded condition
  ordering; correct-vs-incorrect gaze/pupil contrasts; cumulative slopes vs zero;
  WM–EM correlations). Control / no-difference analyses are two-tailed.
- Multiple comparisons within a family use **Benjamini–Hochberg FDR**.

Helpers: `rm_oneway` (3-level one-way), `rm_2x2`, `rm_3x2`, `gg_eps`
(Greenhouse–Geisser ε), `pstr`, `bh_fdr`, `cohend`.

## Data files expected (from the group-level deposit)

| module            | file(s)                                                                    |
|-------------------|----------------------------------------------------------------------------|
| behavioral        | `results/all_subjs_stats.mat`                                              |
| gaze              | `results/gaze_reinstat_res_ab.mat`, `results/gaze_reinstat_res_full.mat`   |
| cumulative        | `results/cumu_reinstat_aa_cor.mat`, `results/cumu_reinstat_ba.mat`         |
| entropy           | `results/spatial_entropy_results.mat`                                      |
| controls          | `results/image_similarity_AB.csv` (auto-built), `data/eye_movement_data/group_eye_movement_combined.mat` |
| pupil             | `results/pupil_window_means.mat` (see below)                              |

## Pupil note

The 1.0–1.5 s analysis window is chosen by a cluster-based permutation test that
needs the full per-trial pupil time series; that step lives in the source script
`S_pupil_cluster.m`. To export the per-subject window means this package reads —
**without modifying any existing script** — run, in order:

```matlab
S_pupil_cluster            % cluster-based window selection; leaves the window means in the workspace
export_pupil_window_means  % writes results/pupil_window_means.mat
```

`repro_pupil.m` then reproduces the 3×2, same-detection ANOVA, and post-hocs.

## Dependencies

- MATLAB with the **Statistics and Machine Learning Toolbox** (`fitrm`, `ranova`).
- **Image Processing Toolbox** for `make_image_similarity.m` (`rgb2gray`, `imresize`).
- Optional: **bayesFactor** toolbox (Krekelberg) for the reported Bayes factors;
  without it the frequentist statistics still print and BF lines are skipped.

## Figures

Statistics are reproduced here. The published figures are produced by the
figure scripts in `scripts/` (e.g. `S_fig1_gaze_tracks_similarity.m`,
`S_cumu_reinstat.m`); this package intentionally separates statistics from plotting.
