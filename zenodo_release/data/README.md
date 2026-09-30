# Group data

*Episodic memory rescues working memory via pattern separation, pattern completion, and predictive recall.* Bai, Fougnie & Michelmann. Zenodo doi: 10.5281/zenodo.22017191.

Group-level statistical data underlying the manuscript, one file per modality.
No individual-level data (per-fixation eye traces, per-trial pupil series) is
included. Reproduce all reported statistics with the code archive (doi:
10.5281/zenodo.22017336): place its `scripts/` folder next to this `data/` folder
and run `run_all.m`.

## files

`behavior.mat` (`behavior`)
: per-subject summaries. one-back accuracy and rt; two-back discrimination index,
  same-item d', reaction times; recognition d'. n = 31.

`gaze.mat` (`gaze`)
: `reinst_ab`  a1-b1 gaze similarity per trial with subsequent b2 accuracy.
  `reinst_full`  b2-b1 and a2-b1 gaze similarity per trial with accuracy.
  `cumu_aa`, `cumu_ba`  cumulative gaze similarity by fixation number, per subject.
  `entropy`  per-subject spatial gaze entropy (a/b x compared/isolated).
  `img_sim`  objective a-b image similarity per pair (base_id, bin, pix_corr).
  `control`  per two-back a-b trial: subj_id, pair_key, condition, b2 correct.
  `bin_lookup`  per one-back b trial: subj_id, trial_id, lure bin (l1/l2).

`pupil.mat` (`pupil`)
: per-subject baseline-corrected pupil time-courses. `pup` a-b lure and a-a target
  by condition; `corr` / `incorr` a-b lure by accuracy; `t` time vector; `subj`
  subject ids; `win` analysis window [1.0 1.5] s; `clusters` precomputed
  cluster-based permutation results (omnibus and pairwise), rows [t0 t1 mass p]. n = 29.
