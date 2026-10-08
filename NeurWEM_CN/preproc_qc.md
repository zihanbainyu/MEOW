# NeurWEM_CN preprocessing QC — 6 subjects

fMRIPrep 25.2.1, MNI152NLin2009cAsym 2mm, TR 1.5 s. nback = 1back×4 + 2back×4 (8 runs); MST = 1 run. Folders now uniform `s101`–`s106` (`s102_full` renamed to `s102`).

| sub | runs | MST | mean FD (nback) | coverage I/U | status | verdict |
|-----|------|-----|-----------------|--------------|--------|---------|
| 101 | 8 + MST | yes | 0.13–0.44 | 0.85 clean | not processed | **MST unusable** (FD 0.96); nback usable, 2back run-04 marginal (FD 0.60) |
| 102 | 9 (old TASK naming) | yes | 0.08–0.16 | 0.89 clean | in analysis/data | clean, reference subject |
| 103 | 8 + MST | yes | 0.13–0.25 | 0.87 clean | not processed | usable, ready |
| 104 | 8 | **no** | 0.17–0.30 | 0.84 clean | not processed | usable, nback-only |
| 105 | 8 + MST | yes | 0.10–0.28 | 0.14 (drift) | in analysis/data | usable **with 1back run-04 excluded** (done → 175k vox) |
| 106 | 8 | **no** | 0.11–0.22 | 0.53* clean | in analysis/data | usable; check 1back run-02 registration |

## Per-subject notes

**101** — all 9 runs present and coverage is clean (I/U 0.85). The problem is motion: MST mean FD = 0.96 (56 vols > 0.5 mm) → MST unusable. 2back run-04 also elevated (FD 0.60, 38 vols > 0.5). The 8 nback runs otherwise sit at 0.13–0.44. The model-RSA searchlight uses nback only, so 101 can still enter that group test if 2back run-04 is dropped or censored; its MST cannot be used.

**102** — gold standard. Old scanner naming (`task-TASK` run-01…09, 9 runs = 8 nback + MST). Mean FD 0.08–0.16, coverage I/U 0.89. Already extracted into analysis/data (betas_wb 790+ trials).

**103** — full acquisition (8 nback + MST). Motion good (0.13–0.25; both run-04s slightly higher, still fine). Coverage clean (I/U 0.87). Not yet through LS-S / analysis-space.

**104** — **no MST scan** (8 nback runs only). Motion good overall; end-of-session drift on both run-04s (1back run-04: 25 vols > 0.5; 2back run-04: 23) but means stay ≤0.30. Coverage clean (I/U 0.84). Nback-only pipeline, same as 106.

**105** — full acquisition. Motion good (0.10–0.28). Per-run masks are all normal size (~265–281k) but spatially offset across runs → whole-set intersection collapses to 63k (I/U 0.14). 1back run-04 is the worst offender; it is already excluded, rebuilding betas_wb to 175k voxels (occipital 100%, hippocampus 97%, DLPFC ~60% from a separate anterior clip). Usable as-is in analysis/data.

**106** — **no MST scan** (8 nback runs only). Motion good (0.11–0.22). Intersection across runs is healthy (263k ≈ 97% of smallest normal run); the low I/U (0.53) is driven by one over-inclusive mask — 1back run-02 at 491k voxels (~2× the others, likely a skull/neck-including brain extraction, not FD: its FD is 0.12). Intersection-based betas_wb is unaffected, but eyeball run-02's registration/mask in `sub-106.html` before fully trusting that run. In analysis/data.

## Group readiness

- **In analysis/data already:** 102, 105, 106.
- **Need processing** (events → LS-S → betas_wb → analysis space): 101, 103, 104.
- **MST present:** 101 (motion-corrupt), 102, 103, 105. **MST absent:** 104, 106.
- For the model-RSA **searchlight** (nback only): 102, 103, 104, 105, 106 are solid; 101 includable with 2back run-04 dropped → up to 6 for the group test.
