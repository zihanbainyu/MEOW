# Hybrid MST–N-back fMRI Task — Procedure

*Draft method description for review. Zihan Bai, Michelmann Lab, NYU.*

## Overview

A single scanning session comprises two parts. **Part 1** interleaves a 1-back and a
2-back object task across four blocks. **Part 2** is a post-task mnemonic-similarity
recognition test (MST). Eye position and pupil size are recorded throughout with an
EyeLink 1000 (1000 Hz). Each run is time-locked to the scanner: it begins on the first
trigger pulse and all events are modeled from that onset.

## Stimuli

Everyday-object photographs from the Stark mnemonic-similarity set, luminance-matched
and resized, presented centrally on a mid-gray background (RGB 124/124/124). The pool is
organized into four pair sets of 120 pairs each — *compared*, *isolated*, *novel*, and
*repeat* — with similar lures drawn from graded similarity bins (higher- and
intermediate-similarity bins for the compared pairs; a separate bin for the repeat
pairs). The fixation target is the combined bullseye-and-crosshair figure recommended for
stable fixation (Thaler et al., 2013).

## Part 1 — Interleaved 1-back and 2-back (4 blocks)

Each block runs a **1-back** encoding run followed by a **2-back** retrieval run over the
same block-specific stimulus set. In both runs the participant judges whether the current
object matches an earlier one — **same** (index finger, button 1), **similar** (middle
finger, button 2), or **new** (no response).

**1-back run — 180 trials/block.** Each object is compared to the one immediately before
it:

- **Repeat** (30/block): an object repeated on consecutive trials → *same*.
- **Compared** (30/block): an object followed immediately by a similar lure → *similar*.
- **Isolated** (30 pairs/block, shown as 60 non-adjacent singletons) → *new*.

Trials are grouped into inseparable two-trial mini-blocks so that no same/similar
comparison — and no mid-run break — ever splits a pair.

**2-back run — ~150 trials/block.** Each object is compared to the one two trials back,
with the same response mapping. The stream is generated from 90 goal pairs per block,
balanced 10/10/10 across the three probe types — **same** (A–A), **similar** (A–B lure),
and **new** (A–N) — within each of three item conditions: previously *compared*,
previously *isolated*, and *novel*. Previously-compared and previously-isolated items
reappear from the same block's 1-back; novel items appear only in the 2-back. The
sequence is optimized so that match trials are aperiodically spaced and no padding trials
are needed.

The response window in both runs is the 1.5 s image presentation.

## Part 2 — Post-task MST (old / similar / new)

A standard mnemonic-similarity recognition test on the *repeat* pool. Each of the 120
repeat pairs had its item A shown twice in the 1-back; its pairmate B was never shown.
Half the pairs are probed with the seen item **A (old)** and half with the unseen
pairmate **B (lure/similar)**; a matched set of 60 fresh foils serves as **new** — **180
trials total (60/60/60)**. The participant judges each object as **OLD** (index finger),
**SIMILAR** (middle finger), or **NEW** (no response). The response window is 2.0 s: the
1.5 s image followed by a 0.5 s blank during which responses are still accepted.

## Trial timing

| Event | 1-back / 2-back | MST |
|---|---|---|
| Fixation (jittered) | 0.75 ± 0.25 s (uniform, 0.5–1.0 s) | 0.75 ± 0.25 s |
| Image | 1.5 s | 1.5 s |
| Responsive blank | — | 0.5 s |
| Response window | 1.5 s (image) | 2.0 s (image + blank) |
| Mean trial | 2.25 s | 2.75 s |

## Run structure and scanner synchronization

Every run follows the same frame: instruction screen → experimenter starts the run →
the run **waits for the scanner's first trigger** ("5") → **6 s lead-in fixation** →
trials → **6 s tail fixation**. The 1-back additionally includes a single **6 s fixation
break near the block midpoint**, placed at a mini-block boundary so it never interrupts a
comparison. Lead-in, tail, and mid-run fixations provide a resting baseline and allow the
BOLD signal to settle before the first and after the last trial.

All within-run timing is locked to the trigger. Independently, the scanner TTL pulses are
routed to the EyeLink Host PC's parallel port and logged in the eye-tracking file as
INPUT events, giving a hardware record of every TR onset for offline synchronization.

Approximate active duration per run (expectation over the fixation jitter; excludes
instructions, between-block rests, and calibration): 1-back ≈ 7.0 min, 2-back ≈ 5.8 min,
MST ≈ 8.5 min. Part 1 totals ≈ 51 min and the full session ≈ 60 min of active task, plus
a 60 s rest after each 2-back run.

## Eye tracking

EyeLink 1000, 1000 Hz. A 9-point calibration is run at the start of the session, and an
optional recalibration is offered before the MST to correct drift accumulated over Part 1.
Fixation onsets, stimulus onsets, responses, and scanner triggers are written to per-run
EDF files.

## Data handling

Per-run behavioral tables (response and reaction time) are written immediately after each
run. The 1-back and 2-back data are consolidated and saved **before** the MST, so Part 1
is preserved even if Part 2 is interrupted. Eye-tracking files are closed on the Host PC
between runs and transferred in a single batch at the end of the session, so data
handling never intrudes on run timing.

---

*Placeholders for the reviewer: display and projector specifications, viewing distance,
button-box mapping at the console, and MRI acquisition parameters (scanner, TR, sequence)
are to be added.*
