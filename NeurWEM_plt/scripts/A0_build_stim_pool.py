#!/usr/bin/env python3
"""
Build a locked stimulus pool for the Hybrid MST N-Back task.

Rationale
---------
A_subject_setup.m used to load the full master pool (180 L1 + 180 L2 pairs)
and, per subject, randomly split each level into thirds -- comp / novel /
UNUSED -- with a clock-seeded RNG. Consequently every subject saw a different,
non-reproducible 2/3 subset of the pairs, and a full third of each level was
never presented to anyone.

This script locks the pool ONCE: it selects a fixed subset of pairs (the exact
number the task actually uses) with a fixed seed, copies them into a new folder
that becomes the pool, and records the selection in a manifest. A_subject_setup.m
then only shuffles condition assignment WITHIN this fixed set, so every subject
sees the same items and the selection is auditable.

What gets locked
----------------
  - L1 pairs: KEEP_PER_LEVEL of the available L1 pairs (A + B)
  - L2 pairs: KEEP_PER_LEVEL of the available L2 pairs (A + B)
Foils, practice pairs, and instruction images are copied through in full: foils
are trial-unique novel items (drawn randomly per subject, variable count due to
2-back padding), and prac_/instr_ are fixed support assets the run scripts load.

Run once from the scripts/ directory:  python3 A0_build_stim_pool.py
Re-running with the same SEED reproduces the identical pool.
"""

import re
import shutil
import random
from pathlib import Path

# ----------------------------------------------------------------------------
SEED = 20260727            # fixed -> reproducible selection; do not change
KEEP_PER_LEVEL = 120       # pairs kept per level (60 compared + 60 novel)
SRC = Path("../stimulus/stim_final")
DST = Path("../stimulus/stim_pool")
# ----------------------------------------------------------------------------


def pair_ids(level):
    """Return sorted integer ids of the A files for a given level ('l1'/'l2')."""
    pat = re.compile(rf"mst_0*(\d+)_A_{level}\.png$")
    ids = []
    for f in SRC.glob(f"mst_*_A_{level}.png"):
        m = pat.match(f.name)
        if m:
            ids.append(int(m.group(1)))
    return sorted(ids)


def copy_pair(stem_id, level):
    """Copy the A and B image of one pair into DST. stem preserves zero-pad."""
    # recover the original zero-padded stem from the A filename on disk
    matches = list(SRC.glob(f"mst_*_A_{level}.png"))
    for f in matches:
        m = re.match(rf"mst_0*{stem_id}_A_{level}\.png$", f.name)
        if m and int(re.match(rf"mst_0*(\d+)_A_{level}\.png$", f.name).group(1)) == stem_id:
            stem = f.name[:-len(f"_A_{level}.png")]  # e.g. 'mst_017'
            shutil.copy2(SRC / f"{stem}_A_{level}.png", DST / f"{stem}_A_{level}.png")
            shutil.copy2(SRC / f"{stem}_B_{level}.png", DST / f"{stem}_B_{level}.png")
            return stem
    raise FileNotFoundError(f"pair {stem_id} ({level}) not found")


def main():
    assert SRC.is_dir(), f"source not found: {SRC.resolve()}"
    DST.mkdir(parents=True, exist_ok=True)

    rng = random.Random(SEED)
    manifest_lines = [f"# stim_pool manifest  seed={SEED}  keep_per_level={KEEP_PER_LEVEL}",
                      "level,id,status"]

    for level in ("l1", "l2"):
        ids = pair_ids(level)
        assert len(ids) >= KEEP_PER_LEVEL, \
            f"{level}: have {len(ids)} pairs, need {KEEP_PER_LEVEL}"
        shuffled = ids[:]
        rng.shuffle(shuffled)
        kept = sorted(shuffled[:KEEP_PER_LEVEL])
        dropped = sorted(shuffled[KEEP_PER_LEVEL:])
        for i in kept:
            copy_pair(i, level)
            manifest_lines.append(f"{level},{i},kept")
        for i in dropped:
            manifest_lines.append(f"{level},{i},dropped")
        print(f"  {level}: kept {len(kept)}, dropped {len(dropped)}")

    # copy support assets in full (foils, practice pairs, instructions)
    n_support = 0
    for f in SRC.iterdir():
        if f.is_file() and f.suffix == ".png" and (
            "_foil" in f.name or f.name.startswith("prac_") or f.name.startswith("instr_")
        ):
            shutil.copy2(f, DST / f.name)
            n_support += 1
    print(f"  support assets copied: {n_support}")

    (DST / "POOL_MANIFEST.csv").write_text("\n".join(manifest_lines) + "\n")
    print(f"Pool built at {DST.resolve()}  (manifest: POOL_MANIFEST.csv)")


if __name__ == "__main__":
    main()
