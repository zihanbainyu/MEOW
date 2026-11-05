#!/usr/bin/env bash
# Refresh the group-level figure inputs from Torch, then regenerate every figure
# locally with figures/make_prelim_panels.m. Behaviour is already local (data/sub*).
# Usage:  bash sync_from_cluster.sh            # just pull the summaries + icons
#         bash sync_from_cluster.sh --extract  # re-run the cluster extraction first
set -e
HOST=zb2164@login.torch.hpc.nyu.edu
A=/scratch/zb2164/CBI/NeurWEM/analysis
SIF=/share/apps/cbi/singularity/nipreps/fmriprep_25.2.1.sif
HERE="$(cd "$(dirname "$0")" && pwd)"

if [ "$1" = "--extract" ]; then
  ssh "$HOST" "cd $A && singularity exec --bind /scratch:/scratch $SIF python3 -W ignore code/extract_2x2.py && singularity exec --bind /scratch:/scratch $SIF python3 -W ignore code/render_icons.py"
fi

mkdir -p "$HERE/data/mstback" "$HERE/data/mst" "$HERE/icons"
scp "$HOST:$A/group/hitmiss/roi_2x2_cells.csv" "$HERE/data/mstback/"
scp "$HOST:$A/group/mst/m3_roi.mat" "$HOST:$A/group/mst/mst_idx.mat" "$HOST:$A/group/mst/mst_behav.mat" "$HOST:$A/group/mst/sep_roi.mat" "$HERE/data/mst/"
scp "$HOST:$A/group/icons/roi_*.png" "$HERE/icons/"

echo "synced. now run:  matlab -batch \"run('scripts/mstback/make_prelim_panels.m')\"   (from $HERE)"
