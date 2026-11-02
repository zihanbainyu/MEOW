#!/bin/bash
#SBATCH --account=torch_pr_1258_general
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=sltfce
#SBATCH --output=/scratch/zb2164/CBI/NeurWEM/analysis/code/sltfce_%j.log
singularity exec --bind /scratch:/scratch /share/apps/cbi/singularity/nipreps/fmriprep_25.2.1.sif python3 /scratch/zb2164/CBI/NeurWEM/analysis/code/sl_tfce.py 101 102 103 104 105 106 107 108 109
