#!/bin/bash
#SBATCH --job-name=droplet_qc_ss
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=80G
#SBATCH --time=02-00:00:00
#SBATCH --output=home/%u/scFound_bench/outputs/qc-%A.out
#SBATCH --error=home/%u/scFound_bench/logs/qc-%A.err

set -euo pipefail

module purge
module load micromamba/rpy-sc

Rscript $HOME/scFound_bench/r_scripts/cell_QC/qc_allsamples.R
