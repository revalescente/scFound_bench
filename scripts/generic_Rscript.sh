#!/bin/bash
#SBATCH --job-name=rscript_exec
#SBATCH --partition=cpu
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --output=/home/%u/scFound_bench/outputs/rscript_exec-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/rscript_exec-%j.err

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# 1. Pulizia e caricamento del modulo
module purge
module load apptainer-bioconductor/3.23

# 2. Esecuzione dello script R
Rscript /home/vreffo/scFound_bench/r_scripts/cellxgene_QC/open_heart.R

echo "End time: $(date)"
