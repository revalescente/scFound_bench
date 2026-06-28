#!/bin/bash
#SBATCH --job-name=dropletqc_array
#SBATCH --partition=cpu
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/home/vreffo/scFound_bench/outputs/be1qc-%A_%a.out
#SBATCH --error=/home/vreffo/scFound_bench/logs/be1qc-%A_%a.err

set -euo pipefail

# Carica il modulo R presente sul tuo cluster
module purge
module load micromamba/rpy-sc

# Array contenente i nomi esatti dei campioni (l'ordine mappa l'ID dell'array 0-7)
SAMPLES=("A549" "CCL-185-IG" "PC9" "CRL5868" "HTB178" "DV90" "HCC78" "PBMCs")

# Seleziona il campione corrente usando il TASK ID di Slurm
CURRENT_SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "====================================================="
echo "Array Job ID : ${SLURM_ARRAY_JOB_ID}"
echo "Task ID      : ${SLURM_ARRAY_TASK_ID}"
echo "Campione     : ${CURRENT_SAMPLE}"
echo "Core allocati: ${SLURM_CPUS_PER_TASK}"
echo "Start        : $(date)"
echo "====================================================="

# Esegui lo script R passando il nome del campione come argomento
Rscript $HOME/scFound_bench/r_scripts/cell_QC/BE1/qc_singlesample.R "${CURRENT_SAMPLE}"

echo "End: $(date)"
echo "-----------------------------------------------------"
