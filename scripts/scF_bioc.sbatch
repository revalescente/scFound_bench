#!/bin/bash
#SBATCH --job-name=scF_bioc
#SBATCH --partition=cpu
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --output=/home/%u/scFound_bench/outputs/slurm-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/slurm-%j.err

# specify the paths of the out and err to have them created in a sensible position

set -euo pipefail

# Specify the PORT number for rstudio
RSTUDIO_PORT="${PORT:-8799}"

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Starting RStudio Server on port: $RSTUDIO_PORT"

# Load the required environment
module purge
module load apptainer-bioconductor/3.23

# Start rstudio server on the compute node
rstudio --ip 0.0.0.0 --port "$RSTUDIO_PORT"

# How to call it: PORT=8888 sbatch --time=00:10:00 --ntasks=1 --cpus-per-task=10 --mem=1G \
# --output=dirpath/slurm-%j.out --error=dirpath/slurm-%j.err filename.sbatch


module load apptainer-bioconductor/3.23

Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("DropletUtils", "scrapper"), ask = FALSE)'
