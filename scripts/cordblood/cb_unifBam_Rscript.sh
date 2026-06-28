#!/bin/bash
#SBATCH --job-name=droplet_qc_ss
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=80G
#SBATCH --time=02-00:00:00
#SBATCH --output=/home/vreffo/scFound_bench/outputs/qc-%A.out
#SBATCH --error=/home/vreffo/scFound_bench/logs/qc-%A.err

set -euo pipefail

module purge
module load micromamba/rpy-sc
#module load samtools

# 1. Unisci tutti i BAM dei singoli SRR in un unico file master.bam
# (-@ 10 usa 10 core per velocizzare il processo)
#samtools merge -@ 10 /projects/shared/intronic_bam/cite_seq_coord_blood/starsolo_out/master_combined.bam /projects/shared/intronic_bam/cite_seq_coord_blood/starsolo_out/SRR580875*/Aligned.sortedByCoord.out.bam

# 2. Crea l'indice del file unito
#samtools index -@ 10 /projects/shared/intronic_bam/cite_seq_coord_blood/starsolo_out/master_combined.bam

Rscript $HOME/scFound_bench/r_scripts/cell_QC/qc_singlesample.R
