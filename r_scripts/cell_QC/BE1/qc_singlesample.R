library(DropletQC)

# 1. Recupera il campione passato come argomento da Bash
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("❌ ERRORE: Nessun campione specificato! Uso: Rscript run_dropletqc_single.R <sample_name>")
}
sample <- args[1]

# 2. Configurazione Percorsi Base (Niente più GTF!)
multi_out <- "/projects/shared/intronic_bam/BE1/cellranger_multi_out/BE1_demux/outs"
out_dir   <- "~/scFound_bench/results"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("\n==================================================\n")
cat(" 🚀 ANALISI DROPQC TRAMITE CELLRANGER TAGS:", sample, "\n")
cat("==================================================\n")

# 3. Definisci i percorsi specifici per questo campione
sample_dir   <- file.path(multi_out, "per_sample_outs", sample, "count")
bam_path     <- file.path(sample_dir, "sample_alignments.bam")
barcode_path <- file.path(sample_dir, "sample_filtered_feature_bc_matrix", "barcodes.tsv.gz")
out_file     <- file.path(out_dir, paste0("nf_BE1_", sample, ".rds"))

# 4. Controlli di sicurezza
if (!file.exists(bam_path)) {
  stop("❌ BAM non trovato al percorso: ", bam_path)
}
if (!file.exists(barcode_path)) {
  stop("❌ Barcodes non trovati al percorso: ", barcode_path)
}

# 5. Esecuzione DropletQC sfruttando i TAG di Cell Ranger
slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))

cat("-> Calcolo della frazione nucleare ultra-rapido (con", slurm_cores, "core)...\n")
nf <- nuclear_fraction_tags(
  bam = bam_path,
  barcodes = barcode_path, # Accetta direttamente il path del file .tsv.gz
  tiles = 1000,
  cores = slurm_cores,
  verbose = TRUE
)

# 6. Salvataggio
saveRDS(nf, file = out_file)
cat("✅ Calcolo completato! File salvato in:", out_file, "\n")
