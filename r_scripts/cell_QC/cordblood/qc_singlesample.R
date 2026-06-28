library(DropletQC)

# 1. Definisci i percorsi
gtf_path  <- "/projects/shared/intronic_bam/ref/ref/refdata-cellranger-hg19-3.0.0/genes/genes.gtf"
bam_path  <- "/projects/shared/intronic_bam/cite_seq_coord_blood/starsolo_out/master_combined.bam"
cb_path   <- "/projects/shared/intronic_bam/cite_seq_coord_blood/starsolo_out/SRR5808750/Solo.out/Gene/raw/barcodes.tsv"
out_dir   <- "~/scFound_bench/results"
out_file  <- file.path(out_dir, "nf_cordblood.rds")

# 2. Controllo difensivo: crea la cartella di output se non esiste
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 3. Esecuzione di DropletQC
nf <- nuclear_fraction_annotation(
  annotation_path = gtf_path,
  bam = bam_path,
  barcodes = cb_path,
  tiles = 1000,
  cores = 10,
  verbose = TRUE
)

# 4. Salvataggio dei risultati
saveRDS(nf, file = out_file)
cat("Analisi completata con successo! File salvato in:", out_file, "\n")
