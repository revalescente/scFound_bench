library(DropletQC)

# 1. Definisci i percorsi base costanti
gtf_path  <- "/projects/shared/intronic_bam/ref/ref/refdata-cellranger-hg19-3.0.0/genes/genes.gtf"
base_dir  <- "/projects/shared/intronic_bam/cite_seq_coord_blood/starsolo_out"
out_dir   <- "~/scFound_bench/results"

# Controllo difensivo: crea la cartella di output se non esiste
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 2. Recupera automaticamente la lista dei campioni (nomi delle sottocartelle SRR...)
samples <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

cat("Campioni individuati per l'analisi:\n")
print(samples)
cat("Totale campioni:", length(samples), "\n\n")

# 3. Ciclo su ogni singolo campione
for (sample in samples) {
  cat("\n" , paste(rep("=", 50), collapse = ""), "\n")
  cat("Inizio elaborazione campione:", sample, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")

  # Ricostruzione dinamica dei percorsi per il campione corrente
  bam_path <- file.path(base_dir, sample, "Aligned.sortedByCoord.out.bam")
  cb_path  <- file.path(base_dir, sample, "Solo.out/Gene/raw/barcodes.tsv")
  out_file <- file.path(out_dir, paste0("nf_", sample, ".rds"))

  # --- CONTROLLI DI SICUREZZA ---

  # 1. Salta il campione se il file di output esiste già (comodo se devi riprendere un job interrotto)
  if (file.exists(out_file)) {
    cat("Risultato già presente per", sample, "-> Salto al prossimo.\n")
    next
  }

  # 2. Verifica che il file BAM esista davvero
  if (!file.exists(bam_path)) {
    warning("ATTENZIONE: File BAM non trovato per ", sample, " all'indirizzo: ", bam_path)
    next
  }

  # 3. Verifica che il file dei barcode esista davvero
  if (!file.exists(cb_path)) {
    warning("ATTENZIONE: File barcodes.tsv non trovato per ", sample, " all'indirizzo: ", cb_path)
    next
  }

  # --- ESECUZIONE ---
  # Usiamo tryCatch per evitare che l'intero ciclo si blocchi se un singolo campione fallisce
  tryCatch({
    cat("Esecuzione di nuclear_fraction_annotation per", sample, "...\n")

    nf <- nuclear_fraction_annotation(
      annotation_path = gtf_path,
      bam = bam_path,
      barcodes = cb_path,
      tiles = 1000,
      cores = 10,
      verbose = TRUE
    )

    # Salvataggio dei risultati con il nome specifico del campione
    saveRDS(nf, file = out_file)
    cat("Campione", sample, "completato con successo! Salvato in:", out_file, "\n")

  }, error = function(e) {
    cat("ERRORE CRITICO durante l'elaborazione del campione", sample, ":\n")
    print(e)
  })
}

cat("\n========================================\n")
cat("Analisi di tutti i campioni completata!\n")
cat("========================================\n")
