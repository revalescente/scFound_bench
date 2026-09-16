library(SingleCellExperiment)
library(anndataR)
library(purrr)
library(future)
library(furrr)

# 1. Configurazione del parallelismo sul cluster
# Usa 8-16 worker: sufficienti a velocizzare senza saturare il disco in scrittura
num_workers <- 16
plan(multisession, workers = num_workers)

data_dir <- "/projects/shared/cellxgene/heart"
output_dir <- "/projects/shared/cellxgene/heart/split_datasets"

# 2. Lettura PARALLELA dei file .h5ad
files <- list.files(data_dir, pattern = "\\.h5ad$", full.names = TRUE)

message(sprintf("⚡ Lettura parallela di %d file .h5ad con %d worker...", length(files), num_workers))

sce_all <- future_map(
  files[c(1,2)],
  \(file) anndataR::read_h5ad(file, as = "SingleCellExperiment"),
  .options = furrr_options(seed = TRUE)
) |> reduce(cbind)

# 3. Splitting per dataset_id
col_indices <- split.default(seq_len(ncol(sce_all)), sce_all$dataset_id)
sce_list <- map(col_indices, ~ sce_all[, .x])

# Convalida
counts_orig <- table(sce_all$dataset_id)
counts_split <- map_int(sce_list, ncol)[names(counts_orig)]

rm(sce_all, col_indices)
gc()

# 4. Salvataggio PARALLELO dei file RDS
if (identical(as.numeric(counts_orig), as.numeric(counts_split))) {

  message(sprintf("✅ CONVALIDA RIUSCITA: %d dataset pronti per il salvataggio.", length(sce_list)))

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("💾 Salvataggio parallelo dei file RDS in corso...")

  future_iwalk(
    sce_list,
    \(sce_obj, dataset_id) {
      file_path <- file.path(output_dir, paste0(dataset_id, ".rds"))

      # compress = FALSE per la massima velocità di scrittura (senza compressione)
      # usa compress = TRUE (o "gzip") se vuoi risparmiare spazio su disco
      saveRDS(sce_obj, file = file_path, compress = FALSE)
    },
    .options = furrr_options(
      seed = TRUE,
      packages = c("SingleCellExperiment", "SummarizedExperiment") # Risolve il warning
    )
  )

  message("🎉 Salvataggio completato con successo!")

} else {
  warning("⚠️ ERRORE DI CONVALIDA: Salvataggio annullato.")
}

# Ripristina l'esecuzione sequenziale standard
plan(sequential)
