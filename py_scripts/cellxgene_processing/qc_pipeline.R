suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(anndataR)
  library(scrapper)
  library(DropletUtils)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Specificare il percorso del file .h5ad")
}

file_path <- args[1]
dataset_id <- tools::file_path_sans_ext(basename(file_path))
tissue_name <- basename(dirname(file_path))
out_dir <- file.path(dirname(file_path), "qc_results")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0(dataset_id, "_qc_metrics.rds"))

if (file.exists(out_file)) {
  cat(sprintf("[SKIP] %s già processato.\n", dataset_id))
  quit(save = "no", status = 0)
}

cat(sprintf("[%s] ELABORAZIONE: %s\n", tissue_name, dataset_id))

tryCatch({
  # 1. Lettura Dataset
  sce <- read_h5ad(file_path, as = "SingleCellExperiment")

  # Gestione nomi degli assay
  if (!"counts" %in% assayNames(sce)) {
    if ("X" %in% assayNames(sce)) {
      assay(sce, "counts") <- assay(sce, "X")
      assay(sce, "X") <- NULL
    } else if (length(assayNames(sce)) > 0) {
      assay(sce, "counts") <- assay(sce, assayNames(sce)[1])
    } else {
      stop("Nessuna matrice di conteggio trovata nell'H5AD")
    }
  }

  feature_names <- if ("feature_name" %in% colnames(rowData(sce))) rowData(sce)$feature_name else rownames(sce)
  is.mito <- grep("^MT-", feature_names, ignore.case = TRUE)

  # 2. QC Rapido
  sce <- quickRnaQc.se(sce, subsets = list(MT = is.mito))
  qc.thresh <- metadata(sce)$qc$thresholds
  mt_thresh <- if (!is.null(qc.thresh$subset.proportion["MT"])) qc.thresh$subset.proportion["MT"] else NA
  umi_summary <- summary(sce$sum)

  # 3. emptyDrops in sicurezza
  empty_pass_pct <- NA
  empty_out <- NULL

  tryCatch({
    lower_bound <- ifelse(!is.null(qc.thresh$sum), as.numeric(qc.thresh$sum), 100)
    e.out <- emptyDrops(counts(sce), lower = lower_bound)

    if (!is.null(e.out) && "FDR" %in% colnames(e.out)) {
      empty_pass_pct <- mean(e.out$FDR <= 0.001, na.rm = TRUE) * 100
      empty_out <- e.out  # 💾 Salva l'oggetto completo e.out
    }
  }, error = function(e) {
    empty_pass_pct <- NA
    empty_out <- NULL
  })

  # 4. MALAT1
  malat1_idx <- which(toupper(feature_names) == "MALAT1")
  malat1_stats <- if (length(malat1_idx) > 0) as.list(summary(counts(sce)[malat1_idx[1], ])) else list()

  # 5. Bias Tipi Cellulari
  cell_type_col <- if ("cell_type" %in% colnames(colData(sce))) "cell_type" else "cell_type_pred"

  if (cell_type_col %in% colnames(colData(sce))) {
    ct_df <- as.data.frame(table(CellType = colData(sce)[[cell_type_col]], Keep = sce$keep, useNA = "ifany")) %>%
      pivot_wider(names_from = Keep, values_from = Freq, values_fill = 0)

    if (!"TRUE" %in% colnames(ct_df)) ct_df[["TRUE"]] <- 0
    if (!"FALSE" %in% colnames(ct_df)) ct_df[["FALSE"]] <- 0

    ct_df <- ct_df %>%
      rename(Kept = `TRUE`, Discarded = `FALSE`) %>%
      mutate(
        Total = Kept + Discarded,
        Kept_Pct = round((Kept / Total) * 100, 2),
        Tissue = tissue_name,
        Dataset_ID = dataset_id
      )
  } else {
    ct_df <- data.frame(
      CellType = "Unknown",
      Kept = sum(sce$keep, na.rm = TRUE),
      Discarded = sum(!sce$keep, na.rm = TRUE),
      Total = ncol(sce),
      Kept_Pct = round(mean(sce$keep, na.rm = TRUE) * 100, 2),
      Tissue = tissue_name,
      Dataset_ID = dataset_id
    )
  }

  # 6. Salvataggio Output
  res <- list(
        summary = data.frame(
        Tissue = tissue_name,
        Dataset_ID = dataset_id,
        Total_Cells = ncol(sce),
        Kept_Cells = sum(sce$keep, na.rm = TRUE),
        Discarded_Cells = sum(!sce$keep, na.rm = TRUE),
        Kept_Pct = round(mean(sce$keep, na.rm = TRUE) * 100, 2),
        MT_Threshold = round(as.numeric(mt_thresh), 4),
        EmptyDrops_Pass_Pct = round(as.numeric(empty_pass_pct), 2),
        Min_UMI = as.numeric(umi_summary["Min."]),
        Median_UMI = as.numeric(umi_summary["Median"]),
        Mean_UMI = as.numeric(umi_summary["Mean"]),
        Max_UMI = as.numeric(umi_summary["Max."]),
        stringsAsFactors = FALSE
        ),
        celltype_bias = ct_df,
        malat1_stats = malat1_stats,
        emptydrops_res = empty_out  # 👈 Oggetto completo (o NULL se fallisce/NA)
  )

  saveRDS(res, out_file)
  cat(sprintf("[SUCCESS] %s salvato.\n", dataset_id))

}, error = function(e) {
  cat(sprintf("[ERROR] %s: %s\n", dataset_id, e$message))
  quit(save = "no", status = 1)
})
