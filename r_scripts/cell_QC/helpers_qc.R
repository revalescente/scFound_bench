library(dplyr)
library(purrr)
library(SingleCellExperiment)

# Caricamento dinamico della Nuclear Fraction
fetch_nuclear_fraction <- function(sce, dataset_name, results_dir = "~/scFound_bench/results") {
  single_path <- file.path(results_dir, paste0("nf_", dataset_name, ".rds"))
  multi_files <- list.files(results_dir, pattern = paste0("^nf_", dataset_name, "_.*\\.rds$"), full.names = TRUE)

  if (file.exists(single_path)) {
    nf <- readRDS(single_path)
    return(nf[colnames(sce), "nuclear_fraction", drop = TRUE])
  }

  if (length(multi_files) > 0) {
    nf_combined <- purrr::map_dfr(multi_files, function(f) {
      s_name <- sub(paste0(".*nf_", dataset_name, "_(.*)\\.rds$"), "\\1", f)
      df <- as.data.frame(readRDS(f))
      df$cell_key <- paste0(rownames(df), "_", s_name)
      df
    })
    return(nf_combined$nuclear_fraction[match(colnames(sce), nf_combined$cell_key)])
  }

  # Fallback: calcolo diretto
  total_reads <- rowSums(
    as.data.frame(colData(sce)[, c("unaligned", "aligned_unmapped", "mapped_to_exon",
                                  "mapped_to_intron", "ambiguous_mapping",
                                  "mapped_to_MT", "mapped_to_ERCC"), drop = FALSE]),
    na.rm = TRUE
  )
  return(sce$mapped_to_intron / total_reads)
}

# Caricatore sicuro di embedding h5ad (scGPT)
add_h5ad_embedding <- function(sce, embed_name, path, obsm_key) {
  if (!file.exists(path)) return(sce)
  ann <- anndataR::read_h5ad(path)
  scgpt_matrix <- as.matrix(ann$obsm[obsm_key][[obsm_key]])
  reducedDim(sce, embed_name) <- scgpt_matrix[colnames(sce), ]
  return(sce)
}

# Caricatore sicuro di embedding CSV (GeneFormer)
add_csv_embedding <- function(sce, embed_name, path, barcode_col = "Barcode") {
  if (!file.exists(path)) return(sce)
  emb <- readr::read_csv(path, show_col_types = FALSE)

  mat <- emb %>%
    filter(as.character(.data[[barcode_col]]) %in% colnames(sce)) %>%
    select(where(is.numeric), -any_of(c("Unnamed: 0", "...1"))) %>%
    as.matrix()

  rownames(mat) <- as.character(emb[[barcode_col]][emb[[barcode_col]] %in% colnames(sce)])
  reducedDim(sce, embed_name) <- mat[colnames(sce), ]
  return(sce)
}
