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
  
  # Se i rownames sono NULL, assegnali da obs_names o rownames(ann$obs)
  if (is.null(rownames(scgpt_matrix))) {
    if (!is.null(ann$obs_names)) {
      rownames(scgpt_matrix) <- as.character(ann$obs_names)
    } else if (!is.null(rownames(ann$obs))) {
      rownames(scgpt_matrix) <- rownames(ann$obs)
    }
  }
  
  # Rimuove la stringa finale "-0" dai rownames degli embeddings
  if (!is.null(rownames(scgpt_matrix))) {
    rownames(scgpt_matrix) <- sub("-0$", "", rownames(scgpt_matrix))
  }
  
  # Identifica le cellule in comune
  common_cells <- intersect(colnames(sce), rownames(scgpt_matrix))
  if (length(common_cells) == 0) {
    warning(paste("Nessun barcode corrispondente trovato per:", embed_name))
    return(sce)
  }
  
  # Mantiene solo le cellule/barcode che matchano in entrambi gli oggetti
  sce <- sce[, common_cells]
  reducedDim(sce, embed_name) <- scgpt_matrix[common_cells, , drop = FALSE]
  
  return(sce)
}

# Caricatore sicuro di embedding CSV (GeneFormer)
add_csv_embedding <- function(sce, embed_name, path, barcode_col = "Barcode") {
  if (!file.exists(path)) return(sce)
  emb <- readr::read_csv(path, show_col_types = FALSE)
  
  # Pulizia barcode: rimuove il trattino e le cifre finali (es. -0, -1, -99)
  emb[[barcode_col]] <- sub("-[0-9]+$", "", as.character(emb[[barcode_col]]))
  
  # Identifica le cellule in comune
  common_cells <- intersect(colnames(sce), emb[[barcode_col]])
  if (length(common_cells) == 0) {
    warning(paste("Nessun barcode corrispondente trovato per:", embed_name))
    return(sce)
  }
  
  # Filtra ed elimina eventuali duplicati di barcode
  emb_filtered <- emb |>
    dplyr::filter(.data[[barcode_col]] %in% common_cells) |>
    dplyr::distinct(.data[[barcode_col]], .keep_all = TRUE)
  
  # Estrazione matrice numerica
  mat <- emb_filtered |>
    dplyr::select(where(is.numeric), -dplyr::any_of(c("Unnamed: 0", "...1"))) |>
    as.matrix()
  
  rownames(mat) <- emb_filtered[[barcode_col]]
  
  # Sincronizza cellule e assegna l'embedding ordinato
  sce <- sce[, common_cells]
  reducedDim(sce, embed_name) <- mat[common_cells, , drop = FALSE]
  
  return(sce)
}
