# ================================================================
#  BENCHMARK: SingleR vs scGPT vs Geneformer — Dataset BE1
#
#  Ground truth  : colData(sce)$Sample  (annotazione manuale)
#  Metodi testati: SingleR (Bioconductor), scGPT, Geneformer
#  Split         : 70% train (reference SingleR) | 30% test
#  Metriche      : Accuracy, F1 macro/weighted, ARI,
#                  F1 per classe, confusion matrix
# ================================================================

# ── Librerie ────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(SingleR)
  library(bluster)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(mclust)      # ARI
  library(caret)       # confusionMatrix, F1
  library(RColorBrewer)
})

set.seed(42)
BPPARAM <- MulticoreParam(workers = 4)

# ================================================================
#  0. CONFIGURAZIONE  ← modifica qui
# ================================================================
SCE_PATH      <- "/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS"
GT_COL        <- "Sample"          # colonna ground truth nel colData

# CSV con predizioni FM: devono avere almeno due colonne:
#   "cell_id"  → uguale ai colnames del SCE
#   "celltype" → etichetta predetta

# scgpt prediction results is in the h5ad: /projects/shared/intronic_bam/datasets/anndata/be1_scGPT_embeddings.h5ad
SCGPT_CSV     <- "/projects/shared/intronic_bam/datasets/anndata/scgpt_predictions.csv"
GENEFORMER_CSV<- "geneformer_predictions.csv"

OUT_DIR       <- "predicts/output_benchmark"
#TRAIN_PROP    <- 0.70              # proporzione training set

# ================================================================
#  1. CARICAMENTO E ISPEZIONE
# ================================================================
cat("=== 1. Caricamento ===\n")
sce <- readRDS(SCE_PATH)
cat("Dimensioni:", nrow(sce), "geni ×", ncol(sce), "cellule\n")

stopifnot(GT_COL %in% colnames(colData(sce)))
sce$gt <- factor(colData(sce)[[GT_COL]])
cat("Tipi cellulari (ground truth):", nlevels(sce$gt), "\n")
print(sort(table(sce$gt), decreasing = TRUE))

# ================================================================
#  2. PRE-PROCESSING (skip se già eseguito)
# ================================================================
# cat("\n=== 2. Pre-processing ===\n")
# 
# if (!"logcounts" %in% assayNames(sce)) {
#   cat("  Normalizzazione scran...\n")
#   qclust <- quickCluster(sce, BPPARAM = BPPARAM)
#   sce    <- computeSumFactors(sce, clusters = qclust, BPPARAM = BPPARAM)
#   sce    <- logNormCounts(sce)
# } else { cat("  logcounts già presenti.\n") }
# 
# if (!"PCA" %in% reducedDimNames(sce)) {
#   cat("  HVG + PCA...\n")
#   dec  <- modelGeneVar(sce, BPPARAM = BPPARAM)
#   hvgs <- getTopHVGs(dec, prop = 0.10)
#   sce  <- runPCA(sce, subset_row = hvgs, ncomponents = 50, BPPARAM = BPPARAM)
# }
# if (!"UMAP" %in% reducedDimNames(sce)) {
#   cat("  UMAP...\n")
#   sce <- runUMAP(sce, dimred = "PCA", n_dimred = 30, BPPARAM = BPPARAM)
# }

# ================================================================
#  3. SPLIT TRAIN / TEST — stratificato per tipo cellulare
# ================================================================
cat("\n=== 3. Split train/test (stratificato) ===\n")

# train_idx <- unlist(lapply(levels(sce$gt), function(ct) {
#   idx    <- which(sce$gt == ct)
#   n_tr   <- max(1L, round(length(idx) * TRAIN_PROP))
#   sample(idx, n_tr)
# }))
# test_idx  <- setdiff(seq_len(ncol(sce)), train_idx)
# 
# sce$split <- "test"
# sce$split[train_idx] <- "train"
# 
# sce_train <- sce[, train_idx]
# sce_test  <- sce[, test_idx]

# Uso train e test definito da scGPT in python
table(sce$is_train)
sce_train <- sce[, sce$is_train == TRUE]
sce_test <- sce[,  sce$is_train != TRUE]

cat(sprintf("  Train: %d cellule | Test: %d cellule\n",
            ncol(sce_train), ncol(sce_test)))

# Verifica bilanciamento classi nel test set
cat("  Distribuzione test set:\n")
print(prop.table(table(sce_test$gt)))
print(prop.table(table(sce_train$gt)))

# ================================================================
#  4. SINGLER — reference da training set (etichette manuali)
# ================================================================
cat("\n=== 4. SingleR ===\n")

# Reference pseudo-bulk dal training con ground truth manuali
ref_singler <- aggregateAcrossCells(
  sce_train,
  ids            = sce_train$gt,
  use.assay.type = "counts",
  statistics     = "sum"
)
ref_singler$sizeFactor <- NULL
ref_singler <- logNormCounts(ref_singler)

cat("  Reference costruito:", ncol(ref_singler), "tipi cellulari\n")

pred_singler <- SingleR(
  test    = sce_test,
  ref     = ref_singler,
  labels  = ref_singler$gt,
  BPPARAM = BPPARAM
)

singler_labels <- pred_singler$labels          # etichette finali
singler_scores <- apply(pred_singler$scores, 1, max)
singler_pruned <- pred_singler$pruned.labels   # NA = bassa confidenza

sce$singler_labels <- NA_character_
sce$singler_scores <- NA_real_
sce$singler_pruned  <- NA_character_

sce$singler_labels[sce$is_train == FALSE] <- singler_labels
sce$singler_scores[sce$is_train == FALSE] <- singler_scores
sce$singler_pruned[sce$is_train == FALSE]  <- singler_pruned

cat("  Cellule rifiutate (pruning):",
    sum(is.na(singler_pruned)), "/", ncol(sce_test), "\n")

# ================================================================
#  5. CARICAMENTO PREDIZIONI FOUNDATION MODELS
# ================================================================
cat("\n=== 5. Predizioni FM ===\n")

load_fm_preds <- function(path, method_name, test_cell_ids, id_col = "Barcode") {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[match(test_cell_ids, df[[id_col]]), ]
  df$celltype
}

gt_test   <- as.character(sce_test$gt)
cell_ids  <- colnames(sce_test)

scgpt_labels      <- load_fm_preds(SCGPT_CSV,      "scGPT",       cell_ids)
geneformer_labels <- sce$pred_GeneFormer[!is.na(sce$pred_GeneFormer)]

cat("  scGPT caricato:      ", length(scgpt_labels),      "predizioni\n")
cat("  Geneformer caricato: ", length(geneformer_labels), "predizioni\n")

# ================================================================
#  6. ARMONIZZAZIONE ETICHETTE
#     → le etichette dei FM potrebbero avere piccole differenze
#        di capitalizzazione/spazi rispetto al ground truth
# ================================================================
harmonize <- function(x, ref_labels) {
  x <- trimws(x)
  # Normalizza maiuscole: cerca match case-insensitive
  lookup <- setNames(ref_labels, tolower(ref_labels))
  x_low  <- tolower(x)
  mapped <- lookup[x_low]
  # Dove non trovato, mantieni l'originale (comparazione fallirà)
  ifelse(is.na(mapped), x, mapped)
}

ref_labels        <- levels(sce$gt)
singler_labels    <- harmonize(singler_labels,    ref_labels)
scgpt_labels      <- harmonize(scgpt_labels,      ref_labels)
geneformer_labels <- harmonize(geneformer_labels, ref_labels)

# ================================================================
#  7. METRICHE PER METODO
# ================================================================
cat("\n=== 6. Calcolo metriche ===\n")

all_labels <- sort(unique(c(gt_test, singler_labels,
                             scgpt_labels, geneformer_labels)))


compute_metrics <- function(pred, truth, method, all_lev) {
  pred  <- factor(pred,  levels = all_lev)
  truth <- factor(truth, levels = all_lev)
  
  # Accuracy globale
  acc <- mean(pred == truth, na.rm = TRUE)
  
  # Confusion matrix via caret
  cm  <- confusionMatrix(pred, truth)
  by_class <- as.data.frame(cm$byClass)   # per classe
  
  # F1 macro (media non pesata su classi)
  f1_per_class <- by_class$F1
  f1_macro     <- mean(f1_per_class, na.rm = TRUE)
  
  # F1 weighted (pesato per frequenza classe)
  class_counts <- table(truth)
  weights      <- as.numeric(class_counts) / sum(class_counts)
  f1_weighted  <- sum(f1_per_class * weights, na.rm = TRUE)
  
  # ARI
  ari <- adjustedRandIndex(as.integer(pred), as.integer(truth))
  
  list(
    method      = method,
    accuracy    = acc,
    f1_macro    = f1_macro,
    f1_weighted = f1_weighted,
    ari         = ari,
    f1_per_class= setNames(f1_per_class,
                           gsub("Class: ", "", rownames(by_class))),
    cm          = cm$table
  )
}

metrics <- list(
  SingleR    = compute_metrics(singler_labels,    gt_test, "SingleR",    all_labels),
  scGPT      = compute_metrics(scgpt_labels,      gt_test, "scGPT",      all_labels),
  Geneformer = compute_metrics(geneformer_labels, gt_test, "Geneformer", all_labels)
)


# Tabella riepilogativa
summary_df <- do.call(rbind, lapply(metrics, function(m) {
  data.frame(
    Metodo      = m$method,
    Accuracy    = round(m$accuracy,    4),
    F1_macro    = round(m$f1_macro,    4),
    F1_weighted = round(m$f1_weighted, 4),
    ARI         = round(m$ari,         4)
  )
})) %>% arrange(desc(F1_macro))

cat("\n┌─── RIEPILOGO METRICHE ────────────────────────────────┐\n")
print(summary_df, row.names = FALSE)
cat("└───────────────────────────────────────────────────────┘\n")

# F1 per classe (tabella lunga)
f1_class_df <- do.call(rbind, lapply(metrics, function(m) {
  data.frame(
    Metodo   = m$method,
    CellType = names(m$f1_per_class),
    F1       = m$f1_per_class
  )
})) %>% filter(!is.na(F1))

# ================================================================
#  8. VISUALIZZAZIONI
# ================================================================
cat("\n=== 7. Visualizzazioni ===\n")
dir.create(OUT_DIR, showWarnings = FALSE)

METHOD_COLS <- c(SingleR    = "#4e79a7",
                 scGPT      = "#59a14f",
                 Geneformer = "#e15759")


## 8a. Barplot metriche globali
df_glob <- summary_df %>%
  pivot_longer(cols = c(Accuracy, F1_macro, F1_weighted, ARI),
               names_to = "Metrica", values_to = "Valore")

p_global <- ggplot(df_glob, aes(x = Metrica, y = Valore, fill = Metodo)) +
  geom_col(position = "dodge", alpha = 0.88, width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", Valore)),
            position = position_dodge(0.65), vjust = -0.4, size = 3) +
  scale_fill_manual(values = METHOD_COLS) +
  scale_y_continuous(limits = c(0, 1.05), labels = percent_format()) +
  labs(title    = "Benchmark BE1 — metriche globali",
       subtitle = "Test set (30%); ground truth = annotazione manuale",
       x = NULL, y = NULL, fill = "Metodo") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "1_metriche_globali.png"),
       p_global, width = 9, height = 6)

## 8b. F1 per classe — dot + segment chart
p_f1class <- ggplot(f1_class_df,
                    aes(x = F1, y = reorder(CellType, F1, mean),
                        colour = Metodo)) +
  geom_line(aes(group = CellType), colour = "grey80", linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = METHOD_COLS) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
  labs(title  = "F1 per tipo cellulare",
       x = "F1-score", y = NULL, colour = "Metodo") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "2_f1_per_celltype.png"),
       p_f1class, width = 9, height = max(6, nlevels(sce$gt) * 0.4 + 3))

## 8c. Confusion matrix normalizzata per ognuno dei tre metodi
plot_cm <- function(cm_table, method_name) {
  mat      <- as.matrix(cm_table)
  mat_norm <- sweep(mat, 2, colSums(mat), "/")  # normalizza per classe vera
  df_cm    <- as.data.frame(as.table(mat_norm))
  colnames(df_cm) <- c("Predetto", "Vero", "Proporzione")
  
  ggplot(df_cm, aes(x = Vero, y = Predetto, fill = Proporzione)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(Proporzione > 0.02,
                                 sprintf("%.2f", Proporzione), "")),
              size = 2.5, colour = "white") +
    scale_fill_gradientn(
      colours = c("#f7fbff", "#6baed6", "#2171b5", "#08306b"),
      limits  = c(0, 1), labels = percent_format()
    ) +
    labs(title = paste("Confusion matrix —", method_name),
         x = "Tipo reale", y = "Tipo predetto", fill = "Proporzione") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

for (m in names(metrics)) {
  p <- plot_cm(metrics[[m]]$cm, m)
  ggsave(file.path(OUT_DIR, paste0("3_confusion_", tolower(m), ".png")),
         p, width = max(8, nlevels(sce$gt) * 0.6),
         height = max(7, nlevels(sce$gt) * 0.55))
}

## 8d. UMAP con predizioni dei tre metodi affiancate
sce_test$pred_singler    <- singler_labels
sce_test$pred_scgpt      <- scgpt_labels
sce_test$pred_geneformer <- geneformer_labels

p_gt  <- plotReducedDim(sce_test, "UMAP", colour_by = "gt",
                        point_size = 0.4, point_alpha = 0.6) +
  ggtitle("Ground truth") + theme_bw(base_size = 9) +
  guides(colour = guide_legend(override.aes = list(size = 2.5),
                               ncol = 1, title = NULL))

p_sr  <- plotReducedDim(sce_test, "UMAP", colour_by = "pred_singler",
                        point_size = 0.4, point_alpha = 0.6) +
  ggtitle("SingleR") + theme_bw(base_size = 9) +
  guides(colour = guide_legend(override.aes = list(size = 2.5),
                               ncol = 1, title = NULL))

p_sg  <- plotReducedDim(sce_test, "UMAP", colour_by = "pred_scgpt",
                        point_size = 0.4, point_alpha = 0.6) +
  ggtitle("scGPT") + theme_bw(base_size = 9) +
  guides(colour = guide_legend(override.aes = list(size = 2.5),
                               ncol = 1, title = NULL))

p_gf  <- plotReducedDim(sce_test, "UMAP", colour_by = "pred_geneformer",
                        point_size = 0.4, point_alpha = 0.6) +
  ggtitle("Geneformer") + theme_bw(base_size = 9) +
  guides(colour = guide_legend(override.aes = list(size = 2.5),
                               ncol = 1, title = NULL))

umap_panel <- ( p_gt | p_sr) / (p_sg | p_gf)
ggsave(file.path(OUT_DIR, "4_umap_confronto.png"),
       umap_panel, width = 16, height = 14)

## 8e. Heatmap F1 per classe × metodo
f1_mat <- f1_class_df %>%
  pivot_wider(names_from = Metodo, values_from = F1) %>%
  textshape::column_to_rownames("CellType") %>%
  as.matrix()

pheatmap(f1_mat,
         cluster_rows = TRUE, cluster_cols = FALSE,
         color    = colorRampPalette(c("#d73027", "#fee090", "#1a9850"))(50),
         breaks   = seq(0, 1, length.out = 51),
         display_numbers = TRUE, number_format = "%.2f",
         number_color = "black", fontsize_number = 8,
         main     = "F1-score per tipo cellulare e metodo",
         filename = file.path(OUT_DIR, "5_heatmap_f1_per_celltype.png"),
         width    = 7, height = max(5, nrow(f1_mat) * 0.4 + 3))

## 8f. Radar chart errori per cellula (errore = predizione ≠ ground truth)
error_df <- data.frame(
  CellType   = gt_test,
  SingleR    = as.integer(singler_labels    != gt_test),
  scGPT      = as.integer(scgpt_labels      != gt_test),
  Geneformer = as.integer(geneformer_labels != gt_test)
) %>%
  group_by(CellType) %>%
  summarise(across(c(SingleR, scGPT, Geneformer), mean), .groups = "drop") %>%
  pivot_longer(-CellType, names_to = "Metodo", values_to = "ErrorRate")

p_err <- ggplot(error_df,
                aes(x = reorder(CellType, ErrorRate, mean),
                    y = ErrorRate, fill = Metodo)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = METHOD_COLS) +
  scale_y_continuous(labels = percent_format()) +
  coord_flip() +
  labs(title = "Tasso di errore per tipo cellulare",
       x = NULL, y = "Errore (%)", fill = "Metodo") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "6_error_rate_per_celltype.png"),
       p_err, width = 10, height = max(6, length(unique(gt_test)) * 0.5 + 3))

cat("  Grafici salvati in", OUT_DIR, "\n")

# ================================================================
#  9. SALVATAGGIO RISULTATI
# ================================================================
cat("\n=== 8. Salvataggio ===\n")

write.csv(summary_df,
          file.path(OUT_DIR, "metriche_globali.csv"),
          row.names = FALSE)

write.csv(f1_class_df,
          file.path(OUT_DIR, "f1_per_celltype.csv"),
          row.names = FALSE)

write.csv(
  data.frame(
    cell_id          = colnames(sce_test),
    ground_truth     = gt_test,
    pred_singler     = singler_labels,
    pred_scgpt       = scgpt_labels,
#    pred_geneformer  = geneformer_labels,
    singler_score    = round(singler_scores, 4),
    singler_rejected = is.na(singler_pruned)
  ),
  file.path(OUT_DIR, "predictions_test_set.csv"),
  row.names = FALSE
)

saveRDS(sce_test, file.path(OUT_DIR, "sce_test_annotated.rds"))

# ================================================================
#  10. RIEPILOGO FINALE
# ================================================================
cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║              RIEPILOGO BENCHMARK BE1                ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Test set     : %d cellule                           \n",
            ncol(sce_test)))
cat(sprintf("║  Tipi cellulari: %d                                  \n",
            nlevels(sce$gt)))
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-12s Acc=%.3f  F1macro=%.3f  ARI=%.3f\n",
            "SingleR:",
            metrics$SingleR$accuracy,
            metrics$SingleR$f1_macro,
            metrics$SingleR$ari))
cat(sprintf("║  %-12s Acc=%.3f  F1macro=%.3f  ARI=%.3f\n",
            "scGPT:",
            metrics$scGPT$accuracy,
            metrics$scGPT$f1_macro,
            metrics$scGPT$ari))
cat(sprintf("║  %-12s Acc=%.3f  F1macro=%.3f  ARI=%.3f\n",
            "Geneformer:",
            metrics$Geneformer$accuracy,
            metrics$Geneformer$f1_macro,
            metrics$Geneformer$ari))
cat("╚══════════════════════════════════════════════════════╝\n")
cat("\nOutput in ./", OUT_DIR, "/\n")
