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
  library(paletteer)
})

OUT_DIR       <- "predicts/output_benchmark"

sce_test <- readRDS("~/scFound_bench/r_scripts/predicts/output_benchmark/sce_test_annotated.rds")

gt_test           <- as.character(sce_test$gt)
singler_labels    <- as.character(sce_test$pred_singler)
scgpt_labels      <- as.character(sce_test$pred_scgpt)
geneformer_labels <- as.character(sce_test$pred_geneformer)

ref_labels        <- levels(sce_test$gt)

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
  
  cm       <- confusionMatrix(pred, truth)
  by_class <- as.data.frame(cm$byClass)
  ct_names <- gsub("Class: ", "", rownames(by_class))
  
  # Macro-averaged metrics (unweighted mean across classes)
  f1_mac   <- mean(by_class$F1,          na.rm = TRUE)
  prec_mac <- mean(by_class$Precision,   na.rm = TRUE)
  sens_mac <- mean(by_class$Sensitivity, na.rm = TRUE)
  spec_mac <- mean(by_class$Specificity, na.rm = TRUE)
  
  list(
    method       = method,
    f1_macro     = f1_mac,
    precision    = prec_mac,
    sensitivity  = sens_mac,
    specificity  = spec_mac,
    f1_per_class   = setNames(by_class$F1,          ct_names),
    prec_per_class = setNames(by_class$Precision,   ct_names),
    sens_per_class = setNames(by_class$Sensitivity, ct_names),
    spec_per_class = setNames(by_class$Specificity, ct_names),
    cm             = cm$table
  )
}

metrics <- list(
  SingleR    = compute_metrics(singler_labels,    gt_test, "SingleR",    all_labels),
  scGPT      = compute_metrics(scgpt_labels,      gt_test, "scGPT",      all_labels),
  Geneformer = compute_metrics(geneformer_labels, gt_test, "Geneformer", all_labels)
)

summary_df <- do.call(rbind, lapply(metrics, function(m) {
  data.frame(
    Method      = m$method,
    F1_macro    = round(m$f1_macro,    4),
    Precision   = round(m$precision,   4),
    Sensitivity = round(m$sensitivity, 4),
    Specificity = round(m$specificity, 4)
  )
})) %>% arrange(desc(F1_macro))

cat("\n┌─── BENCHMARK SUMMARY ─────────────────────────────────┐\n")
print(summary_df, row.names = FALSE)
cat("└───────────────────────────────────────────────────────┘\n")

# Per-class metrics (long format)
f1_class_df <- do.call(rbind, lapply(metrics, function(m) {
  data.frame(
    Method   = m$method,
    CellType = names(m$f1_per_class),
    F1       = m$f1_per_class
  )
})) %>% filter(!is.na(F1))

# ================================================================
#  8. VISUALIZZAZIONI
# ================================================================
cat("\n=== 7. Visualisations ===\n")
dir.create(OUT_DIR, showWarnings = FALSE)

METHOD_COLS <- c(SingleR    = "#FED789FF",
                 scGPT      = "#023743FF",
                 Geneformer = "#72874EFF")

## 8a. Global metrics bar chart
df_glob <- summary_df %>%
  pivot_longer(cols = c(F1_macro, Precision, Sensitivity, Specificity),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("F1_macro", "Precision",
                                    "Sensitivity", "Specificity")))

p_global <- ggplot(df_glob, aes(x = Metric, y = Value, fill = Method)) +
  geom_col(position = "dodge", alpha = 0.88, width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", Value)),
            position = position_dodge(0.65), vjust = -0.4, size = 4) +
  scale_fill_manual(values = METHOD_COLS) +
  scale_y_continuous(limits = c(0, 1.08), labels = percent_format()) +
  labs(title    = "BE1 cell-type annotation benchmark",
       subtitle = "Test set (30%); ground truth = manual annotation",
       x = NULL, y = NULL, fill = "Method") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        axis.text         = element_text(size = 15))
p_global

ggsave(file.path(OUT_DIR, "1_global_metrics.png"),
       p_global, width = 9, height = 6)

## 8b. F1 per cell type — dot + segment chart
p_f1class <- ggplot(f1_class_df,
                    aes(x = F1, y = reorder(CellType, F1, mean),
                        colour = Method)) +
  geom_line(aes(group = CellType), colour = "grey82", linewidth = 0.8) +
  geom_point(size = 3.5, alpha = 0.92) +
  scale_colour_manual(values = METHOD_COLS) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
  labs(title  = "F1-score by cell type",
       x = "F1-score", y = NULL, colour = "Method") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank())
ggsave(file.path(OUT_DIR, "2_f1_per_celltype.pdf"),
       p_f1class, width = 9, height = max(6, nlevels(sce$gt) * 0.4 + 3))

## 8c. Normalised confusion matrix for each method
plot_cm <- function(cm_table, method_name) {
  mat      <- as.matrix(cm_table)
  mat_norm <- sweep(mat, 2, colSums(mat), "/")  # normalise by true class
  df_cm    <- as.data.frame(as.table(mat_norm))
  colnames(df_cm) <- c("Predicted", "True", "Proportion")
  
  ggplot(df_cm, aes(x = True, y = Predicted, fill = Proportion)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(Proportion > 0.02,
                                 sprintf("%.2f", Proportion), "")),
              size = 2.5, colour = "white", fontface = "bold") +
    scale_fill_gradientn(
      colours = c("#fafafa", "#b8c9e8", "#5b7ec9", "#1a3a6b"),
      limits  = c(0, 1), labels = percent_format(), name = "Proportion"
    ) +
    labs(title = paste("Confusion matrix —", method_name),
         x = "True label", y = "Predicted label") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

for (m in names(metrics)) {
  p <- plot_cm(metrics[[m]]$cm, m)
  ggsave(file.path(OUT_DIR, paste0("3_confusion_", tolower(m), ".png")),
         p, width = max(8, nlevels(sce_test$gt) * 0.6),
         height = max(7, nlevels(sce_test$gt) * 0.55))
}

## 8d. UMAP: ground truth + three methods side by side
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

umap_panel <- (p_gt | p_sr) / (p_sg | p_gf)
ggsave(file.path(OUT_DIR, "4_umap_comparison.pdf"),
       umap_panel, width = 16, height = 14)

## 8e. Heatmap: F1 per cell type × method
f1_mat <- f1_class_df |>
  pivot_wider(names_from = Method, values_from = F1) |>
  textshape::column_to_rownames("CellType") |>
  as.matrix()

ufp_cols <- paletteer_d("rcartocolor::BluYl")

pheatmap(f1_mat,
         cluster_rows = TRUE, cluster_cols = FALSE,
         color = colorRampPalette(ufp_cols)(50),
         breaks   = seq(0, 1, length.out = 51),
         display_numbers = TRUE, number_format = "%.2f",
         number_color = "white", fontsize_number = 15,
         main     = "F1-score by cell type and method",
         filename = file.path(OUT_DIR, "5_heatmap_f1_per_celltype.png"),
         width    = 7, height = max(5, nrow(f1_mat) * 0.4 + 3))



## 8f. Error rate per cell type
error_df <- data.frame(
  CellType   = gt_test,
  SingleR    = as.integer(singler_labels    != gt_test),
  scGPT      = as.integer(scgpt_labels      != gt_test),
  Geneformer = as.integer(geneformer_labels != gt_test)
) %>%
  group_by(CellType) %>%
  summarise(across(c(SingleR, scGPT, Geneformer), mean), .groups = "drop") %>%
  pivot_longer(-CellType, names_to = "Method", values_to = "ErrorRate")

p_err <- ggplot(error_df,
                aes(x = reorder(CellType, ErrorRate, mean),
                    y = ErrorRate, fill = Method)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = METHOD_COLS) +
  scale_y_continuous(labels = percent_format()) +
  coord_flip() +
  labs(title = "Error rate by cell type",
       x = NULL, y = "Error rate", fill = "Method") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank())
ggsave(file.path(OUT_DIR, "6_error_rate_per_celltype.pdf"),
       p_err, width = 10, height = max(6, length(unique(gt_test)) * 0.5 + 3))

cat("  Plots saved to", OUT_DIR, "\n")

# ================================================================
#  9. SAVE RESULTS
# ================================================================
cat("\n=== 9. Saving results ===\n")

write.csv(summary_df,
          file.path(OUT_DIR, "global_metrics.csv"),
          row.names = FALSE)

write.csv(per_class_df,
          file.path(OUT_DIR, "per_celltype_metrics.csv"),
          row.names = FALSE)

write.csv(
  data.frame(
    cell_id          = colnames(sce_test),
    ground_truth     = gt_test,
    pred_singler     = singler_labels,
    pred_scgpt       = scgpt_labels,
    pred_geneformer  = geneformer_labels,
    singler_score    = round(singler_scores, 4),
    singler_rejected = is.na(singler_pruned)
  ),
  file.path(OUT_DIR, "predictions_test_set.csv"),
  row.names = FALSE
)

saveRDS(sce_test, file.path(OUT_DIR, "sce_test_annotated.rds"))

# ================================================================
#  10. FINAL SUMMARY
# ================================================================
cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║                  BE1 BENCHMARK SUMMARY                      ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Test set   : %d cells                                       \n",
            ncol(sce_test)))
cat(sprintf("║  Cell types : %d                                              \n",
            nlevels(sce$gt)))
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-12s F1=%.3f  Prec=%.3f  Sens=%.3f  Spec=%.3f\n",
            "SingleR:",
            metrics$SingleR$f1_macro,
            metrics$SingleR$precision,
            metrics$SingleR$sensitivity,
            metrics$SingleR$specificity))
cat(sprintf("║  %-12s F1=%.3f  Prec=%.3f  Sens=%.3f  Spec=%.3f\n",
            "scGPT:",
            metrics$scGPT$f1_macro,
            metrics$scGPT$precision,
            metrics$scGPT$sensitivity,
            metrics$scGPT$specificity))
cat(sprintf("║  %-12s F1=%.3f  Prec=%.3f  Sens=%.3f  Spec=%.3f\n",
            "Geneformer:",
            metrics$Geneformer$f1_macro,
            metrics$Geneformer$precision,
            metrics$Geneformer$sensitivity,
            metrics$Geneformer$specificity))
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\nOutput in ./", OUT_DIR, "/\n")
