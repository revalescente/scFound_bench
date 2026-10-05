## Uso ------------------------------------------------------------------------
source("plots/umap_plots.R")

# sce <- readRDS("~/datasets/sc_mixology/sc_mix.RDS")
reducedDimNames(sce)           # -> verifica i nomi esatti
colnames(colData(sce))         # -> scegli la variabile di colore

COLOR_BY <- "Sample"           # es. "cell_line" / "celltype"
OUT      <- "plots/umap"

## palette condivisa: stessa cell line = stesso colore ovunque
pal <- make_palette(colData(sce)[[COLOR_BY]])

## --- 1. pannello comparativo (la figura da slide) ---------------------------
p_panel <- panel_umap(
  sce,
  dimreds  = c("UMAP", "UMAP_scGPT", "UMAP_GeneFormer"),
  titles   = c("PCA (baseline)", "scGPT", "GeneFormer"),
  color_by = COLOR_BY,
  palette  = pal
)
save_slide_png(p_panel, file.path(OUT, "umap_panel_embeddings.png"),
               width = 13.33, height = 5.2)

## --- 2. singole UMAP a tutta slide ------------------------------------------
for (d in c("UMAP", "UMAP_scGPT", "UMAP_GeneFormer")) {
  p <- plot_umap(sce, d, COLOR_BY, palette = pal,
                 title = d, subtitle = paste(ncol(sce), "cellule"),
                 point_size = 0.5, legend = FALSE)
  save_slide_png(p, file.path(OUT, paste0("umap_", d, ".png")),
                 width = 7, height = 6)
}

## --- 3. versione con predizioni zero-shot (dopo compare_embeddings) ---------
# p <- plot_umap(sce, "UMAP_scGPT", "pred_scGPT", palette = pal,
#                title = "scGPT", subtitle = "cell type predetti (zero-shot, k=10)")
# save_slide_png(p, file.path(OUT, "umap_scGPT_pred.png"), width = 7, height = 6)

## --- 4. celle sbagliate evidenziate (utile per discutere gli errori) --------
# sce$correct <- ifelse(is.na(sce$pred_scGPT), NA,
#                       ifelse(as.character(sce$pred_scGPT) == as.character(sce[[COLOR_BY]]),
#                              "Corretto", "Errato"))
# p <- plot_umap(sce, "UMAP_scGPT", "correct",
#                palette = c(Corretto = "grey80", Errato = "#D55E00"),
#                title = "scGPT - errori di annotazione", label = FALSE, legend = TRUE)
# save_slide_png(p, file.path(OUT, "umap_scGPT_errors.png"), width = 7, height = 6)
