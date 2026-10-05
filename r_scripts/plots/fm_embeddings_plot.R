library(SingleCellExperiment)
library(scater)
library(patchwork)
library(Nebulosa)


SCE_PATH      <- "/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS"
sce_be1 <- readRDS(SCE_PATH)

## --- parametri ---------------------------------------------------------------
DR      <- c(scGPT = "UMAP_scGPT", GeneFormer = "UMAP_GeneFormer")  # controlla reducedDimNames(sce_be1)
CT      <- "Sample"        # tipo cellulare in BE1 (cambia se hai una colonna dedicata)
STATUS  <- "cell_status"
GENE    <- "MALAT1"

LEGEND_LABELS <- c(
  Sample      = "Cell types",
  cell_status = "Cell Status"
)

## se i rownames sono Ensembl ID, mappa il simbolo:
# GENE <- rownames(sce_be1)[match("MALAT1", rowData(sce_be1)$gene_names)]

## palette colorblind-safe, fissa per label (stesso colore nei due pannelli)
ok <- c("#0072B2","#F69A00","#009E73","#CC79A7","#D55E00","#56B4E9","#F0E442","#000000")
pal_ct     <- setNames(ok[seq_along(levels(factor(sce_be1[[CT]])))],     levels(factor(sce_be1[[CT]])))
pal_status <- setNames(ok[seq_along(levels(factor(sce_be1[[STATUS]])))], levels(factor(sce_be1[[STATUS]])))

theme_umap <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size = 9, colour = "grey40"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.key = element_blank(), plot.background = element_rect(fill = "white", colour = NA))

## --- helper: subset temporaneo, solo per il pannello corrente ----------------
.drop_na_cells <- function(sce, dimred, colour_by, exprs_values = NULL) {
  keep <- stats::complete.cases(reducedDim(sce, dimred)[, 1:2, drop = FALSE])
  
  if (colour_by %in% colnames(colData(sce))) {
    keep <- keep & !is.na(sce[[colour_by]])
  } else if (colour_by %in% rownames(sce)) {
    keep <- keep & !is.na(assay(sce, exprs_values %||% "logcounts")[colour_by, ])
  } else {
    stop("'", colour_by, "' non e' ne' una colonna di colData ne' un gene")
  }
  
  n_na <- sum(!keep)
  if (n_na > 0) message("[", dimred, " | ", colour_by, "] rimosse ", n_na,
                        " celle con NA (", round(100 * n_na / ncol(sce), 2), "%)")
  sce[, keep]
}
`%||%` <- function(a, b) if (is.null(a)) b else a



## --- helper di plot -----------------------------------------------------------
umap_cat <- function(sce, dimred, colour_by, palette, title = NULL,
                     legend_title = NULL) {
  s <- .drop_na_cells(sce, dimred, colour_by)
  s[[colour_by]] <- droplevels(factor(s[[colour_by]]))   # niente categorie vuote in legenda
  
  ## etichetta legenda: esplicita > mappa > nome della colonna
  if (is.null(legend_title))
    legend_title <- if (colour_by %in% names(LEGEND_LABELS)) LEGEND_LABELS[[colour_by]] else colour_by
  
  plotReducedDim(s, dimred = dimred, colour_by = colour_by,
                 point_size = 0.2, point_alpha = 1, rasterise = FALSE) +
    scale_colour_manual(values = palette, name = legend_title, drop = FALSE,
                        aesthetics = c("colour", "fill")) +
    guides(colour = guide_legend(title = legend_title,
                                 override.aes = list(size = 4, alpha = 1))) +
    ggtitle(title) + theme_umap
}

umap_gene <- function(sce, dimred, gene, title = NULL) {
  s <- .drop_na_cells(sce, dimred, gene, exprs_values = "logcounts")
  
  plotReducedDim(s, dimred = dimred, colour_by = gene,
                 by_exprs_values = "logcounts", order_by = gene,
                 point_size = 0.2, point_alpha = 1, rasterise = FALSE) +
    scale_colour_viridis_c(name = paste0(gene, "\n(logcounts)"),
                           aesthetics = c("colour", "fill")) +
    ggtitle(title) + theme_umap
}


## --- pannelli: colonna 1 = scGPT, colonna 2 = GeneFormer ---------------------
p1 <- umap_cat (sce_be1, DR["scGPT"],      CT,     pal_ct)
p2 <- umap_cat (sce_be1, DR["GeneFormer"], CT,     pal_ct)
(p3 <- umap_cat (sce_be1, DR["scGPT"],      STATUS, pal_status))
(p4 <- umap_cat (sce_be1, DR["GeneFormer"], STATUS, pal_status))
(p5 <- umap_gene(sce_be1, DR["scGPT"],      GENE))
(p6 <- umap_gene(sce_be1, DR["GeneFormer"], GENE))




# separate per FM
fig2a <- (p1 | p3 | p5) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "BE1 - scGPT embeddings UMAP",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

fig2a

fig2b <- (p2 | p4 | p6) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "BE1 - Geneformer embeddings UMAP",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

fig2b

ggsave("plots/umap/be1_scGPT_embs.png", fig2a,
       width = 22, height = 5, dpi = 300, bg = "white")
ggsave("plots/umap/be1_geneformer_embs.png", fig2b,
       width = 22, height = 5, dpi = 300, bg = "white")


# ------------- let's try Nebulosa

plot_density(sce_be1, "MALAT1", reduction = "UMAP_scGPT")
plot_density(sce_be1, "MALAT1", reduction = "UMAP_GeneFormer", size = 1) +
  theme(plot.title = element_blank())


# sistemazione umap malat1

umap_gene <- function(sce, dimred, gene, title = NULL) {
  s <- .drop_na_cells(sce, dimred, gene, exprs_values = "logcounts")
  
  ## i punti scuri (bassi) vengono disegnati per ultimi -> stanno sopra
  s$.plot_order <- -as.numeric(logcounts(s)[gene, ])
  
  plotReducedDim(s, dimred = dimred, colour_by = gene,
                 by_exprs_values = "logcounts", order_by = ".plot_order",
                 point_size = 0.6, point_alpha = 0.7, rasterise = FALSE) +
    scale_colour_viridis_c(name = paste0(gene, "\n(logcounts)"),
                           aesthetics = c("colour", "fill")) +
    ggtitle(title) + theme_umap
}
(p5 <- umap_gene(sce_be1, DR["scGPT"],      GENE))
(p6 <- umap_gene(sce_be1, DR["GeneFormer"], GENE))


# bene!!

