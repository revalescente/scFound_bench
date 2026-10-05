## =============================================================================
## UMAP "da slide" per embedding di foundation models
## PCA / scGPT / GeneFormer  ->  PNG 300 dpi
##
## Scelte di design:
##  - palette Okabe-Ito (colorblind-safe, 8 colori max)
##  - colore legato all'entita' (stessa label = stesso colore in tutti i pannelli)
##  - niente griglia/assi/tick: solo due frecce d'orientamento in basso a sinistra
##  - label direttamente sui cluster (niente legenda gigante) quando <= 8 gruppi
##  - punti rasterizzati: PNG leggero anche con 26k celle
## =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
## opzionale ma consigliato: install.packages("ggrastr")
.has_rastr <- requireNamespace("ggrastr", quietly = TRUE)

## --- Palette -----------------------------------------------------------------
## Okabe-Ito: 8 tinte con separazione garantita per deuteranopia/protanopia
OKABE_ITO <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
               "#D55E00", "#56B4E9", "#F0E442", "#000000")

#' Palette fissa per le label: stesso colore in ogni pannello e in ogni figura
make_palette <- function(levels, palette = OKABE_ITO) {
  levels <- sort(unique(as.character(levels)))
  if (length(levels) > length(palette)) {
    warning("piu' di ", length(palette), " gruppi: considera un raggruppamento 'Other'")
    palette <- grDevices::colorRampPalette(palette)(length(levels))
  }
  setNames(palette[seq_along(levels)], levels)
}

## --- Tema --------------------------------------------------------------------
theme_umap <- function(base_size = 13, base_family = "") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size * 1.15,
                                     hjust = 0, margin = margin(b = 2)),
      plot.subtitle   = element_text(colour = "grey35", size = base_size * 0.85,
                                     hjust = 0, margin = margin(b = 6)),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = base_size * 0.9),
      legend.text     = element_text(size = base_size * 0.85, colour = "grey20"),
      legend.key      = element_blank(),
      plot.margin     = margin(10, 10, 10, 10),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

## --- Frecce d'orientamento (stile "DimPlot minimal") -------------------------
.axis_arrows <- function(df, lab_x = "UMAP 1", lab_y = "UMAP 2",
                         frac = 0.18, size = 3) {
  rx <- range(df$x); ry <- range(df$y)
  x0 <- rx[1] - diff(rx) * 0.02
  y0 <- ry[1] - diff(ry) * 0.02
  list(
    annotate("segment", x = x0, xend = x0 + diff(rx) * frac, y = y0, yend = y0,
             arrow = arrow(length = unit(2, "mm"), type = "closed"),
             linewidth = 0.5, colour = "grey30"),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + diff(ry) * frac,
             arrow = arrow(length = unit(2, "mm"), type = "closed"),
             linewidth = 0.5, colour = "grey30"),
    annotate("text", x = x0 + diff(rx) * frac / 2, y = y0 - diff(ry) * 0.035,
             label = lab_x, size = size, colour = "grey30", hjust = 0.5, vjust = 1),
    annotate("text", x = x0 - diff(rx) * 0.035, y = y0 + diff(ry) * frac / 2,
             label = lab_y, size = size, colour = "grey30", angle = 90,
             hjust = 0.5, vjust = 0),
    coord_cartesian(clip = "off")
  )
}

## =============================================================================
## Plot di una singola UMAP
## =============================================================================
#' @param sce SingleCellExperiment
#' @param dimred nome in reducedDims(): "UMAP", "UMAP_scGPT", "UMAP_GeneFormer"
#' @param color_by colonna di colData (categoriale) oppure un gene (se is_gene=TRUE)
#' @param palette named vector colore->label (usa make_palette() per condividerla)
#' @param label TRUE = scrive il nome del gruppo sul centroide
#' @param legend TRUE = mostra anche la legenda
#' @param point_size / point_alpha estetica dei punti
#' @param shuffle mescola le righe: evita che l'ultimo gruppo copra gli altri
plot_umap <- function(sce,
                      dimred,
                      color_by,
                      palette      = NULL,
                      title        = dimred,
                      subtitle     = NULL,
                      label        = TRUE,
                      legend       = FALSE,
                      point_size   = 0.35,
                      point_alpha  = 0.75,
                      shuffle      = TRUE,
                      seed         = 42,
                      raster_dpi   = 300) {

  stopifnot(dimred %in% reducedDimNames(sce))
  stopifnot(color_by %in% colnames(colData(sce)))

  xy <- as.matrix(reducedDim(sce, dimred))[, 1:2, drop = FALSE]
  df <- data.frame(x = xy[, 1], y = xy[, 2],
                   grp = as.character(colData(sce)[[color_by]]),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$grp), ]
  if (shuffle) { set.seed(seed); df <- df[sample(nrow(df)), ] }

  if (is.null(palette)) palette <- make_palette(df$grp)
  df$grp <- factor(df$grp, levels = names(palette))

  ## centroidi (mediana = robusta agli outlier) per le label dirette
  cen <- aggregate(cbind(x, y) ~ grp, data = df, FUN = stats::median)

  pt_layer <- if (.has_rastr) {
    ggrastr::rasterise(
      geom_point(size = point_size, alpha = point_alpha, shape = 16, stroke = 0),
      dpi = raster_dpi
    )
  } else {
    geom_point(size = point_size, alpha = point_alpha, shape = 16, stroke = 0)
  }

  p <- ggplot(df, aes(x, y, colour = grp)) +
    pt_layer +
    scale_colour_manual(values = palette, drop = FALSE, name = NULL) +
    .axis_arrows(df) +
    labs(title = title, subtitle = subtitle) +
    theme_umap()

  if (label) {
    p <- p +
      ggrepel::geom_label_repel(
        data = cen, aes(x, y, label = grp),
        inherit.aes = FALSE,
        fill = alpha("white", 0.78), label.size = 0, label.padding = unit(1.4, "mm"),
        colour = "grey10", fontface = "bold", size = 3.6,
        box.padding = 0.3, min.segment.length = 0.25,
        segment.colour = "grey45", segment.size = 0.3, seed = seed
      )
  }

  if (legend) {
    p <- p + guides(colour = guide_legend(
      override.aes = list(size = 3.2, alpha = 1))) +
      theme(legend.position = "right")
  } else {
    p <- p + theme(legend.position = "none")
  }
  p
}

## =============================================================================
## Pannello comparativo: PCA | scGPT | GeneFormer
## =============================================================================
#' @param dimreds vettore di reducedDimNames da affiancare
#' @param titles etichette dei pannelli (default = dimreds)
panel_umap <- function(sce,
                       dimreds  = c("UMAP", "UMAP_scGPT", "UMAP_GeneFormer"),
                       titles   = c("PCA (baseline)", "scGPT", "GeneFormer"),
                       color_by = "Sample",
                       palette  = NULL,
                       nrow     = 1,
                       legend   = TRUE,
                       ...) {

  if (is.null(palette))
    palette <- make_palette(colData(sce)[[color_by]])

  ps <- Map(function(d, t) {
    plot_umap(sce, dimred = d, color_by = color_by, palette = palette,
              title = t, legend = FALSE, label = TRUE, ...)
  }, dimreds, titles)

  out <- patchwork::wrap_plots(ps, nrow = nrow)

  if (legend) {
    ## una sola legenda condivisa a destra
    leg <- plot_umap(sce, dimreds[1], color_by, palette, legend = TRUE)
    out <- out + patchwork::plot_layout(guides = "collect") &
      scale_colour_manual(values = palette, drop = FALSE, name = NULL) &
      guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1))) &
      theme(legend.position = "right")
  }
  out & theme(plot.background = element_rect(fill = "white", colour = NA))
}

## =============================================================================
## Salvataggio PNG pronto per le slide (16:9)
## =============================================================================
#' @param width/height in pollici; 300 dpi -> PNG nitido anche a schermo intero
save_slide_png <- function(plot, file, width = 13.33, height = 5.2, dpi = 300) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(file, plot, width = width, height = height, dpi = dpi,
         units = "in", bg = "white")
  message("salvato: ", normalizePath(file))
  invisible(file)
}
