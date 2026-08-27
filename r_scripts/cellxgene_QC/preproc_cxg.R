library(SingleCellExperiment)
library(anndataR)
library(scrapper)
library(scater)
library(DropletUtils)


output_dir <- "/projects/shared/cellxgene/heart/split_datasets"
files <- list.files(output_dir, pattern = "\\.rds$", full.names = TRUE)
files

sce <- readRDS("/projects/shared/cellxgene/heart/split_datasets/1c739a3e-c3f5-49d5-98e0-73975e751201.rds")
sce

table(sce$dataset_id)
table(sce$cell_type)

assay(sce, "counts") <- assay(sce, "X")
assay(sce, "X") <- NULL
sce

is.mito <- grep("^MT-", rowData(sce)$feature_name, ignore.case = TRUE)
is.ribo <- grep("^RP[SL]", rowData(sce)$feature_name, ignore.case = TRUE)

# 3. Visualizza i primi geni trovati per verifica
head(rowData(sce)$feature_name[is.ribo])
head(rowData(sce)$feature_name[is.mito])

sce.qc <- quickRnaQc.se(
  sce, 
  subsets=list(MT=is.mito)
  )
sce.qc
colData(sce.qc)
table(sce.qc$keep)

# rimuoviamo solo cellule con alto contenuto di mitocondri
(qc.thresh <- metadata(sce.qc)$qc$thresholds)
table(sce.qc$subset.proportion.MT > qc.thresh$subset.proportion["MT"])
table(sce.qc$subset.proportion.RB > qc.thresh$subset.proportion["RB"])

# grafici esplorativi
gridExtra::grid.arrange(
  plotColData(sce.qc, y="sum", colour_by="keep", point_size = 0.3) +
    geom_hline(yintercept=qc.thresh$sum, linetype="dashed", color="red") +
    scale_y_log10() +
    ggtitle("Total count"),
  plotColData(sce.qc, y="detected", colour_by="keep", point_size = 0.3) +
    geom_hline(yintercept=qc.thresh$detected, linetype="dashed", color="red") +
    scale_y_log10() +
    ggtitle("Detected features"),
  plotColData(sce.qc, y="subset.proportion.MT", colour_by="keep", point_size = 0.3) + 
    geom_hline(yintercept=qc.thresh$subset.proportion["MT"], linetype="dashed", color="red") +
    ggtitle("Mito prop"),
  ncol=3
)

gridExtra::grid.arrange(
  # 1. Violin Plot - Total count
  plotColData(sce.qc, x = "keep", y = "sum", colour_by = "keep", point_size = 0) +
    scale_y_log10() +
    ggtitle("Total count"),
  
  # 2. Violin Plot - Detected features
  plotColData(sce.qc, x = "keep", y = "detected", colour_by = "keep", point_size = 0) +
    scale_y_log10() +
    ggtitle("Detected features"),
  
  ncol = 2
)

plotColData(sce.qc, x="sum", y="subset.proportion.MT", colour_by="keep") +
  geom_hline(yintercept=qc.thresh$subset.proportion["MT"], linetype="dashed", color="red")

plotColData(sce.qc, x="sum", y="subset.proportion.RB", colour_by="keep") +
  geom_hline(yintercept=qc.thresh$subset.proportion["RB"], linetype="dashed", color="red")

# c'è bias verso un tipo cellulare? 
table(sce.qc$cell_type, sce.qc$keep)

sce.qc$high.ribo <- (sce.qc$subset.proportion.RB > qc.thresh$subset.proportion["RB"])

table(sce.qc$cell_type, sce.qc$high.ribo)

# gli UMI counts come si distribuiscono? gli emptyDrops voglio al massimo 100 UMI
summary(sce.qc$sum)

# verify empty droplet
set.seed(100)
e.out <- emptyDrops(counts(sce), lower = qc.thresh$sum)

# See ?emptyDrops for an explanation of why there are NA values.
summary(e.out$FDR <= 0.001)


# MALAT1
idx <- which(rowData(sce.qc)$feature_name == "MALAT1")

if (length(idx) > 0) {
  sce.qc$malat1_sum <- counts(sce)[idx, ]
  cat("Gene trovato alla riga:", idx, "\n")
} else {
  message("MALAT1 non trovato nei rownames. Verifica se stai usando gli ID Ensembl.")
}

summary(sce.qc$malat1_sum)

plotColData(sce.qc, x = "keep", y = "malat1_sum", colour_by = "keep", point_size = 0) +
  scale_y_log10() +
  ggtitle("MALAT1 counts")


plotColData(sce.qc, x = "sum", y = "malat1_sum", colour_by = "keep", point_size = 0.3) +
  ggtitle("MALAT1 counts")
