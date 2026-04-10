library(MultiAssayExperiment)
library(SingleCellMultiModal)
library(SingleCellExperiment)
library(ExperimentHub)

sce <- CITEseq(DataType="cord_blood", modes="*", dry.run=FALSE, version="1.0.0",
               DataClass="SingleCellExperiment")

gene_all <- rownames(sce)
gene_m <- gene_all[grep("HUMAN_", gene_all)]

sce_m <- sce[gene_m, ]
rownames(sce_m) <- sub("^HUMAN_", "", rownames(sce_m))

class(counts(sce_m))

counts(sce_m) <- as.matrix(counts(sce_m))
sce <- sce_m
sce
# create RData object ####
saveRDS(sce, file = "~/datasets/cord_blood/cord_blood_sce.RDS")
