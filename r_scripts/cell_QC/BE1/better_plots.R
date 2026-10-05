suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(anndataR)
})

SCE_PATH      <- "/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS"

sce <- readRDS(SCE_PATH)

sce

