suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scrapper)
})

SCE_PATH      <- "/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS"

sce <- readRDS(SCE_PATH)
sce

table(sce$cell_status)
table(sce$discard)

summary(sce$subsets_Mito_percent)

is_mito <- grepl("^Mt-", rownames(sce), ignore.case = TRUE)
table(is_mito)

sce <- quickRnaQc.se(
  sce, 
  subsets=list(MT=is_mito)
)

table(sce$keep, sce$cell_status, useNA = "always")

sce$low_qc <- !sce$keep | sce$cell_status == "empty_droplet"

