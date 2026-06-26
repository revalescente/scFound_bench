library(DropletQC)
library(purrr)

base_dir  <- "/home/vreffo/scFound_bench/results"
samples <- list.files(base_dir, full.names = FALSE, recursive = FALSE)

list_nf <- map(samples, \(sample){
  nf <- readRDS(paste0("~/scFound_bench/results/",sample))
  nf
})
names(list_nf) <- stringr::str_extract(samples, "SRR\\d+")
str(list_nf)

head(list_nf[[1]])

sce_cb <- readRDS("/projects/shared/intronic_bam/datasets/cord_blood/cord_blood_sce.RDS")
sce_cb
sapply(list_nf, dim)
