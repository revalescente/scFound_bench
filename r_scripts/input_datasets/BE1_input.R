
library(DropletUtils)
library(SingleCellExperiment)
library(purrr)

# Data creation
data_dir <- "/home/revalescente/datasets/BE1/BE1run12"
(dataset_list <- list.files(data_dir))

sce_list <- map(dataset_list, \(dataset) {
  current_data_dir <- file.path(data_dir, dataset)
  
  read10xCounts(
    samples = current_data_dir,
    sample.names = dataset,
    col.names = TRUE)
})
names(sce_list) <- dataset_list
sce <- do.call(cbind, sce_list)

# Saving the object ####

saveRDS(sce, file = "/home/revalescente/datasets/BE1/BE1_sce.RDS")
