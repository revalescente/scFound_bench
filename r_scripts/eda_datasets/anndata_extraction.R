library(SingleCellExperiment)
library(anndataR)
library(purrr)


# To save the raw objects in h5ad

sce_be1 <- readRDS("~/datasets/BE1/BE1_sce.RDS")
sce_mix <- readRDS("~/datasets/sc_mixology/sc_mix.RDS")
sce_cb <- readRDS("~/datasets/cord_blood/cb_sce.RDS")


iwalk(list("be1" = sce_be1, "sc_mix" = sce_mix, "cord_blood" = sce_cb), \(sce, name) {
  print(name)
  if ("ID" %in% colnames(rowData(sce))) {
    rownames(sce) <- rowData(sce)$ID
  } else {
    rownames(sce) <- make.unique(rownames(sce))
  }
  write_h5ad(sce, path = paste0("~/datasets/anndata/", name, ".h5ad"), compression = "gzip")
}) 

# To save train and test dataset splitted

convert_train_h5ad <- function(
    input_dirs = c(
      "~/datasets/BE1",
      "~/datasets/sc_mixology",
      "~/datasets/cord_blood"),
    out_dir = "~/datasets/anndata/"
    ) {

  # Create the output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Loop through each folder
  for (folder in input_dirs) {
    if (!dir.exists(folder)) {
      warning(paste("Directory does not exist, skipping:", folder))
      next
    }
    
    # List all .rds files in the folder
    all_files <- list.files(path = folder, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
    
    # Filter for files that have 'train' or 'test' in their names
    target_files <- all_files[grepl("train|test", basename(all_files), ignore.case = TRUE)]
    
    # Process each matching file
    for (file_path in target_files) {
      # Extract the file name without the extension
      name <- tools::file_path_sans_ext(basename(file_path))
      
      message(paste("Processing:", name, "from", folder))
      
      # Read the RDS file
      sce <- readRDS(file_path)
      
      # Apply your rowname logic
      if ("ID" %in% colnames(rowData(sce))) {
        rownames(sce) <- rowData(sce)$ID
      } else {
        rownames(sce) <- make.unique(rownames(sce))
      }
      
      # Save the object as .h5ad
      out_path <- paste0(out_dir, name, ".h5ad")
      write_h5ad(sce, path = out_path, compression = "gzip")
      
      message(paste("Successfully saved to:", out_path))
    }
  }
}

convert_train_h5ad()
