library(purrr)
library(DropletQC)

# 1. Define the paths
base_dir <- "~/cite_seq_coord_blood_raw_data"
gtf_path <- "~/ref/refdata-cellranger-hg19-3.0.0/genes/genes.gtf"
cores    <- 4  # Match the 4 cores available on your machine

# 2. Find all SRR directories in the base folder
srr_dirs <- list.dirs(base_dir, recursive = FALSE)
srr_dirs <- srr_dirs[grepl("SRR", basename(srr_dirs))] # Only keep folders starting with "SRR"

message("Found ", length(srr_dirs), " potential SRR folders.")

# 3. Loop through each folder using purrr::walk
walk(srr_dirs, function(srr_path) {
  
  srr_id <- basename(srr_path)
  
  # Define the paths to the BAM, barcodes, and where we will save the output
  bam_file     <- file.path(srr_path, "STARsolo_output", "Aligned.sortedByCoord.out.bam")
  barcode_file <- file.path(srr_path, "STARsolo_output", "Solo.out", "Gene", "filtered", "barcodes.tsv")
  output_csv   <- file.path(srr_path, paste0(srr_id, "_nuclear_fraction.csv"))
  
  # A. Check if we already processed this one!
  if (file.exists(output_csv)) {
    message("\n[SKIP] Already processed: ", srr_id)
    return(NULL) # Move to the next one
  }
  
  # B. Check if it's actually an RNA folder (meaning the BAM exists)
  if (!file.exists(bam_file) || !file.exists(barcode_file)) {
    message("\n[SKIP] No BAM/barcodes found (likely an ADT folder): ", srr_id)
    return(NULL) # Move to the next one
  }
  
  message("\n==================================================")
  message(">>> Running DropletQC on RNA sample: ", srr_id)
  message("==================================================")
  
  # C. Run DropletQC safely using tryCatch (in case one BAM is corrupted)
  tryCatch({
    nf <- nuclear_fraction_annotation(
      annotation_path = gtf_path, 
      bam = bam_file,
      barcodes = barcode_file,
      tiles = 1, 
      cores = cores, 
      verbose = TRUE
    )
    
    # Save the result to a CSV file inside the SRR folder
    write.csv(nf, file = output_csv, row.names = TRUE)
    message("✅ Successfully saved: ", output_csv)
    
  }, error = function(e) {
    message("❌ ERROR processing ", srr_id, ": ", e$message)
  })
})

message("\n🎉 DropletQC pipeline finished!")