library(purrr)
library(quarto)

datasets <- c("BE1","cordblood", "sc_mixology")
output_dir <- "~/scFound_bench/r_scripts/cell_QC/QC_quarto"

# tutti
purrr::walk(datasets, \(ds) {
  # Creiamo un file qmd specifico per il dataset per separare la cache
  qmd_name <- sprintf("QC_main_%s.qmd", ds)
  qmd_path <- file.path(output_dir, qmd_name)
  file.copy(file.path(output_dir, "QC_main.qmd"), qmd_path, overwrite = TRUE)
  
  quarto::quarto_render(
    input = qmd_path,
    execute_params = list(dataset = ds)
  )
})

# singolo
# Esempio per singolo con cache separata
ds_single <- "BE1"
qmd_single <- file.path(output_dir, sprintf("QC_main_%s.qmd", ds_single))
file.copy(file.path(output_dir, "QC_main.qmd"), qmd_single, overwrite = TRUE)
quarto::quarto_render(
  input = qmd_single,
  execute_params = list(dataset = ds_single)
)
