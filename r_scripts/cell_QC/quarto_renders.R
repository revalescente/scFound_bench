library(purrr)
library(quarto)

datasets <- c("BE1", "cordblood", "sc_mixology")
output_dir <- getwd() # or specify a folder like "~/scFound_bench/reports"

purrr::walk(datasets, \(ds) {
  quarto::quarto_render(
    input = "~/scFound_bench/r_scripts/cell_QC/QC_main.qmd",
    execute_params = list(dataset = ds),
    output_file = file.path(output_dir, paste0("QC_", ds, ".html"))
  )
})
