library(purrr)
library(quarto)

datasets <- c("BE1","cordblood", "sc_mixology")
output_dir <- "~/scFound_bench/r_scripts/cell_QC/QC_quarto"

# tutti
purrr::walk(datasets, \(ds) {
  quarto::quarto_render(
    input = file.path(output_dir, "QC_main.qmd"),
    output_file = sprintf("QC_main_%s.html", ds),
    execute_params = list(dataset = ds)
  )
})

# singolo
quarto::quarto_render(
  input = file.path(output_dir, "QC_main.qmd"),
  execute_params = list(dataset = "BE1")
)
