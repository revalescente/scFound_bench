datasets <- c("BE1", "cordblood", "sc_mixology")

purrr::walk(datasets, \(ds) {
  quarto::quarto_render(
    input = "~/scFound_bench/r_scripts/cell_QC/QC_main.qmd",
    execute_params = list(dataset = ds),
    output_file = paste0("QC_", ds, ".html")
  )
})
