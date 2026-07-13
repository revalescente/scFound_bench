qmdfiles <- c("~/scFound_bench/r_scripts/cell_QC/BE1/QC_be1.qmd",
  "~/scFound_bench/r_scripts/cell_QC/cordblood/QC_cordblood.qmd",
  "~/scFound_bench/r_scripts/cell_QC/mixology/QC_mixology.qmd")

purrr::map(qmdfiles, \(file) {quarto::quarto_render(input = file)})
