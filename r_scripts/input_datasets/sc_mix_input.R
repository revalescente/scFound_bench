### create sc_mix data ####

# download the data from here:
# https://github.com/LuyiTian/sc_mixology/blob/master/data/sincell_with_class_5cl.RData

# 5 cell line 
load("~/datasets/sc_mixology/sincell_with_class_5cl.RData")

sce_sc_10x_5cl_qc

# create RData object ####

saveRDS(sce_sc_10x_5cl_qc, file = "~/datasets/sc_mixology/sc_mix.RDS")


# save h5ad object #####

# library(zellkonverter)
# out_path <- tempfile(pattern = ".h5ad")

# writeH5AD(sce_sc_10x_5cl_qc, file = "sc_mixolgy_10x_5cl.h5ad")