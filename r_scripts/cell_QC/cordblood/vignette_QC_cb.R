library(DropletQC)
library(purrr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scuttle)
library(scran)
library(scater)


nf <- readRDS("~/scFound_bench/results/nf_cordblood.rds")
sce_cb <- readRDS("/projects/shared/intronic_bam/datasets/cord_blood/cord_blood_sce.RDS")
sce_cb

# qc metrics
# Identifica i geni mitocondriali basandoti sui rownames
is_mito <- grepl("^MT-", rownames(sce_cb), ignore.case = TRUE)
summary(is_mito) # Giusto per verificare quanti ne trova

# Calcola e aggiunge automaticamente le colonne al colData
sce_cb <- addPerCellQC(sce_cb, subsets = list(Mito = is_mito))

hist(nf$nuclear_fraction)
dim(nf)



# Estraiamo solo le righe dei barcode veri
nf_filtrato <- nf[colnames(sce_cb), , drop = FALSE]
hist(nf_filtrato$nuclear_fraction)

# aggiungo nuclear fractiono ai colData
sce_cb$nuclear_fraction_QC <- nf_filtrato$nuclear_fraction

table(sce_cb$celltype)

# Distribution of nuclear fraction per cell type
p <- ggplot(as.data.frame(colData(sce_cb)), aes(x=nuclear_fraction_QC))
p + geom_density() +
  facet_wrap(~celltype) +
  labs(title = "Distribution of nuclear fraction for each cell type - Cord Blood dataset")

# Distribution of UMI counts per 
p <- p + aes(y=log10(sum)) +
  aes(colour=celltype) +
  geom_point(alpha = 0.7)
p


# Let's identify empty droplet
# Get data frame with the nuclear fraction in the first column and umi counts in
# the second
nf_umi <- data.frame(nf=sce_cb$nuclear_fraction_QC,
                         umi=sce_cb$sum)

# Run identify_empty_drops
cb_ed <- identify_empty_drops(nf_umi=nf_umi)
head(cb_ed)
table(cb_ed$cell_status)

# This function tries to identify the population of empty droplets, 
# but can fail if the population is very small or there are none. 
# To check if the population of empty droplets has been identified correctly 
# it can be instructive to visualise where the cut-off has been drawn:

cb_ed <- identify_empty_drops(nf_umi=nf_umi, include_plot = TRUE)
# Useful parameters: nf_rescue = 0.05 (default), umi_rescue = 1000 (default, log10 = 3)

table(cb_ed$celltype, cb_ed$cell_status) # tutte natural killers (?)

# add reduced dims
sce_cb <- logNormCounts(sce_cb)
dec <- modelGeneVar(sce_cb)
hvgs <- getTopHVGs(dec, prop = 0.1)

sce_cb <- runPCA(sce_cb, subset_row = hvgs, ncomponents = 50)
sce_cb <- runUMAP(sce_cb, dimred = "PCA")
sce_cb

umapp <- plotUMAP(sce_cb, colour_by = "celltype")
pcap <- plotPCA(sce_cb, colour_by = "celltype")

umapb <- plotUMAP(sce_cb, colour_by = "nuclear_fraction_QC")
pcab <- plotPCA(sce_cb, colour_by = "nuclear_fraction_QC")

pcap | pcab
umapp | umapb

# check MALAT1 and nuclear fraction relationship ----
# Non c'è MALAT1
umapm <- plotUMAP(sce_cb, colour_by = "MALAT1")
pcam <- plotPCA(sce_cb, colour_by = "MALAT1")

pcab | pcam
umapb | umapm

# scatterplot 
plotExpression(sce_cb, 
               features = "MALAT1", 
               x = "nuclear_fraction_QC", 
               other_fields = "celltype",
               exprs_values = "logcounts") +
  facet_wrap(~Sample) +                   # <--- Divide il grafico in pannelli per ogni Sample
  aes(color = Sample) +                   # (Opzionale) Colora i punti in base al Sample
  theme_minimal() +
  labs(title = "MALAT1 vs Nuclear Fraction per cell type",
       y = "Espressione di MALAT1 (logcounts)")


# identify damaged cells
cb_ed$celltype <- sce_cb$celltype
cb_ed_dc <- identify_damaged_cells(cb_ed, verbose = FALSE, output_plots = TRUE)
table(cb_ed_dc[[1]]$celltype, cb_ed_dc[[1]]$cell_status)

cb_ed_dc$plots$`Natural Killers`

saveRDS(sce_cb, "/projects/shared/intronic_bam/datasets/cord_blood/cord_blood_sce.RDS")
