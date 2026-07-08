library(DropletQC)
library(purrr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scuttle)
library(scran)
library(scater)

samples <- c("A549", "CCL-185-IG", "PC9", "CRL5868", "HTB178", "DV90", "HCC78", "PBMCs")

# read nuclear fraction list
nf_list <- imap(samples, \(sample, name) {
  s_path <- paste0("~/scFound_bench/results/nf_BE1_", sample,".rds")
  nf <- readRDS(s_path)
  nf
})
names(nf_list) <- samples

# read sce object of be1
sce_be1 <- readRDS("/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS")
sce_be1

# abbina nf al sce
# 1. Combiniamo la lista in un unico data frame
nf_combined <- do.call(rbind, lapply(names(nf_list), function(s) {
  df <- as.data.frame(nf_list[[s]])
  # Crea la chiave identica ai colnames di sce_be1 (es. "AAACCCAAGCCTATCA-1_A549")
  df$cell_key <- paste0(rownames(df), "_", s)
  return(df)
}))

# 2. Mappiamo i valori nel colData usando il match con i colnames di sce_be1
# Questo approccio è sicuro al 100% anche se mancano cellule o l'ordine è diverso
colData(sce_be1)$nuclear_fraction <- nf_combined$nuclear_fraction[match(colnames(sce_be1), nf_combined$cell_key)]

# 3. (Opzionale) Verifica che il caricamento sia andato a buon fine
head(colData(sce_be1)[, c("Sample", "Barcode", "nuclear_fraction")])

table(sce_be1$Sample)

# Distribution of nuclear fraction per cell type
p <- ggplot(as.data.frame(colData(sce_be1)), aes(x=nuclear_fraction))
p + geom_density() +
  facet_wrap(~Sample) +
  labs(title = "Distribution of nuclear fraction - BE1 dataset by cell identity")

# Distribution of UMI counts per 
p + aes(y=log10(sum)) +
  aes(colour=Sample) +
  geom_point(alpha = 0.7)

# Let's identify empty droplet
# Get data frame with the nuclear fraction in the first column and umi counts in
# the second
nf_umi <- data.frame(nf=sce_be1$nuclear_fraction,
                     umi=sce_be1$sum)
table(is.na(nf_umi$nf))
nf_umi$nf[is.na(nf_umi$nf)] <- 0

# Run identify_empty_drops
cb_be1 <- identify_empty_drops(nf_umi=nf_umi)
head(cb_be1)
table(cb_be1$cell_status)

# This function tries to identify the population of empty droplets, 
# but can fail if the population is very small or there are none. 
# To check if the population of empty droplets has been identified correctly 
# it can be instructive to visualise where the cut-off has been drawn:

cb_be1 <- identify_empty_drops(nf_umi=nf_umi, include_plot = TRUE)

# Useful parameters: nf_rescue = 0.05 (default), umi_rescue = 1000 (default, log10 = 3)


# add reduced dims
sce_be1 <- logNormCounts(sce_be1)
dec <- modelGeneVar(sce_be1)
hvgs <- getTopHVGs(dec, prop = 0.1)

sce_be1 <- runPCA(sce_be1, subset_row = hvgs, ncomponents = 50)
sce_be1 <- runUMAP(sce_be1, dimred = "PCA")
sce_be1

umapp <- plotUMAP(sce_be1, colour_by = "Sample")
pcap <- plotPCA(sce_be1, colour_by = "Sample")

umapb <- plotUMAP(sce_be1, colour_by = "nuclear_fraction")
pcab <- plotPCA(sce_be1, colour_by = "nuclear_fraction")

pcap | pcab
umapp | umapb

# check MALAT1 and nuclear fraction relationship ----
umapm <- plotUMAP(sce_be1, colour_by = "MALAT1")
pcam <- plotPCA(sce_be1, colour_by = "MALAT1")

pcab | pcam
umapb | umapm

# scatterplot 
plotExpression(sce_be1, 
               features = "MALAT1", 
               x = "nuclear_fraction", 
               other_fields = "Sample",
               exprs_values = "logcounts") +
  facet_wrap(~Sample) +                   # <--- Divide il grafico in pannelli per ogni Sample
  aes(color = Sample) +                   # (Opzionale) Colora i punti in base al Sample
  theme_minimal() +
  labs(title = "MALAT1 vs Nuclear Fraction per Sample",
       y = "Espressione di MALAT1 (logcounts)")

# identify damaged cells
cb_ed$celltype <- sce_be1$celltype
cb_ed_dc <- identify_damaged_cells(cb_ed, verbose = FALSE, output_plots = TRUE)
table(cb_ed_dc[[1]]$cell_status)

wrap_plots(cb_ed_dc[[2]], nrow = 8)

saveRDS(sce_be1, file = "/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS")
