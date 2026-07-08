library(DropletQC)
library(purrr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scuttle)
library(scran)
library(scater)


sc_mix <- readRDS("/projects/shared/intronic_bam/datasets/sc_mixology/sc_mix.RDS")

hist(prop.table(sc_mix$mapped_to_exon))

summary(sc_mix$mapped_to_exon)
head(colData(sc_mix))

# calcoliamo percentuale di sequence esoniche

# Calcoliamo il totale delle reads grezze per ogni cellula
total_sequenced_reads <- colData(sc_mix)$unaligned + 
  colData(sc_mix)$aligned_unmapped + 
  colData(sc_mix)$mapped_to_exon + 
  colData(sc_mix)$mapped_to_intron + 
  colData(sc_mix)$ambiguous_mapping + 
  colData(sc_mix)$mapped_to_MT + 
  colData(sc_mix)$mapped_to_ERCC

# Calcoliamo la percentuale e salviamola nel colData
colData(sc_mix)$nuclear_fraction <- (colData(sc_mix)$mapped_to_exon / total_sequenced_reads)

# Distribution of nuclear fraction
p <-ggplot(as.data.frame(colData(sc_mix)), aes(x=nuclear_fraction))
p + geom_density() +
  facet_wrap(~cell_line) +
  labs(title = "Distribution of nuclear fraction for each cell type - Mixology dataset")+
  xlim(c(0,1))

# Distribution of UMI counts
p + aes(y=log10(sum)) +
  aes(colour=cell_line) +
  geom_point(alpha = 0.7)


# Let's identify empty droplet ----
# Get data frame with the nuclear fraction in the first column and umi counts in
# the second
nf_umi <- data.frame(nf=sc_mix$nuclear_fraction,
                     umi=sc_mix$sum)
head(nf_umi)
# Run identify_empty_drops
cb_mix <- identify_empty_drops(nf_umi=nf_umi)
head(cb_mix)
table(cb_mix$cell_status)

# This function tries to identify the population of empty droplets, 
# but can fail if the population is very small or there are none. 
# To check if the population of empty droplets has been identified correctly 
# it can be instructive to visualise where the cut-off has been drawn:

cb_mix <- identify_empty_drops(nf_umi=nf_umi, include_plot = TRUE)

# Useful parameters: nf_rescue = 0.05 (default), umi_rescue = 1000 (default, log10 = 3)


# add reduced dims
sc_mix <- logNormCounts(sc_mix)
dec <- modelGeneVar(sc_mix)
hvgs <- getTopHVGs(dec, prop = 0.1)

sc_mix <- runPCA(sc_mix, subset_row = hvgs, ncomponents = 50)
sc_mix <- runUMAP(sc_mix, dimred = "PCA")
sc_mix

umapp <- plotUMAP(sc_mix, colour_by = "cell_line")
pcap <- plotPCA(sc_mix, colour_by = "cell_line")

umapb <- plotUMAP(sc_mix, colour_by = "nuclear_fraction")
pcab <- plotPCA(sc_mix, colour_by = "nuclear_fraction")

pcap | pcab
umapp | umapb


# identify damaged cells
cb_mix$cell_line <- sc_mix$cell_line
cb_mix_dc <- identify_damaged_cells(cb_mix, verbose = FALSE, output_plots = TRUE)
table(cb_mix_dc[[1]]$cell_status)

wrap_plots(cb_mix_dc[[2]], nrow = 5)

# Qui i valori sono altissimi... 

# salvo oggetto

saveRDS(sc_mix, "/projects/shared/intronic_bam/datasets/sc_mixology/sc_mix.RDS")
