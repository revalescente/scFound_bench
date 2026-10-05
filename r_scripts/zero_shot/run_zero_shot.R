## Esempio d'uso -----------------------------------------------------------
source("zero_shot/zero_shot_annotation.R")

SCE_PATH      <- "/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS"
sce_be1 <- readRDS(SCE_PATH)

reducedDimNames(sce_be1)   # controlla i nomi: "PCA", "scGPT", "Geneformer"...

## Se il flag train/test non e' ancora in colData, ricostruiscilo dallo split
## usato per generare i .h5ad, es.:
# train_ids <- colnames(readRDS("~/datasets/BE1/train_BE1_sce.RDS"))
# sce_be1$is_train <- colnames(sce_be1) %in% train_ids

## --- singolo embedding ---------------------------------------------------
res_gf <- zero_shot_annotate(
  sce_be1,
  embed_name    = "GeneFormer",
  cell_type_key = "Sample",     # "celltype" per cordblood, "cell_line" per mixology
  train_col     = "is_train",
  kappa         = 10,
  metric        = "euclidean"   # = faiss IndexFlatL2 dei notebook
)

res_gf$metrics
head(res_gf$per_class)

## --- confronto scGPT vs Geneformer vs PCA --------------------------------
cmp <- compare_embeddings(
  sce_be1,
  embeddings    = c("scGPT", "GeneFormer", "PCA"),
  cell_type_key = "Sample",
  train_col     = "is_train",
  kappa         = 10,
  BPPARAM       = BiocParallel::MulticoreParam(4)  # SnowParam(4) su Windows
)

cmp$metrics          # tabella comparativa
sce_be1 <- cmp$sce   # ora contiene pred_scGPT / conf_scGPT / pred_Geneformer ...

## --- sensibilita' a k ----------------------------------------------------
sweep <- sweep_kappa(sce_be1, c("scGPT"),
                     cell_type_key = "Sample", kappas = c(5, 10, 25, 50))
sweep
