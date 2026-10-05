## =============================================================================
## Zero-shot cell type annotation from foundation-model embeddings
##
## Porting in R della procedura usata nei notebook scGPT
## (py_scripts/scGPT/zero_shot_annotation/zero_shot_*.ipynb):
##   1. split reference (train) / query (test)
##   2. kNN sugli embedding con distanza L2  -> equivalente a faiss.IndexFlatL2
##   3. majority voting sui kappa vicini
##   4. metriche macro-averaged (accuracy, precision, recall, F1)
##   5. confusion matrix normalizzata per riga (true label)
##
## Funziona con qualunque embedding salvato in reducedDims(sce):
## "scGPT", "Geneformer", "PCA", ...
##
## Autore: Valerio Reffo
## =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(BiocNeighbors)   # kNN esatto/approssimato, multicore
  library(BiocParallel)
  library(pheatmap)
})

## -----------------------------------------------------------------------------
## Helper: majority voting sui vicini
## -----------------------------------------------------------------------------
#' @param nn_idx matrice n_query x kappa di indici nel reference
#' @param ref_labels vettore character con le label del reference
#' @return data.frame con pred, n_votes, confidence
.majority_vote <- function(nn_idx, ref_labels) {
  lvls <- sort(unique(ref_labels))
  # matrice di label dei vicini
  nn_lab <- matrix(ref_labels[nn_idx], nrow = nrow(nn_idx), ncol = ncol(nn_idx))

  # conteggio voti per classe (veloce: tabulate per riga su codici interi)
  codes <- matrix(match(nn_lab, lvls), nrow = nrow(nn_lab))
  counts <- t(apply(codes, 1, function(r) tabulate(r, nbins = length(lvls))))
  if (length(lvls) == 1L) counts <- matrix(counts, ncol = 1L)

  win <- max.col(counts, ties.method = "first")  # come pandas value_counts()[0]
  n_votes <- counts[cbind(seq_len(nrow(counts)), win)]

  data.frame(
    pred       = lvls[win],
    n_votes    = n_votes,
    confidence = n_votes / ncol(nn_idx),
    stringsAsFactors = FALSE
  )
}

## -----------------------------------------------------------------------------
## Helper: metriche macro-averaged (come sklearn con average="macro")
## -----------------------------------------------------------------------------
.classification_metrics <- function(actual, preds,
                                    labels = NULL,
                                    zero_division = 0) {
  actual <- as.character(actual)
  preds  <- as.character(preds)
  # sklearn di default usa l'unione delle classi presenti in y_true e y_pred
  if (is.null(labels)) labels <- sort(unique(c(actual, preds)))

  per_class <- do.call(rbind, lapply(labels, function(cl) {
    tp <- sum(preds == cl & actual == cl)
    fp <- sum(preds == cl & actual != cl)
    fn <- sum(preds != cl & actual == cl)
    prec <- if ((tp + fp) == 0) zero_division else tp / (tp + fp)
    rec  <- if ((tp + fn) == 0) zero_division else tp / (tp + fn)
    f1   <- if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)
    data.frame(label = cl, support = sum(actual == cl),
               precision = prec, recall = rec, f1 = f1,
               stringsAsFactors = FALSE)
  }))

  list(
    summary = data.frame(
      accuracy  = mean(preds == actual),
      precision = mean(per_class$precision),
      recall    = mean(per_class$recall),
      macro_f1  = mean(per_class$f1),
      weighted_f1 = stats::weighted.mean(per_class$f1, per_class$support)
    ),
    per_class = per_class
  )
}

## -----------------------------------------------------------------------------
## Helper: confusion matrix normalizzata per riga
## -----------------------------------------------------------------------------
.confusion_matrix <- function(actual, preds, normalize = TRUE) {
  actual <- as.character(actual); preds <- as.character(preds)
  lv <- sort(unique(c(actual, preds)))
  cm <- table(Actual = factor(actual, levels = lv),
              Predicted = factor(preds, levels = lv))
  cm <- as.matrix(cm)
  if (normalize) {
    cm <- cm / rowSums(cm)
    cm[is.nan(cm)] <- 0
  }
  cm
}

## =============================================================================
## FUNZIONE PRINCIPALE
## =============================================================================
#' Zero-shot annotation da embeddings
#'
#' @param sce SingleCellExperiment che contiene reference e query
#' @param embed_name nome in reducedDims(sce): "scGPT", "Geneformer", "PCA", ...
#' @param cell_type_key colonna di colData con le label ("Sample", "celltype", "cell_line")
#' @param train_col colonna logica di colData: TRUE = reference, FALSE = query
#' @param kappa numero di vicini (default 10, come nei notebook)
#' @param metric "euclidean" (= faiss IndexFlatL2) oppure "cosine"
#' @param BPPARAM parallelizzazione per la ricerca kNN
#' @param plot se TRUE disegna la confusion matrix
#' @return list(predictions, confidence, metrics, per_class, confusion_matrix, ...)
zero_shot_annotate <- function(sce,
                               embed_name,
                               cell_type_key = "celltype",
                               train_col     = "is_train",
                               kappa         = 10,
                               metric        = c("euclidean", "cosine"),
                               BPPARAM       = BiocParallel::SerialParam(),
                               plot          = TRUE,
                               verbose       = TRUE) {

  metric <- match.arg(metric)
  stopifnot(embed_name %in% reducedDimNames(sce))
  stopifnot(cell_type_key %in% colnames(colData(sce)))
  stopifnot(train_col %in% colnames(colData(sce)))

  if (verbose) message("\n--- Zero-shot annotation | embedding: ", embed_name,
                       " | k = ", kappa, " | metric = ", metric, " ---")

  ## 1. split reference / query -----------------------------------------------
  is_ref   <- as.logical(colData(sce)[[train_col]])
  is_query <- !is_ref
  if (anyNA(is_ref)) stop("NA in colData(sce)$", train_col)

  emb <- as.matrix(reducedDim(sce, embed_name))
  if (anyNA(emb)) stop("NA negli embedding '", embed_name, "'")

  ## la cosine similarity si ottiene come L2 su vettori normalizzati
  if (metric == "cosine") {
    nrm <- sqrt(rowSums(emb^2)); nrm[nrm == 0] <- 1
    emb <- emb / nrm
  }

  ref_embed   <- emb[is_ref, , drop = FALSE]
  query_embed <- emb[is_query, , drop = FALSE]

  ref_labels   <- as.character(colData(sce)[[cell_type_key]][is_ref])
  query_labels <- as.character(colData(sce)[[cell_type_key]][is_query])

  if (verbose) message("Reference: ", nrow(ref_embed), " celle | Query: ",
                       nrow(query_embed), " celle | dim = ", ncol(emb))

  ## 2. kNN (L2 esatto = faiss.IndexFlatL2) -----------------------------------
  knn <- BiocNeighbors::queryKNN(
    X       = ref_embed,
    query   = query_embed,
    k       = min(kappa, nrow(ref_embed)),
    BNPARAM = BiocNeighbors::KmknnParam(),   # esatto
    BPPARAM = BPPARAM
  )

  ## 3. majority voting --------------------------------------------------------
  vote <- .majority_vote(knn$index, ref_labels)
  preds <- vote$pred
  names(preds) <- colnames(sce)[is_query]

  ## 4. metriche ---------------------------------------------------------------
  mt <- .classification_metrics(query_labels, preds)
  if (verbose) {
    print(round(mt$summary, 4))
  }

  ## 5. confusion matrix -------------------------------------------------------
  cm <- .confusion_matrix(query_labels, preds, normalize = TRUE)

  if (plot) {
    cm_plot <- cm
    rownames(cm_plot) <- paste0(rownames(cm_plot), " (True)")
    colnames(cm_plot) <- paste0(colnames(cm_plot), " (Pred)")
    pheatmap::pheatmap(
      cm_plot,
      cluster_rows    = FALSE,
      cluster_cols    = FALSE,
      display_numbers = TRUE,
      number_format   = "%.2f",
      color  = colorRampPalette(c("white", "#6A51A3"))(100),  # cmap "Purples"
      breaks = seq(0, 1, length.out = 101),
      main   = paste0("Confusion matrix - ", embed_name,
                      "  (righe = true, colonne = pred)"),
      angle_col = 45
    )
  }

  list(
    embedding        = embed_name,
    predictions      = preds,
    truth            = setNames(query_labels, colnames(sce)[is_query]),
    confidence       = setNames(vote$confidence, colnames(sce)[is_query]),
    metrics          = mt$summary,
    per_class        = mt$per_class,
    confusion_matrix = cm,
    knn_distances    = knn$distance
  )
}

## =============================================================================
## WRAPPER: confronto fra più embedding (scGPT vs Geneformer vs PCA ...)
## =============================================================================
#' @param embeddings vettore di nomi di reducedDims da confrontare
#' @return list(results = list per embedding, metrics = data.frame comparativo)
compare_embeddings <- function(sce,
                               embeddings    = reducedDimNames(sce),
                               cell_type_key = "celltype",
                               train_col     = "is_train",
                               kappa         = 10,
                               metric        = "euclidean",
                               store_in_sce  = TRUE,
                               plot          = TRUE,
                               BPPARAM       = BiocParallel::SerialParam()) {

  res <- lapply(embeddings, function(e) {
    zero_shot_annotate(sce, embed_name = e, cell_type_key = cell_type_key,
                       train_col = train_col, kappa = kappa, metric = metric,
                       plot = plot, BPPARAM = BPPARAM)
  })
  names(res) <- embeddings

  metrics_df <- do.call(rbind, lapply(embeddings, function(e) {
    cbind(embedding = e, kappa = kappa, res[[e]]$metrics)
  }))
  rownames(metrics_df) <- NULL

  ## opzionale: scrivi le predizioni dentro colData(sce)
  if (store_in_sce) {
    for (e in embeddings) {
      col_pred <- paste0("pred_", e)
      col_conf <- paste0("conf_", e)
      v <- rep(NA_character_, ncol(sce)); names(v) <- colnames(sce)
      v[names(res[[e]]$predictions)] <- res[[e]]$predictions
      cf <- rep(NA_real_, ncol(sce)); names(cf) <- colnames(sce)
      cf[names(res[[e]]$confidence)] <- res[[e]]$confidence
      colData(sce)[[col_pred]] <- factor(v)
      colData(sce)[[col_conf]] <- cf
    }
  }

  list(results = res, metrics = metrics_df, sce = sce)
}

## =============================================================================
## SENSIBILITA' A kappa (utile per il benchmark)
## =============================================================================
sweep_kappa <- function(sce, embeddings, cell_type_key, train_col = "is_train",
                        kappas = c(5, 10, 15, 25, 50), metric = "euclidean") {
  do.call(rbind, lapply(kappas, function(k) {
    do.call(rbind, lapply(embeddings, function(e) {
      r <- zero_shot_annotate(sce, e, cell_type_key, train_col, kappa = k,
                              metric = metric, plot = FALSE, verbose = FALSE)
      cbind(embedding = e, kappa = k, r$metrics)
    }))
  }))
}
