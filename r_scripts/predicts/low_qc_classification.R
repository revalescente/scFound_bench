## ---------------------------------------------------------------------------
## Classificazione low_qc da embeddings di foundation models (GeneFormer, scGPT)
## 5-fold CV stratificata | SOLO Random Forest (ranger) | no class weights
##
## Salva cv_setup.rds: fold, etichette, embeddings e funzione di valutazione,
## cosi' low_qc_svm.R (da lanciare quando vuoi) usa ESATTAMENTE gli stessi fold.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(matrixStats)
  library(ranger)
  library(pROC)
  library(PRROC)
  library(purrr)     # >= 1.0 per list_rbind()
  library(tidyr)     # expand_grid()
  library(scrapper)
})

SCE_PATH      <- "/projects/shared/intronic_bam/datasets/BE1/BE1_sce.RDS"

sce <- readRDS(SCE_PATH)
# is_mito <- grepl("^Mt-", rownames(sce), ignore.case = TRUE)
# sce <- quickRnaQc.se(
#   sce, 
#   subsets=list(MT=is_mito)
# )
# sce$low_qc <- !sce$keep | sce$cell_status == "empty_droplet"

set.seed(42)

## --- Parametri -------------------------------------------------------------
EMBEDDINGS <- c("GeneFormer", "scGPT")   # nomi in reducedDimNames(sce)
K          <- 5                          # fold
RF_THREADS <- parallel::detectCores()    # ranger parallelizza internamente
RF_TREES   <- 5000
OUTDIR     <- "predicts/results_low_qc"
dir.create(OUTDIR, showWarnings = FALSE)

Sys.setenv(OMP_NUM_THREADS = RF_THREADS)

## --- 1. Preparazione dati --------------------------------------------------
stopifnot(all(EMBEDDINGS %in% reducedDimNames(sce)))

keep_cells <- !is.na(sce$low_qc)
message(sprintf("Cellule totali: %d | con low_qc non-NA: %d",
                ncol(sce), sum(keep_cells)))

y <- factor(ifelse(sce$low_qc[keep_cells], "low", "pass"),
            levels = c("pass", "low"))   # "low" = classe positiva
cells <- colnames(sce)[keep_cells]
print(table(y))

X_list <- lapply(EMBEDDINGS, function(nm) {
  m <- as.matrix(reducedDim(sce, nm))[keep_cells, , drop = FALSE]
  storage.mode(m) <- "double"

  ## ranger ricava le variabili da colnames(x): se sono NULL fallisce con
  ## "No covariates found." anche a matrice piena (caso scGPT).
  if (is.null(colnames(m)))
    colnames(m) <- paste0(nm, "_", seq_len(ncol(m)))

  na_cells <- rowSums(is.na(m)) > 0
  if (any(na_cells))
    warning(sprintf("%s: %d celle con NA nell'embedding", nm, sum(na_cells)))

  ## namespace esplicito: DelayedArray/sparseMatrixStats esportano colVars
  v  <- matrixStats::colVars(m, na.rm = TRUE)
  ok <- !is.na(v) & v > 0
  if (sum(ok) == 0)
    stop(sprintf("%s: nessuna colonna con varianza > 0 (embedding degenere?)", nm))
  message(sprintf("%s: %d/%d colonne tenute", nm, sum(ok), ncol(m)))

  m[, ok, drop = FALSE]
})
names(X_list) <- EMBEDDINGS
print(sapply(X_list, dim))
print(sapply(X_list, function(m) head(colnames(m), 3)))

## --- 2. Fold stratificati --------------------------------------------------
make_folds <- function(y, k) {
  fold <- integer(length(y))
  for (lev in levels(y)) {
    idx <- which(y == lev)
    fold[sample(idx)] <- rep_len(seq_len(k), length(idx))
  }
  fold
}
fold_id <- make_folds(y, K)
print(table(fold_id, y))

## --- 3. Metriche -----------------------------------------------------------
eval_pred <- function(truth, score, thr = 0.5) {
  ok <- !is.na(score) & !is.na(truth)
  if (any(!ok)) warning(sprintf("eval_pred: %d predizioni NA scartate", sum(!ok)))
  truth <- factor(truth[ok], levels = c("pass", "low")); score <- score[ok]

  pred <- factor(ifelse(score >= thr, "low", "pass"), levels = levels(truth))
  cm <- table(pred = pred, truth = truth)
  TP <- cm["low", "low"];  FP <- cm["low", "pass"]
  FN <- cm["pass", "low"]; TN <- cm["pass", "pass"]

  sens <- TP / (TP + FN); spec <- TN / (TN + FP); prec <- TP / (TP + FP)
  mcc  <- (TP * TN - FP * FN) /
    sqrt(prod(c(TP + FP, TP + FN, TN + FP, TN + FN)) + 1e-12)

  auroc <- as.numeric(pROC::auc(pROC::roc(truth, score, levels = c("pass", "low"),
                                          direction = "<", quiet = TRUE)))
  auprc <- PRROC::pr.curve(scores.class0 = score[truth == "low"],
                           scores.class1 = score[truth == "pass"])$auc.integral

  data.frame(accuracy = (TP + TN) / sum(cm),
             balanced_accuracy = (sens + spec) / 2,
             sensitivity = sens, specificity = spec,
             precision = prec, f1 = 2 * prec * sens / (prec + sens),
             mcc = mcc, auroc = auroc, auprc = auprc,
             TP = TP, FP = FP, FN = FN, TN = TN)
}

## --- 4. CV: Random Forest --------------------------------------------------
## Un fit = una chiamata di funzione: Xtr/Xte/rf restano LOCALI e vengono
## liberati a ogni iterazione, l'ambiente globale resta pulito.
## Sequenziale: ranger usa gia' RF_THREADS thread per fit.

fit_rf_fold <- function(embedding, fold) {
  X  <- X_list[[embedding]]
  tr <- fold_id != fold; te <- fold_id == fold

  ## scaling stimato SOLO sul training del fold (per la RF e' ininfluente,
  ## ma tiene il preprocessing identico a quello dell'SVM)
  ctr <- colMeans(X[tr, , drop = FALSE])
  sdv <- matrixStats::colSds(X[tr, , drop = FALSE]); sdv[sdv == 0] <- 1
  Xtr <- scale(X[tr, , drop = FALSE], center = ctr, scale = sdv)
  Xte <- scale(X[te, , drop = FALSE], center = ctr, scale = sdv)
  stopifnot(ncol(Xtr) > 0, ncol(Xtr) == ncol(Xte), !is.null(colnames(Xtr)))

  t0 <- Sys.time()
  rf <- ranger::ranger(x = Xtr, y = y[tr], num.trees = RF_TREES,
                       mtry = floor(sqrt(ncol(Xtr))), probability = TRUE,
                       num.threads = RF_THREADS, importance = "impurity",
                       seed = 100 + fold)
  p  <- predict(rf, Xte, num.threads = RF_THREADS)$predictions[, "low"]
  tt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("RF %-12s fold %d: %.1f min | OOB err %.4f",
                  embedding, fold, tt, rf$prediction.error))

  list(metrics = cbind(embedding = embedding, model = "RF", fold = fold,
                       eval_pred(y[te], p), time_min = tt,
                       stringsAsFactors = FALSE),
       oof = data.frame(cell = cells[te], embedding = embedding, fold = fold,
                        truth = y[te], p_RF = p, stringsAsFactors = FALSE),
       importance = sort(ranger::importance(rf), decreasing = TRUE)[1:50])
}

rf_res <- tidyr::expand_grid(embedding = EMBEDDINGS, fold = seq_len(K)) |>
  purrr::pmap(fit_rf_fold)

metrics_rf <- purrr::list_rbind(purrr::map(rf_res, "metrics"))
oof        <- purrr::list_rbind(purrr::map(rf_res, "oof"))

## --- 5. Aggregazione -------------------------------------------------------
summary_rf <- aggregate(
  cbind(balanced_accuracy, sensitivity, specificity, auroc, auprc, mcc, time_min) ~
    embedding + model,
  data = metrics_rf, FUN = function(x) c(mean = mean(x), sd = sd(x)))
print(summary_rf)

## Metriche sulle predizioni out-of-fold aggregate (una curva per embedding)
oof_perf <- purrr::list_rbind(purrr::map(EMBEDDINGS, function(emb) {
  d <- oof[oof$embedding == emb, ]
  cbind(embedding = emb, model = "RF", eval_pred(d$truth, d$p_RF))
}))
print(oof_perf)

write.csv(metrics_rf, file.path(OUTDIR, "rf_metrics_per_fold.csv"), row.names = FALSE)
write.csv(oof_perf,   file.path(OUTDIR, "rf_metrics_pooled_oof.csv"), row.names = FALSE)
saveRDS(list(metrics = metrics_rf, oof = oof, res = rf_res),
        file.path(OUTDIR, "rf_cv_results.rds"))

## --- 6. Setup condiviso per l'SVM (e per la baseline QC) -------------------
saveRDS(list(EMBEDDINGS = EMBEDDINGS, K = K, OUTDIR = OUTDIR,
             keep_cells = keep_cells, cells = cells,
             y = y, fold_id = fold_id, X_list = X_list,
             eval_pred = eval_pred, oof_rf = oof, perf_rf = oof_perf),
        file.path(OUTDIR, "cv_setup.rds"))
message("Salvato ", file.path(OUTDIR, "cv_setup.rds"),
        " -> usalo con low_qc_svm.R quando vuoi aggiungere l'SVM.")

## --- 7. Curve ROC / PR out-of-fold ----------------------------------------
pdf(file.path(OUTDIR, "rf_roc_pr_curves.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
cols <- setNames(c("#1f77b4", "#d62728")[seq_along(EMBEDDINGS)], EMBEDDINGS)

first <- TRUE
for (emb in EMBEDDINGS) {
  d <- oof[oof$embedding == emb, ]
  r <- pROC::roc(d$truth, d$p_RF, levels = c("pass", "low"),
                 direction = "<", quiet = TRUE)
  plot(r, add = !first, col = cols[emb], legacy.axes = TRUE,
       main = "ROC - Random Forest (out-of-fold)")
  first <- FALSE
}
legend("bottomright", names(cols), col = cols, lwd = 2, bty = "n")

first <- TRUE
for (emb in EMBEDDINGS) {
  d <- oof[oof$embedding == emb, ]
  pr <- PRROC::pr.curve(scores.class0 = d$p_RF[d$truth == "low"],
                        scores.class1 = d$p_RF[d$truth == "pass"], curve = TRUE)
  if (first) {
    plot(pr$curve[, 1], pr$curve[, 2], type = "l", ylim = c(0, 1), col = cols[emb],
         xlab = "Recall", ylab = "Precision", main = "Precision-Recall (OOF)")
    first <- FALSE
  } else lines(pr$curve[, 1], pr$curve[, 2], col = cols[emb])
}
legend("bottomleft", names(cols), col = cols, lwd = 2, bty = "n")
dev.off()

## --- 8. Top feature per embedding ------------------------------------------
imp <- purrr::list_rbind(purrr::map2(rf_res, seq_along(rf_res), function(z, i) {
  data.frame(embedding = z$metrics$embedding[1], fold = z$metrics$fold[1],
             feature = names(z$importance), importance = as.numeric(z$importance),
             stringsAsFactors = FALSE)
}))
imp_avg <- aggregate(importance ~ embedding + feature, data = imp, FUN = mean)
imp_top <- do.call(rbind, lapply(split(imp_avg, imp_avg$embedding), function(d)
  head(d[order(-d$importance), ], 15)))
print(imp_top)
write.csv(imp_avg, file.path(OUTDIR, "rf_importance.csv"), row.names = FALSE)

## --- 9. (opzionale) predizioni OOF nel colData -----------------------------
## for (emb in EMBEDDINGS) {
##   v <- rep(NA_real_, ncol(sce)); d <- oof[oof$embedding == emb, ]
##   v[match(d$cell, colnames(sce))] <- d$p_RF
##   colData(sce)[[paste0("p_lowqc_RF_", emb)]] <- v
## }
