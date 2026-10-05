## ---------------------------------------------------------------------------
## Baseline: regressione logistica sulle 4 metriche QC che DEFINISCONO low_qc
## Da eseguire dopo low_qc_classification.R (riusa: keep_cells, y, fold_id,
## eval_pred, oof, EMBEDDINGS, OUTDIR)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(splines)
  library(pROC)
  library(PRROC)
})

## --- 0. Identificazione delle colonne QC -----------------------------------
cd <- colData(sce)
message("Colonne disponibili in colData:")
print(colnames(cd))

## Auto-guess: CONTROLLA e sovrascrivi a mano se sbagliato ------------------

QC_VARS <- c(
  umi      = "sum",
  genes    = "detected",
  mito     = "subsets_Mito_percent",
  intronic = "nuclear_fraction"
)
print(QC_VARS)

## >>> Se l'auto-guess sbaglia, decommenta e correggi: <<<
# QC_VARS <- c(umi = "sum", genes = "detected",
#              mito = "subsets_Mito_percent", intronic = "intronic_fraction")

stopifnot(!any(is.na(QC_VARS)), all(QC_VARS %in% colnames(cd)))

## --- 1. Design matrix ------------------------------------------------------
qc <- data.frame(
  umi      = log1p(cd[[QC_VARS["umi"]]][keep_cells]),       # log per la coda lunga
  genes    = log1p(cd[[QC_VARS["genes"]]][keep_cells]),
  mito     = cd[[QC_VARS["mito"]]][keep_cells],
  intronic = cd[[QC_VARS["intronic"]]][keep_cells]
)
qc$y <- y

message(sprintf("NA nelle QC vars: %d cellule", sum(!complete.cases(qc))))
print(summary(qc[, 1:4]))

## --- patch: gestione NA in intronic ---------------------------------------
qc$intronic_na <- as.integer(is.na(qc$intronic))
qc$intronic[is.na(qc$intronic)] <- 0

## Se le soglie di `keep` sono state calcolate per campione (es. isOutlier con
## batch=), aggiungi il campione come covariata: scommenta le due righe.
# qc$sample <- droplevels(factor(cd$Sample[keep_cells]))
# USE_SAMPLE <- TRUE
USE_SAMPLE <- FALSE

f_lin <- y ~ umi + genes + mito + intronic

## Versione flessibile: low_qc nasce da SOGLIE, quindi la relazione con le
## covariate è fortemente non lineare. Le spline danno il tetto realistico
## raggiungibile con queste 4 informazioni.
f_spl <- y ~ ns(umi, 5) + ns(genes, 5) + ns(mito, 5) + ns(intronic, 5)

## --- 2. Stessi fold della CV precedente ------------------------------------
oof_glm <- do.call(rbind, lapply(seq_len(K), function(f) {
  tr <- fold_id != f; te <- fold_id == f

  m_lin <- glm(f_lin, data = qc[tr, ], family = binomial())
  m_spl <- glm(f_spl, data = qc[tr, ], family = binomial())
  
  
  data.frame(cell    = colnames(sce)[keep_cells][te],
             fold    = f,
             truth   = qc$y[te],
             p_GLM   = predict(m_lin, qc[te, ], type = "response"),
             p_GLMns = predict(m_spl, qc[te, ], type = "response"),
             stringsAsFactors = FALSE)
}))

perf_glm <- rbind(
  cbind(embedding = "QC-only", model = "GLM",
        eval_pred(factor(oof_glm$truth, levels = c("pass","low")), oof_glm$p_GLM)),
  cbind(embedding = "QC-only", model = "GLM-splines",
        eval_pred(factor(oof_glm$truth, levels = c("pass","low")), oof_glm$p_GLMns))
)
print(perf_glm)

## --- 2b. Baseline casuale (no-information rate) ----------------------------
## Score uniforme + soglia alla prevalenza = tiro a caso stratificato.
## Ripetuto N_REP volte perche' una singola estrazione e' rumorosa.
set.seed(1)
N_REP <- 50
prev  <- mean(y == "low")

chance_rep <- do.call(rbind, lapply(seq_len(N_REP), function(r)
  eval_pred(y, runif(length(y)), thr = 1 - prev)))

perf_chance <- cbind(embedding = "Chance", model = "Random stratificato",
                     as.data.frame(t(colMeans(chance_rep))))
print(round(rbind(media = colMeans(chance_rep),
                  sd    = apply(chance_rep, 2, sd)), 4))

message(sprintf("Attesi teorici: AUROC 0.500 | AUPRC %.3f (= prevalenza) | bal.acc 0.500 | MCC 0",
                prev))

## --- 2c. Baseline PCA: stesso gioco, ma con una riduzione classica ---------
## Domanda: i foundation model battono una semplice PCA sull'espressione?
N_PC  <- min(50, ncol(reducedDim(sce, "PCA")))
X_pca <- as.matrix(reducedDim(sce, "PCA"))[keep_cells, seq_len(N_PC), drop = FALSE]
colnames(X_pca) <- paste0("PC", seq_len(N_PC))
message(sprintf("PCA: %d componenti su %d disponibili",
                N_PC, ncol(reducedDim(sce, "PCA"))))

oof_pca <- do.call(rbind, lapply(seq_len(K), function(f) {
  tr <- fold_id != f; te <- fold_id == f
  
  ctr <- colMeans(X_pca[tr, , drop = FALSE])
  sdv <- matrixStats::colSds(X_pca[tr, , drop = FALSE]); sdv[sdv == 0] <- 1
  Ztr <- as.data.frame(scale(X_pca[tr, , drop = FALSE], center = ctr, scale = sdv))
  Zte <- as.data.frame(scale(X_pca[te, , drop = FALSE], center = ctr, scale = sdv))
  
  m_glm <- glm(y ~ ., data = cbind(Ztr, y = y[tr]), family = binomial())
  m_rf  <- ranger::ranger(x = as.matrix(Ztr), y = y[tr], num.trees = RF_TREES,
                          mtry = floor(sqrt(N_PC)), probability = TRUE,
                          num.threads = RF_THREADS, seed = 300 + f)
  
  data.frame(cell      = colnames(sce)[keep_cells][te],
             fold      = f,
             truth     = y[te],
             p_PCA_GLM = predict(m_glm, Zte, type = "response"),
             p_PCA_RF  = predict(m_rf, as.matrix(Zte),
                                 num.threads = RF_THREADS)$predictions[, "low"],
             stringsAsFactors = FALSE)
}))

perf_pca <- rbind(
  cbind(embedding = "PCA", model = sprintf("GLM (%d PC)", N_PC),
        eval_pred(oof_pca$truth, oof_pca$p_PCA_GLM)),
  cbind(embedding = "PCA", model = sprintf("RF (%d PC)", N_PC),
        eval_pred(oof_pca$truth, oof_pca$p_PCA_RF))
)
print(perf_pca)

## --- 3. Confronto con i modelli sugli embeddings ---------------------------
## oof_perf viene dallo script precedente (predizioni out-of-fold aggregate)
compare <- rbind(oof_perf, perf_glm, perf_pca, perf_chance)
compare <- compare[order(-compare$auprc), ]
print(compare[, c("embedding", "model", "precision", "f1",
                  "sensitivity", "specificity")])
write.csv(compare, file.path(OUTDIR, "comparison_embeddings_vs_qc.csv"),
          row.names = FALSE) 

## --- 4. Test statistico sulle ROC (DeLong, curve appaiate) -----------------
## Stesse cellule -> curve correlate: roc.test con method="delong" è appropriato.
oof_all <- merge(oof_glm[, c("cell", "p_GLM", "p_GLMns")],
                 oof, by = "cell")
roc_of <- function(score, truth)
  pROC::roc(truth, score, levels = c("pass", "low"), direction = "<", quiet = TRUE)

delong <- do.call(rbind, lapply(EMBEDDINGS, function(emb) {
  d <- oof_all[oof_all$embedding == emb, ]
  tr <- factor(d$truth, levels = c("pass", "low"))
  base <- roc_of(d$p_GLMns, tr)
  do.call(rbind, lapply(c("p_RF", "p_SVM"), function(m) {
    tt <- pROC::roc.test(roc_of(d[[m]], tr), base, method = "delong",
                         paired = TRUE)
    data.frame(embedding = emb, model = sub("p_", "", m),
               auc_model = as.numeric(tt$estimate[1]),
               auc_qc_splines = as.numeric(tt$estimate[2]),
               delta = as.numeric(diff(rev(tt$estimate))),
               p_value = tt$p.value)
  }))
}))
print(delong)
write.csv(delong, file.path(OUTDIR, "delong_tests.csv"), row.names = FALSE)

## --- 5. Valore incrementale: embedding + QC vs QC da solo ------------------
## Risponde a: "gli embeddings aggiungono segnale OLTRE alle 4 metriche?"
incr <- do.call(rbind, lapply(EMBEDDINGS, function(emb) {
  X <- X_list[[emb]]
  do.call(rbind, lapply(seq_len(K), function(f) {
    tr <- fold_id != f; te <- fold_id == f
    Xtr <- cbind(X[tr, , drop = FALSE], as.matrix(qc[tr, 1:4]))
    Xte <- cbind(X[te, , drop = FALSE], as.matrix(qc[te, 1:4]))
    rf <- ranger(x = Xtr, y = y[tr], num.trees = RF_TREES,
                 probability = TRUE, num.threads = RF_THREADS, seed = 200 + f)
    data.frame(embedding = emb, fold = f,
               cell = colnames(sce)[keep_cells][te], truth = y[te],
               p = predict(rf, Xte, num.threads = RF_THREADS)$predictions[, "low"],
               stringsAsFactors = FALSE)
  }))
}))

perf_incr <- do.call(rbind, lapply(EMBEDDINGS, function(emb) {
  d <- incr[incr$embedding == emb, ]
  cbind(embedding = emb, model = "RF (embedding + QC)",
        eval_pred(factor(d$truth, levels = c("pass", "low")), d$p))
}))
print(perf_incr)

write.csv(rbind(compare, perf_incr),
          file.path(OUTDIR, "comparison_full.csv"), row.names = FALSE)

## --- 6. Curve ROC di confronto --------------------------------------------
pdf(file.path(OUTDIR, "roc_embeddings_vs_qc.pdf"), width = 6, height = 6)
d <- oof_all[oof_all$embedding == EMBEDDINGS[1], ]
tr <- factor(d$truth, levels = c("pass", "low"))
plot(roc_of(d$p_GLMns, tr), col = "black", lwd = 2, legacy.axes = TRUE,
     main = "low_qc: embeddings vs 4 metriche QC (out-of-fold)")
plot(roc_of(d$p_GLM, tr), add = TRUE, col = "grey60", lwd = 2, lty = 2)
pal <- c("#1f77b4", "#7fbfe8", "#d62728", "#f2a0a1"); j <- 1
for (emb in EMBEDDINGS) for (m in c("p_RF", "p_SVM")) {
  dd <- oof_all[oof_all$embedding == emb, ]
  plot(roc_of(dd[[m]], factor(dd$truth, levels = c("pass", "low"))),
       add = TRUE, col = pal[j], lwd = 2); j <- j + 1
}
legend("bottomright", bty = "n", lwd = 2, cex = 0.8,
       col = c("black", "grey60", pal), lty = c(1, 2, rep(1, 4)),
       legend = c("QC splines", "QC lineare",
                  paste(rep(EMBEDDINGS, each = 2), c("RF", "SVM"))))
dev.off()
