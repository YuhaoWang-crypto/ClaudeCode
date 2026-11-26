# Locked expression prognostic panel: train penalized Cox on 60% of a TCGA cohort, lock it,
# score the held-out 40% and the other cancers. Tests the matrix-row boundary: a locked panel
# yields a per-sample score, but population value (C-index, incremental value) is optimistic and
# does not transfer across cohorts.
suppressMessages({library(survival); library(glmnet); library(data.table)}); set.seed(1)
T <- Sys.getenv("TB_DATA"); OUT <- "results"
expr <- as.data.frame(fread(file.path(T, "expr_sub.tsv.gz"))); rownames(expr) <- expr[[1]]; expr[[1]] <- NULL
v <- apply(expr, 1, var); expr <- expr[order(-v)[1:2500], ]                 # top 2500 variable genes
clin <- as.data.frame(fread(file.path(T, "clin_sub.tsv"))); rownames(clin) <- clin$sample
cidx <- function(time, status, risk) {
  if (length(unique(status)) < 2 || sd(risk, na.rm = TRUE) == 0 || sum(is.finite(risk)) < 10) return(NA)
  as.numeric(survival::concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance) }
train_panel <- function(samp) {
  X <- t(as.matrix(expr[, samp])); mu <- colMeans(X); sdv <- apply(X, 2, sd); sdv[sdv == 0] <- 1
  Xs <- scale(X, mu, sdv); y <- Surv(clin[samp, "OS.time"], clin[samp, "OS"])
  cv <- cv.glmnet(Xs, y, family = "cox", alpha = 0.3, nfolds = 5)
  b <- as.numeric(coef(cv, s = "lambda.min")); names(b) <- colnames(Xs)
  list(mu = mu, sdv = sdv, b = b, genes = names(b)[b != 0]) }
score <- function(m, samp) as.numeric(scale(t(as.matrix(expr[, samp])), m$mu, m$sdv) %*% m$b)

cohorts <- c("LUAD", "LUSC", "KIRC", "BRCA"); res <- list(); transfer <- list(); nulls <- c()
for (co in cohorts) {
  s <- clin$sample[clin$cancer == co]; ev <- clin[s, "OS"]
  for (rep in 1:3) {
    tr <- c(sample(s[ev == 1], round(.6 * sum(ev == 1))), sample(s[ev == 0], round(.6 * sum(ev == 0))))
    te <- setdiff(s, tr); m <- train_panel(tr)
    rk_tr <- score(m, tr); rk_te <- score(m, te); cl <- clin[te, ]
    ok <- !is.na(cl$age) & !is.na(cl$stage_num)
    cclin <- coxph(Surv(OS.time, OS) ~ age + stage_num, clin[tr, ])
    c_clin <- cidx(cl$OS.time[ok], cl$OS[ok], predict(cclin, cl[ok, ]))
    trd <- data.frame(OS.time = clin[tr, "OS.time"], OS = clin[tr, "OS"], age = clin[tr, "age"], st = clin[tr, "stage_num"], panel = rk_tr)
    cj <- coxph(Surv(OS.time, OS) ~ age + st + panel, trd)
    ted <- data.frame(OS.time = cl$OS.time, OS = cl$OS, age = cl$age, st = cl$stage_num, panel = rk_te)[ok, ]
    c_both <- cidx(ted$OS.time, ted$OS, predict(cj, ted))
    res[[length(res) + 1]] <- data.frame(cancer = co, rep = rep, n_train = length(tr), n_test = length(te),
      panel_size = length(m$genes), c_train = cidx(clin[tr, "OS.time"], clin[tr, "OS"], rk_tr),
      c_test_panel = cidx(cl$OS.time, cl$OS, rk_te), c_test_clinical = c_clin, c_test_both = c_both)
    for (co2 in setdiff(cohorts, co)) { s2 <- clin$sample[clin$cancer == co2]
      transfer[[length(transfer) + 1]] <- data.frame(train = co, test = co2, rep = rep,
        c_index = cidx(clin[s2, "OS.time"], clin[s2, "OS"], score(m, s2))) }
  }
  # one permuted-label sanity panel per cohort
  tr <- c(sample(s[ev == 1], round(.6 * sum(ev == 1))), sample(s[ev == 0], round(.6 * sum(ev == 0)))); te <- setdiff(s, tr)
  clp <- clin; clp[tr, "OS"] <- sample(clp[tr, "OS"]); clp[tr, "OS.time"] <- sample(clp[tr, "OS.time"])
  assign("clin_bak", clin); clin <<- clp; mp <- train_panel(tr); clin <<- clin_bak
  nulls <- c(nulls, cidx(clin[te, "OS.time"], clin[te, "OS"], score(mp, te)))
}
R <- do.call(rbind, res); TR <- do.call(rbind, transfer)
fwrite(R, file.path(OUT, "prognostic_splits.csv")); fwrite(TR, file.path(OUT, "prognostic_transfer.csv"))
fwrite(data.frame(cancer = cohorts, null_permuted_c_index = round(nulls, 3)), file.path(OUT, "prognostic_null.csv"))
print(aggregate(cbind(c_train, c_test_panel, c_test_clinical, c_test_both, panel_size) ~ cancer, R, function(x) round(mean(x), 3)))
cat("null permuted C-index (want ~0.5):", round(nulls, 3), "\n")
print(aggregate(c_index ~ train + test, TR, function(x) round(mean(x), 3)))
cat("WROTE_OK\n")
