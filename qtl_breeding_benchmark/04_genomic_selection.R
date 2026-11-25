# Genomic selection: train marker weights on public lines, predict GEBV for lines the model never saw.
# Wheat: CIMMYT 599 lines x 1279 DArT markers x 4 environments (ships with BGLR).
# Rice : Rice Diversity Panel, 413 accessions x ~36K SNPs (prepared by 02_rice_gwas.sh).
suppressMessages({library(rrBLUP); library(BGLR)}); set.seed(1)
W <- file.path(Sys.getenv("QTL_WORK"), "rice"); D <- Sys.getenv("QTL_DATA"); OUT <- "results"
res <- list(); add <- function(...) { res[[length(res) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
  write.csv(do.call(rbind, res), file.path(OUT, "genomic_selection.csv"), row.names = FALSE) }
gblup <- function(K, y, test) { yy <- y; yy[test] <- NA
  m <- mixed.solve(yy, K = K); as.numeric(m$beta) + m$u[test] }
bayesB <- function(X, y, test) { yy <- y; yy[test] <- NA
  f <- BGLR(yy, ETA = list(list(X = X, model = "BayesB")), nIter = 3000, burnIn = 1000, verbose = FALSE,
            saveAt = file.path(tempdir(), "b_")); f$yHat[test] }
folds <- function(n, k = 5) sample(rep(1:k, length.out = n))
cc <- function(a, b) if (sum(complete.cases(a, b)) > 9) cor(a, b, use = "complete.obs") else NA

# ================= wheat =================
data(wheat); X <- scale(wheat.X); K <- tcrossprod(X) / ncol(X); n <- nrow(X); E <- paste0("env", 1:4)
for (rep in 1:3) { fo <- folds(n)
  pg <- matrix(NA, n, 4); pb <- matrix(NA, n, 4)
  for (k in 1:5) { te <- which(fo == k)
    for (e in 1:4) { pg[te, e] <- gblup(K, wheat.Y[, e], te); if (rep == 1) pb[te, e] <- bayesB(X, wheat.Y[, e], te) } }
  for (e in 1:4) { # per-fold correlations, then averaged
    add(data = "wheat", scenario = "unseen lines, same environment", model = "GBLUP", trait = E[e], target = E[e], rep = rep,
        r = mean(sapply(1:5, function(k) cc(pg[fo == k, e], wheat.Y[fo == k, e]))), n_test = n / 5)
    if (rep == 1) add(data = "wheat", scenario = "unseen lines, same environment", model = "BayesB", trait = E[e], target = E[e], rep = rep,
        r = mean(sapply(1:5, function(k) cc(pb[fo == k, e], wheat.Y[fo == k, e]))), n_test = n / 5)
    for (e2 in setdiff(1:4, e)) add(data = "wheat", scenario = "unseen lines, different environment", model = "GBLUP", trait = E[e], target = E[e2],
        rep = rep, r = mean(sapply(1:5, function(k) cc(pg[fo == k, e], wheat.Y[fo == k, e2]))), n_test = n / 5) } }
cat("wheat done\n")

# ================= rice =================
raw <- read.table(file.path(W, "panel.raw"), header = TRUE, check.names = FALSE)
M <- as.matrix(raw[, -(1:6)]); M <- apply(M, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
bim <- read.table(file.path(W, "panel.bim")); Xs <- scale(M); Xs[is.na(Xs)] <- 0; K <- tcrossprod(Xs) / ncol(Xs)
ph <- read.delim(file.path(D, "rice/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"), check.names = FALSE)
names(ph) <- trimws(names(ph)); ph <- ph[match(raw$IID, ph$NSFTVID), ]
gp <- read.csv(file.path(D, "rice/germplasm.csv"), skip = 1); sub <- as.character(gp$Sub.population[match(raw$IID, gp$NSFTV.ID)])
main <- c("IND", "TEJ", "TRJ", "AUS"); n <- nrow(M)
within_sub <- function(pred, y, idx) { v <- sapply(main, function(s) cc(pred[idx & sub == s], y[idx & sub == s])); mean(v, na.rm = TRUE) }
top_markers <- function(y, rows, k = 10) { # greedy: best marginal SNP, then next best >= 1 Mb away
  r2 <- as.numeric(cor(M[rows, ], y[rows], use = "pairwise.complete.obs"))^2; r2[is.na(r2)] <- 0; sel <- c()
  for (j in order(-r2)) { if (length(sel) == k) break
    if (all(bim$V1[sel] != bim$V1[j] | abs(bim$V4[sel] - bim$V4[j]) > 1e6)) sel <- c(sel, j) }; sel }
traits <- c("Flowering time at Arkansas", "Plant height", "Seed length", "Amylose content")
for (tr in traits) { y <- ph[[tr]]; ok <- !is.na(y) & !is.na(sub)
  for (rep in 1:3) { fo <- rep(NA, n); fo[ok] <- folds(sum(ok)); p <- pk <- pl <- rep(NA, n)
    leak <- top_markers(y, which(ok))                      # WRONG on purpose: markers chosen with test lines included
    for (k in 1:5) { te <- which(fo == k); tr_rows <- which(ok & fo != k)
      yy <- y; yy[!ok] <- NA; p[te] <- gblup(K, yy, c(te, which(!ok)))[seq_along(te)]
      s <- top_markers(y, tr_rows); d <- data.frame(y = y, M[, s]); pk[te] <- predict(lm(y ~ ., d[tr_rows, ]), d[te, ])
      d <- data.frame(y = y, M[, leak]); pl[te] <- predict(lm(y ~ ., d[tr_rows, ]), d[te, ]) }
    fr <- function(v) mean(sapply(1:5, function(k) cc(v[which(fo == k)], y[which(fo == k)])))
    add(data = "rice", scenario = "unseen lines, random split (whole panel r)", model = "GBLUP", trait = tr, target = tr, rep = rep, r = fr(p), n_test = sum(ok) / 5)
    add(data = "rice", scenario = "unseen lines, random split (r within subpopulation)", model = "GBLUP", trait = tr, target = tr, rep = rep, r = within_sub(p, y, ok), n_test = NA)
    add(data = "rice", scenario = "10 top markers chosen in training fold", model = "OLS on 10 SNPs", trait = tr, target = tr, rep = rep, r = fr(pk), n_test = sum(ok) / 5)
    add(data = "rice", scenario = "10 top markers chosen in training fold (r within subpopulation)", model = "OLS on 10 SNPs", trait = tr, target = tr, rep = rep, r = within_sub(pk, y, ok), n_test = NA)
    add(data = "rice", scenario = "10 top markers chosen on ALL lines (leaky)", model = "OLS on 10 SNPs", trait = tr, target = tr, rep = rep, r = fr(pl), n_test = sum(ok) / 5)
    if (tr == "Flowering time at Arkansas") for (t2 in c("Flowering time at Faridpur", "Flowering time at Aberdeen"))
      add(data = "rice", scenario = "unseen lines, different environment", model = "GBLUP", trait = tr, target = t2, rep = rep,
          r = mean(sapply(1:5, function(k) cc(p[which(fo == k)], ph[[t2]][which(fo == k)]))), n_test = NA) }
  for (s in main) { te <- which(ok & sub == s); yy <- y; yy[!ok] <- NA   # whole subpopulation never seen in training
    pr <- gblup(K, yy, c(te, which(!ok)))[seq_along(te)]
    add(data = "rice", scenario = "unseen subpopulation", model = "GBLUP", trait = tr, target = s, rep = 1, r = cc(pr, y[te]), n_test = length(te)) }
  cat(tr, "done\n") }
a <- do.call(rbind, res); print(aggregate(r ~ data + scenario + model + trait, a, function(v) round(mean(v), 3)))
