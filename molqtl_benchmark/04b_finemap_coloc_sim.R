# Known-truth simulations on real GEUVADIS genotypes (real LD, simulated phenotypes).
# A. SuSiE fine-mapping: do 95% credible sets contain the causal variant, and what happens when the
#    LD reference comes from the wrong ancestry?
# B. coloc: how often does it call "shared causal variant" when the two traits truly share one,
#    and when they have two different causal variants in LD?
suppressMessages({library(susieR); library(coloc); library(data.table)}); set.seed(1)
W <- file.path(Sys.getenv("MQ_WORK"), "regions"); files <- sort(list.files(W, full.names = TRUE))
fm <- list(); cl <- list()
zscore <- function(X, y) { n <- length(y); r <- as.numeric(cor(X, y)); r * sqrt((n - 2) / (1 - r^2)) }
eval_cs <- function(fit, causal, label, k, R) { cs <- fit$sets$cs
  hit <- if (length(cs)) sapply(cs, function(s) any(causal %in% s)) else logical(0)
  # how many causal variants are tagged (r2 > 0.8) by the top-PIP SNP of some CS
  data.frame(method = label, k = k, n_cs = length(cs), cs_with_causal = sum(hit), causal_found = sum(sapply(causal, function(c) any(sapply(cs, function(s) c %in% s)))),
             median_cs_size = if (length(cs)) median(sapply(cs, length)) else NA, top_pip_is_causal = which.max(fit$pip) %in% causal) }
for (f in files) {
  d <- fread(f, header = FALSE); set <- d[[1]]; X <- as.matrix(d[, -1]); storage.mode(X) <- "double"
  eur <- set != "YRI"; tr <- set == "EUR_train"; te <- set == "EUR_test"; yri <- set == "YRI"
  keep <- apply(X[tr, ], 2, sd) > 0 & apply(X[te, ], 2, sd) > 0 & apply(X[yri, ], 2, sd) > 0; X <- X[, keep]; p <- ncol(X)
  Xs <- scale(X[eur, ]); R_in <- cor(X[tr, ]); R_eur <- cor(X[te, ]); R_yri <- cor(X[yri, ])
  # ---------- A. fine-mapping, GWAS-like study = EUR_train (n = 261)
  for (k in 1:3) { causal <- sample(p, k); b <- rnorm(k) ; g <- as.numeric(scale(X[tr, causal, drop = FALSE]) %*% b)
    y <- as.numeric(scale(g)) * sqrt(0.4) + rnorm(sum(tr)) * sqrt(0.6); z <- zscore(X[tr, ], y)
    fit <- tryCatch(susie(X[tr, ], y, L = 10), error = function(e) NULL)
    if (!is.null(fit)) fm[[length(fm) + 1]] <- eval_cs(fit, causal, "individual-level data", k)
    for (nm in c("in-sample LD", "external LD, same ancestry (112 EUR)", "external LD, other ancestry (89 YRI)")) {
      R <- list(R_in, R_eur, R_yri)[[match(nm, c("in-sample LD", "external LD, same ancestry (112 EUR)", "external LD, other ancestry (89 YRI)"))]]
      fit <- tryCatch(suppressWarnings(susie_rss(z, R, n = sum(tr), L = 10)), error = function(e) NULL)
      if (!is.null(fit)) fm[[length(fm) + 1]] <- eval_cs(fit, causal, paste("summary stats +", nm), k) } }
  # ---------- B. colocalisation: eQTL on real people (n = 373), GWAS summary statistics drawn with the same LD
  n <- sum(eur); R <- cor(X[eur, ]); maf <- colMeans(X[eur, ]) / 2; maf <- pmin(maf, 1 - maf)
  gwas <- function(c2, zexp = 7) as.numeric(crossprod(Xs, rnorm(n)) / sqrt(n - 1)) + R[, c2] * zexp   # z ~ N(R*lambda, R)
  eqtl <- function(c1, pve) { y <- Xs[, c1] * sqrt(pve) + rnorm(n) * sqrt(1 - pve); r <- as.numeric(cor(Xs, y)); list(b = r, v = (1 - r^2) / (n - 2)) }
  run <- function(scen, c1, c2, pve = 0.15) { e <- if (is.na(c1)) eqtl(1, 0) else eqtl(c1, pve); z <- gwas(c2); N <- 50000
    s <- suppressWarnings(suppressMessages(coloc.abf(list(beta = e$b, varbeta = e$v, N = n, sdY = 1, type = "quant", snp = 1:p, MAF = maf),
      list(beta = z / sqrt(N), varbeta = rep(1 / N, p), N = N, sdY = 1, type = "quant", snp = 1:p, MAF = maf))))$summary
    cl[[length(cl) + 1]] <<- data.frame(scenario = scen, H3 = s["PP.H3.abf"], H4 = s["PP.H4.abf"]) }
  sink(tempfile())
  c1 <- sample(p, 1); r2 <- R[c1, ]^2
  run("same causal variant, strong eQTL", c1, c1); run("same causal variant, weak eQTL (4% of variance)", c1, c1, 0.04)
  for (b in list(c("different variants, r2 < 0.1", 0, .1), c("different variants, r2 0.3-0.7", .3, .7), c("different variants, r2 0.8-0.95", .8, .95))) {
    cand <- setdiff(which(r2 >= as.numeric(b[2]) & r2 <= as.numeric(b[3])), c1); if (length(cand)) run(b[1], c1, cand[sample.int(length(cand), 1)]) }
  run("GWAS signal only, no eQTL", NA, c1)
  sink()
}
fm <- do.call(rbind, fm); cl <- do.call(rbind, cl)
A <- do.call(rbind, lapply(split(fm, list(fm$method, fm$k), drop = TRUE), function(d) data.frame(method = d$method[1], n_causal = d$k[1], runs = nrow(d),
  cs_per_run = round(mean(d$n_cs), 2), coverage = round(sum(d$cs_with_causal) / sum(d$n_cs), 3), false_cs_per_run = round(mean(d$n_cs - d$cs_with_causal), 3),
  power = round(sum(d$causal_found) / sum(d$k), 3), median_cs_size = median(d$median_cs_size, na.rm = TRUE), top_pip_is_causal = round(mean(d$top_pip_is_causal), 3))))
B <- do.call(rbind, lapply(split(cl, cl$scenario), function(d) data.frame(scenario = d$scenario[1], runs = nrow(d), mean_PP_H4 = round(mean(d$H4), 3),
  called_shared_H4_gt_0.8 = round(mean(d$H4 > .8), 3), called_distinct_H3_gt_0.8 = round(mean(d$H3 > .8), 3))))
rownames(A) <- rownames(B) <- NULL; A <- A[order(A$n_causal, A$method), ]
write.csv(A, "results/finemap_sim.csv", row.names = FALSE); write.csv(B, "results/coloc_sim.csv", row.names = FALSE); print(A); print(B)
