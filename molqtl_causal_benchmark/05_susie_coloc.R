suppressMessages({library(susieR); library(coloc); library(data.table)}); set.seed(1); options(warn=-1)
W <- file.path(Sys.getenv("MQ_WORK"), "regions"); files <- sort(list.files(W, full.names = TRUE))[1:60]
res <- list()
zfun <- function(X, y) { n <- length(y); r <- as.numeric(cor(X, y)); r * sqrt((n - 2)/(1 - r^2 + 1e-12)) }
pick <- function(R, c1, lo, hi) { r2 <- R[c1,]^2; cand <- setdiff(which(r2 >= lo & r2 <= hi), c1); if (length(cand)) cand[sample.int(length(cand),1)] else NA }
for (f in files) {
  d <- tryCatch(fread(f, header=FALSE), error=function(e) NULL); if (is.null(d)) next
  X <- as.matrix(d[,-1]); storage.mode(X) <- "double"; X <- X[, apply(X,2,sd)>0, drop=FALSE]
  p <- ncol(X); n <- nrow(X); if (p < 50 || p > 1200) next
  Xs <- scale(X); R <- cor(X); dimnames(R) <- list(as.character(1:p), as.character(1:p)); maf <- pmin(colMeans(X)/2, 1-colMeans(X)/2)
  sim <- function(cs, pve=0.25) { g <- Xs[,cs,drop=FALSE] %*% rnorm(length(cs)); as.numeric(scale(g)*sqrt(pve)+rnorm(n)*sqrt(1-pve)) }
  D <- function(z) list(beta=z/sqrt(n), varbeta=rep(1/n,p), N=n, sdY=1, type="quant", snp=as.character(1:p), MAF=maf, LD=R)
  run <- function(scen, cs1, cs2, truth) {
    z1 <- zfun(Xs, sim(cs1)); z2 <- zfun(Xs, sim(cs2))
    d1 <- D(z1); d2 <- D(z2)
    abf <- tryCatch(coloc.abf(d1, d2)$summary["PP.H4.abf"], error=function(e) NA)
    s1 <- tryCatch(runsusie(d1, max_iter=100), error=function(e) NULL)
    s2 <- tryCatch(runsusie(d2, max_iter=100), error=function(e) NULL)
    sh <- NA
    if (!is.null(s1) && !is.null(s2)) { cs <- tryCatch(coloc.susie(s1,s2)$summary, error=function(e) NULL)
      if (!is.null(cs) && nrow(cs)>0) sh <- max(cs$PP.H4.abf, na.rm=TRUE) else sh <- 0 }
    data.frame(scenario=scen, truth=truth, abf_H4=as.numeric(abf), susie_H4=as.numeric(sh))
  }
  c1 <- sample(p,1)
  res[[length(res)+1]] <- run("shared_single", c1, c1, "shared")
  cB <- pick(R, c1, 0.2, 0.7); if (!is.na(cB)) res[[length(res)+1]] <- run("two_distinct", c1, cB, "distinct")
  cC <- pick(R, c1, 0.0, 0.1); if (!is.na(cC)) res[[length(res)+1]] <- run("shared_plus_extra", c(c1,cC), c1, "shared")
}
R2 <- do.call(rbind, res); R2 <- R2[!is.na(R2$abf_H4),]
write.csv(R2, "results/susie_coloc_raw.csv", row.names=FALSE)
cat("n runs:", nrow(R2), "\n")
for (sc in c("shared_single","two_distinct","shared_plus_extra")) { d <- R2[R2$scenario==sc,]; if(nrow(d)==0) next
  cat(sprintf("%-18s truth=%-8s | abf: shared-call %.2f, meanH4 %.2f | susie: shared-call %.2f, meanH4 %.2f | n=%d\n",
    sc, d$truth[1], mean(d$abf_H4>0.8,na.rm=T), mean(d$abf_H4,na.rm=T), mean(d$susie_H4>0.8,na.rm=T), mean(d$susie_H4,na.rm=T), nrow(d))) }
cat("WROTE_OK\n")
