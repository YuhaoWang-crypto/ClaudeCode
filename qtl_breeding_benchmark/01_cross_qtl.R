# Cross QTL mapping with R/qtl (HMM genotype probabilities + Haley-Knott interval mapping).
# Data: Attie B6 x BTBR F2 intercross (544 mice), from rqtl/qtl2data.
# (a) known-answer test: two Mendelian coat traits whose causal genes are known;
# (b) discovery/validation: map insulin QTL in half the cross, test them in the unseen half.
suppressMessages(library(qtl)); set.seed(1)
D <- Sys.getenv("QTL_DATA"); OUT <- "results"
g <- read.csv(file.path(D, "b6btbr/b6btbr_geno.csv"), comment.char = "#", na.strings = c("-", "NA"), check.names = FALSE)
p <- read.csv(file.path(D, "b6btbr/b6btbr_pheno.csv"), comment.char = "#")
gm <- read.csv(file.path(D, "b6btbr/b6btbr_gmap.csv"), comment.char = "#")
pm <- read.csv(file.path(D, "b6btbr/b6btbr_pmap.csv"), comment.char = "#")
gm <- gm[gm$chr != "X", ]; stopifnot(all(g[[1]] == p[[1]]))
code <- c(BB = "A", BR = "H", RR = "B")
G <- as.data.frame(lapply(g[, gm$marker], function(x) code[as.character(x)]))
csv <- rbind(c(rep("", 3), gm$chr), c(rep("", 3), gm$pos),
             cbind(p$log10_insulin_10wk, p$agouti_tan, p$tufted, as.matrix(G)))
colnames(csv) <- c("insulin", "agouti", "tufted", gm$marker); csv[is.na(csv)] <- "-"
f <- tempfile(fileext = ".csv"); write.csv(csv, f, row.names = FALSE, quote = FALSE)
cross <- read.cross("csv", file = f, genotypes = c("A", "H", "B"), na.strings = "-", crosstype = "f2")
cross <- calc.genoprob(jittermap(cross), step = 1, error.prob = 0.002)
pmb <- setNames(pm$pos, pm$marker)

# ---- (a) Mendelian traits: is the LOD peak at the known gene?
known <- list(agouti = c(chr = 2, Mb = 154.8, gene = "a (nonagouti)"), tufted = c(chr = 17, Mb = 27.3, gene = "Itpr3"))
res_a <- do.call(rbind, lapply(names(known), function(tr) {
  s <- scanone(cross, pheno.col = tr, method = "hk", model = "binary")
  pk <- s[which.max(s$lod), ]; ci <- lodint(s, chr = as.character(pk$chr), drop = 1.5)
  mk <- find.marker(cross, pk$chr, pk$pos); fl <- find.marker(cross, ci$chr, ci$pos)
  data.frame(trait = tr, known_gene = known[[tr]]["gene"], known_chr = known[[tr]]["chr"], known_Mb = as.numeric(known[[tr]]["Mb"]),
             peak_chr = as.character(pk$chr), peak_cM = pk$pos, peak_LOD = pk$lod, peak_marker_Mb = pmb[mk],
             ci_lo_Mb = min(pmb[fl]), ci_hi_Mb = max(pmb[fl]),
             second_best_chr_LOD = max(tapply(s$lod, s$chr, max)[names(tapply(s$lod, s$chr, max)) != as.character(pk$chr)]))
}))
write.csv(res_a, file.path(OUT, "cross_known_answer.csv"), row.names = FALSE)

# ---- (b) insulin: discovery half -> validation half, 20 random splits
n <- nind(cross); rows <- list(); full <- scanone(cross, pheno.col = "insulin", method = "hk")
write.csv(data.frame(chr = full$chr, cM = full$pos, lod = full$lod), file.path(OUT, "cross_insulin_fullscan.csv"), row.names = FALSE)
thr <- summary(scanone(subset(cross, ind = sample(n, n %/% 2)), pheno.col = "insulin", method = "hk", n.perm = 500), alpha = 0.05)[1]
for (r in 1:20) {
  d <- sample(n, n %/% 2); cd <- subset(cross, ind = d); cv <- subset(cross, ind = setdiff(1:n, d))
  s <- scanone(cd, pheno.col = "insulin", method = "hk"); pk <- summary(s, threshold = thr)
  for (i in seq_len(nrow(pk))) {
    mk <- find.marker(cd, pk$chr[i], pk$pos[i])
    gd <- pull.geno(fill.geno(cd, method = "argmax"))[, mk]; gv <- pull.geno(fill.geno(cv, method = "argmax"))[, mk]
    fd <- summary(lm(cd$pheno$insulin ~ gd)); fv <- summary(lm(cv$pheno$insulin ~ gv))
    rows[[length(rows) + 1]] <- data.frame(split = r, chr = as.character(pk$chr[i]), cM = pk$pos[i], lod_disc = pk$lod[i],
      r2_disc = fd$r.squared, r2_valid = fv$r.squared, p_valid = fv$coefficients[2, 4],
      same_direction = sign(fd$coefficients[2, 1]) == sign(fv$coefficients[2, 1]))
  }
}
b <- do.call(rbind, rows); b$threshold <- thr
write.csv(b, file.path(OUT, "cross_discovery_validation.csv"), row.names = FALSE)
print(res_a); print(thr); print(aggregate(cbind(r2_disc, r2_valid, replicated = p_valid < .05 & same_direction) ~ chr, b, mean)); print(table(b$chr))
