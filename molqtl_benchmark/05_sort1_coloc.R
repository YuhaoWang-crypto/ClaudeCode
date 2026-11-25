# Real known-answer locus: 1p13 (SORT1 / PSRC1 / CELSR2) and LDL cholesterol.
# Known biology (Musunuru et al. 2010): rs12740374 creates a liver-specific enhancer; the LDL-lowering allele
# raises hepatic SORT1 expression. Expectation: colocalisation in liver, not in non-hepatic tissues,
# and a negative expression -> LDL direction.
suppressMessages(library(coloc)); D <- Sys.getenv("MQ_DATA")
gw <- read.delim(file.path(D, "ldl_sort1.tsv"), header = FALSE)[, c(1:8, 10)]
names(gw) <- c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "n")
gw$id <- paste(gw$pos, pmin(gw$ea, gw$oa), pmax(gw$ea, gw$oa)); gw <- gw[!duplicated(gw$id) & !is.na(gw$beta) & gw$se > 0, ]
tissues <- c(liver = "QTD000266", blood = "QTD000356", LCL = "QTD000221", muscle = "QTD000281", adipose = "QTD000116")
gene <- c(SORT1 = "ENSG00000134243", PSRC1 = "ENSG00000134222", CELSR2 = "ENSG00000143126")
out <- list()
for (t in names(tissues)) {
  e <- read.delim(file.path(D, sprintf("eqtl_%s_sort1.tsv", tissues[t])), header = FALSE)
  names(e)[c(1, 3:5, 8:11, 14, 19)] <- c("gene", "pos", "ref", "alt", "maf", "p", "beta", "se", "an", "rsid")
  for (g in names(gene)) {
    x <- e[e$gene == gene[g] & !is.na(e$se) & e$se > 0, ]; if (nrow(x) < 50) next
    x$id <- paste(x$pos, pmin(x$ref, x$alt), pmax(x$ref, x$alt)); x <- x[!duplicated(x$id), ]
    m <- merge(x, gw, by = "id", suffixes = c(".e", ".g")); n_e <- max(x$an) / 2
    r <- suppressWarnings(coloc.abf(
      list(beta = m$beta.e, varbeta = m$se.e^2, N = n_e, MAF = m$maf, type = "quant", snp = m$id),
      list(beta = m$beta.g, varbeta = m$se.g^2, N = median(m$n), MAF = m$maf, type = "quant", snp = m$id)))$summary
    le <- m[which.min(m$p.e), ]                       # lead eQTL variant; align GWAS effect to the eQTL alt allele
    bg <- ifelse(le$ea == le$alt, le$beta.g, -le$beta.g)
    out[[length(out) + 1]] <- data.frame(tissue = t, gene = g, n_eqtl = n_e, snps = nrow(m), min_p_eqtl = min(m$p.e), min_p_gwas = min(m$p.g),
      PP_H3_distinct = round(r["PP.H3.abf"], 3), PP_H4_shared = round(r["PP.H4.abf"], 3),
      lead_eqtl = le$rsid, lead_eqtl_p = le$p.e, lead_eqtl_gwas_p = le$p.g,
      wald_ratio_LDL_per_SD_expression = round(bg / le$beta.e, 3), wald_se = round(abs(le$se.g / le$beta.e), 3))
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, "results/sort1_coloc.csv", row.names = FALSE); print(res[, c(1:3, 5:10, 12)])
