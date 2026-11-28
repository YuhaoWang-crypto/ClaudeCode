# For each LDL locus, run coloc (eQTL vs LDL GWAS) for EVERY gene with an eQTL, in liver and blood.
# Question: does colocalization uniquely identify the validated causal gene, or light up many genes?
suppressMessages(library(coloc)); options(warn=-1)
W <- Sys.getenv("CW_DIR"); OUT <- "results"
loci <- jsonlite::fromJSON(file.path(OUT, "causal_loci.json"))
Nqtl <- c(liver = 208, blood = 670); Ngwas <- 1320000
res <- list()
for (g in names(loci)) {
  gw <- tryCatch(read.delim(file.path(W, "loci", paste0(g, ".gwas.tsv"))), error = function(e) NULL)
  if (is.null(gw)) next
  gw$key <- paste(gw$pos, pmin(gw$ea, gw$oa), pmax(gw$ea, gw$oa)); gw <- gw[!duplicated(gw$key) & !is.na(gw$beta) & gw$se > 0, ]
  for (tis in c("liver", "blood")) {
    f <- file.path(W, "loci", paste0(g, ".", tis, ".tsv")); if (!file.exists(f)) next
    eq <- read.delim(f); eq <- eq[!is.na(eq$se) & eq$se > 0 & eq$maf > 0 & eq$maf < 1, ]
    eq$key <- paste(eq$pos, pmin(eq$ref, eq$alt), pmax(eq$ref, eq$alt))
    for (gene in unique(eq$gene_id)) {
      e <- eq[eq$gene_id == gene, ]; e <- e[!duplicated(e$key), ]
      m <- merge(e, gw, by = "key", suffixes = c(".e", ".g")); if (nrow(m) < 50) next
      r <- tryCatch(suppressWarnings(coloc.abf(
        list(beta = m$beta.e, varbeta = m$se.e^2, N = Nqtl[tis], MAF = m$maf, type = "quant", snp = m$key),
        list(beta = as.numeric(m$beta.g), varbeta = as.numeric(m$se.g)^2, N = Ngwas, MAF = m$maf, type = "quant", snp = m$key)))$summary,
        error = function(e) NULL)
      if (is.null(r)) next
      res[[length(res) + 1]] <- data.frame(locus = g, tissue = tis, gene_id = gene,
        gene_symbol = e$gene[1], n_snps = nrow(m), min_p_eqtl = min(m$p.e), min_p_gwas = min(as.numeric(m$p.g)),
        PP_H4 = as.numeric(r["PP.H4.abf"]), PP_H3 = as.numeric(r["PP.H3.abf"]),
        is_causal = as.integer(gene == loci[[g]]$ensg))
    }
  }
}
R <- do.call(rbind, res); write.csv(R, file.path(OUT, "coloc_all.csv"), row.names = FALSE)
cat("total gene-locus-tissue coloc tests:", nrow(R), "\n")
print(table(R$tissue))
