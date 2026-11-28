# pQTL (deCODE/AGES SomaScan, PMID 35078996) vs LDL GWAS colocalization.
# Tests whether pQTL recovers the protein-mechanism causal genes that cis-eQTL coloc missed.
suppressMessages({ .libPaths(c(Sys.getenv("RLIB"), .libPaths())); library(coloc) })
args <- commandArgs(trailingOnly=TRUE)
PQ <- args[1]; LOCI <- args[2]; OUT <- args[3]
N_PQTL <- 35559    # Ferkingstad deCODE plasma proteomics
N_GWAS <- 1320016  # GLGC LDL-C
genes <- sub("\\.pqtl\\.tsv$","", list.files(PQ, pattern="\\.pqtl\\.tsv$"))
rows <- list()
for (g in genes) {
  pf <- file.path(PQ, paste0(g,".pqtl.tsv"))
  gf <- file.path(LOCI, paste0(g,".gwas.tsv"))
  if (!file.exists(gf)) { message(g," no gwas"); next }
  pq <- read.table(pf, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  gw <- read.table(gf, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  pq <- pq[!is.na(pq$hm_beta) & !is.na(pq$standard_error) & !is.na(pq$hm_effect_allele_frequency) & pq$standard_error>0, ]
  pq$key <- paste(pq$hm_pos, pmin(pq$hm_effect_allele,pq$hm_other_allele), pmax(pq$hm_effect_allele,pq$hm_other_allele), sep="_")
  gw <- gw[!is.na(gw$beta) & !is.na(gw$se) & !is.na(gw$eaf) & gw$se>0, ]
  gw$key <- paste(gw$pos, pmin(gw$ea,gw$oa), pmax(gw$ea,gw$oa), sep="_")
  pq <- pq[!duplicated(pq$key), ]; gw <- gw[!duplicated(gw$key), ]
  m <- intersect(pq$key, gw$key)
  if (length(m) < 50) { message(g," only ",length(m)," shared snps"); next }
  pq <- pq[match(m, pq$key), ]; gw <- gw[match(m, gw$key), ]
  # align pQTL effect to GWAS effect allele
  flip <- ifelse(pq$hm_effect_allele == gw$ea, 1, -1)
  pq_beta <- pq$hm_beta * flip
  maf <- pmin(pq$hm_effect_allele_frequency, 1-pq$hm_effect_allele_frequency)
  ok <- maf>0.001 & maf<0.999
  m<-m[ok]; pq<-pq[ok,]; gw<-gw[ok,]; pq_beta<-pq_beta[ok]; maf<-maf[ok]
  D1 <- list(beta=pq_beta, varbeta=pq$standard_error^2, snp=m, MAF=maf, N=N_PQTL, type="quant")
  gmaf <- pmin(gw$eaf,1-gw$eaf); gmaf[gmaf<=0|gmaf>=1]<-NA
  D2 <- list(beta=gw$beta, varbeta=gw$se^2, snp=m, MAF=ifelse(is.na(gmaf),maf,gmaf), N=N_GWAS, type="quant")
  r <- tryCatch(coloc.abf(D1, D2), error=function(e){message(g," coloc err ",conditionMessage(e)); NULL})
  if (is.null(r)) next
  pp <- r$summary
  rows[[g]] <- data.frame(gene=g, nsnp=length(m),
    pqtl_minP=min(pchisq((pq_beta/pq$standard_error)^2,1,lower.tail=FALSE)),
    PP_H0=pp[["PP.H0.abf"]], PP_H1=pp[["PP.H1.abf"]], PP_H2=pp[["PP.H2.abf"]],
    PP_H3=pp[["PP.H3.abf"]], PP_H4=pp[["PP.H4.abf"]])
  message(sprintf("%-7s nsnp=%d pQTL_minP=%.1e PP.H4=%.3f", g, length(m), rows[[g]]$pqtl_minP, pp[["PP.H4.abf"]]))
}
res <- do.call(rbind, rows)
write.csv(res, OUT, row.names=FALSE)
cat("written ", OUT, "\n")
