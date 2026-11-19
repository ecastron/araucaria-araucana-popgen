#!/usr/bin/env Rscript
# sNMF (LEA) -- the SAME algorithm Varas-Myrik et al. (2022) used -- on OUR assemblies.
# Design separates ALGORITHM from DATA: sNMF is run on our working set AND on a
# published-style filtered set; ADMIXTURE is run on the latter for the cross-comparison.
suppressPackageStartupMessages(library(LEA))
args <- commandArgs(TRUE); geno <- args[1]; tag <- args[2]
P <- "/home/ecastron/araucaria"; OUT <- file.path(P,"results/Q1/snmf")
set.seed(2026)
cat(sprintf("\n########## sNMF: %s ##########\n", tag))
pr <- snmf(geno, K = 1:6, entropy = TRUE, repetitions = 10, project = "new",
           CPU = 8, alpha = 100)
ce <- sapply(1:6, function(k) min(cross.entropy(pr, K = k)))
cat(sprintf("%-4s %14s\n","K","min cross-entropy"))
for (k in 1:6) cat(sprintf("%-4d %14.5f %s\n", k, ce[k], ifelse(k==which.min(ce),"<- BEST","")))
best <- which.min(ce)
cat(sprintf("\nsNMF selects K = %d\n", best))
if (best >= 2) {
  q <- Q(pr, K = 2, run = which.min(cross.entropy(pr, K = 2)))
  s <- readLines(paste0(geno,".samples"))
  pop <- read.delim(file.path(P,"metadata/popmap_177.tsv"), header=FALSE,
                    col.names=c("id","stand"), stringsAsFactors=FALSE)
  st <- pop$stand[match(s, pop$id)]
  cat("\nK=2 ancestry, mean of cluster-1 coefficient per stand:\n")
  m <- tapply(q[,1], st, mean)
  for (nm in names(sort(m, decreasing=TRUE))) cat(sprintf("  %-4s %.3f  (n=%d)\n", nm, m[nm], sum(st==nm,na.rm=TRUE)))
}
saveRDS(list(ce=ce, best=best), file.path(OUT, paste0(tag,"_snmf.rds")))
