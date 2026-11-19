#!/usr/bin/env Rscript
# What each candidate sample floor actually yields, so the threshold is chosen against counts
# rather than intuition. Two floors are NOT the same thing here and must not be conflated:
#   TOTAL reads      -- the MRA 2020 criterion (>=1000). 91.6% of a 16S library can be host.
#   BACTERIAL reads  -- the same floor applied after removing chloroplast + mitochondria.
# Q6 needs >=25 genotyped trees to fit; Q5 needs >=20 libraries. Both are marked.

P <- "/home/ecastron/araucaria"
st16 <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx16 <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org <- rep(FALSE, ncol(st16))
for (rk in colnames(tx16)) org <- org | grepl("Chloroplast|Mitochondri", tx16[[rk]], ignore.case = TRUE)
org[is.na(org)] <- FALSE
stI <- readRDS(file.path(P, "microbiome/ITS/seqtab_nochim.rds"))
man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
sink(file.path(P, "results/Q5/threshold_table.txt"), split = TRUE)

rad <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
geno <- unique(rad$tree_id)

blk <- function(cnt, lab, tot) {
  i <- match(rownames(cnt), man$run_accession)
  comp <- man$compartment[i]; tid <- man$tree_id[i]
  d <- rowSums(cnt)
  cat(sprintf("\n%s\n", lab))
  cat(sprintf("  %-9s %-16s", "floor", "metric"))
  for (cc in c("needle","branch","root")) cat(sprintf("%-22s", cc)); cat("\n")
  for (fl in c(1000, 500, 250, 100, 0)) {
    cat(sprintf("  >=%-7d %-16s", fl, if (tot) "TOTAL reads" else "BACTERIAL reads"))
    for (cc in c("needle","branch","root")) {
      s <- which(comp == cc & d >= fl)
      ng <- length(unique(na.omit(tid[s][tid[s] %in% geno])))
      cat(sprintf("%3d libs / %2d trees%-4s", length(s), ng,
                  if (ng >= 25) " OK" else ""))
    }
    cat("\n")
  }
}
blk(st16, "===== 16S, MRA criterion: TOTAL reads, organelles retained =====", TRUE)
blk(st16[, !org, drop = FALSE], "===== 16S, organelles REMOVED: bacterial reads only =====", FALSE)
blk(stI, "===== ITS (uncontaminated) =====", TRUE)
cat("\nQ6 needs >=25 genotyped trees; Q5 needs >=20 libraries.\n")
sink()
