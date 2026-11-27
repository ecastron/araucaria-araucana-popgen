#!/usr/bin/env Rscript
# Revision (25 Nov 2026): sNMF on the 163 analyzed trees, with all SNPs and
# with one SNP per RAD locus. The published Figure 2E used results/Q1/snmf/ours91.geno, built from
# all 177 libraries (relatives included); this rerun matches the 163-tree set of the paper.
# Settings as scripts/run_snmf.R: K = 1 to 6, 10 repetitions, alpha = 100, seed 2026.
# Inputs : results/revision/{all163,thin1}.geno (+ .samples), written by structure.py
# Outputs: results/revision/snmf_<tag>.tsv (cross-entropy by K), snmf_<tag>_K2.Q (best run at K = 2),
#          printed per-stand mean Coastal ancestry
suppressPackageStartupMessages(library(LEA))
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/revision")
pop <- read.delim(file.path(P, "metadata/popmap_177.tsv"), header = FALSE,
                  col.names = c("id", "stand"), stringsAsFactors = FALSE)
coast <- c("VA", "TR", "PN")
for (tag in c("all163", "thin1")) {
  set.seed(2026)
  geno <- file.path(OUT, paste0(tag, ".geno"))
  pr <- snmf(geno, K = 1:6, entropy = TRUE, repetitions = 10, project = "new",
             CPU = 8, alpha = 100)
  ce <- sapply(1:6, function(k) min(cross.entropy(pr, K = k)))
  write.table(data.frame(K = 1:6, min_cross_entropy = ce),
              file.path(OUT, paste0("snmf_", tag, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  q <- Q(pr, K = 2, run = which.min(cross.entropy(pr, K = 2)))
  s <- readLines(paste0(geno, ".samples")); st <- pop$stand[match(s, pop$id)]
  # orient so column 1 is the Coastal cluster
  if (mean(q[st %in% coast, 1]) < mean(q[st %in% coast, 2])) q <- q[, 2:1]
  write.table(data.frame(id = s, stand = st, coastal = q[, 1], andean = q[, 2]),
              file.path(OUT, paste0("snmf_", tag, "_K2.Q")), sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("\n== %s: best K by cross-entropy = %d\n", tag, which.min(ce)))
  print(round(ce, 5))
  m <- tapply(q[, 1], st, mean)
  cat("mean Coastal ancestry per stand:\n"); print(round(sort(m, decreasing = TRUE), 3))
  cat(sprintf("Andean ancestry in Coastal trees: %.3f\n", mean(q[st %in% coast, 2])))
}
