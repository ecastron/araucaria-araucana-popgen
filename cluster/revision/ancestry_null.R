#!/usr/bin/env Rscript
# Revision (26 Nov 2026): sNMF K = 2 on each simulated replicate written by
# ancestry_null.py, with the Figure 2E settings (LEA::snmf, 10 repetitions, alpha = 100,
# best run by cross-entropy, seed 2026), plus the observed 163-tree matrix (all SNPs and an
# evenly spaced 20,000-SNP subset).
# Simulated individuals are grouped by stand in observed2.pkl stand order; the stand labels
# below are read from results/revision/ancnull/stand_order.tsv (written by the sbatch wrapper).
# Per dataset: mean Andean component in Coastal trees, mean Coastal component per Coastal
# stand, mean Coastal component in Andean trees and its range over Andean stands.
# Usage: ancestry_null.R TASK NTASKS
# Output: results/revision/ancnull/snmf_task<TASK>.tsv
suppressPackageStartupMessages(library(LEA))
a <- commandArgs(TRUE); task <- as.integer(a[1]); ntask <- as.integer(a[2])
P <- "/home/ecastron/araucaria"; D <- file.path(P, "results/revision/ancnull")
so <- read.delim(file.path(D, "stand_order.tsv"))           # stand, n, coastal (grouped order)
sim_st <- rep(so$stand, so$n)
coast <- so$stand[so$coastal]
pop <- read.delim(file.path(P, "metadata/popmap_177.tsv"), header = FALSE, col.names = c("id", "stand"))
obs_st <- pop$stand[match(readLines(file.path(P, "results/revision/all163.geno.samples")), pop$id)]

summ <- function(q, st) {
  if (mean(q[st %in% coast, 1]) < mean(q[st %in% coast, 2])) q <- q[, 2:1]   # column 1 = Coastal
  ma <- tapply(q[!st %in% coast, 1], st[!st %in% coast], mean)
  c(andean_in_coastal = mean(q[st %in% coast, 2]),
    setNames(sapply(coast, function(s) mean(q[st == s, 1])), paste0("coastal_", coast)),
    coastal_in_andean = mean(q[!st %in% coast, 1]),
    andean_stand_min = min(ma), andean_stand_max = max(ma))
}
files <- sort(list.files(D, pattern = "^H[23]_[0-9]+\\.geno$"))
rk <- as.integer(sub("^H[23]_([0-9]+)\\.geno$", "\\1", files))
files <- files[rk %% ntask == task - 1]   # same split as the Python step
if (task == 1) files <- c(files, "obs20k.geno")
res <- list()
for (f in files) {
  set.seed(2026)
  pr <- snmf(file.path(D, f), K = 2, entropy = TRUE, repetitions = 10, project = "new",
             CPU = 4, alpha = 100)
  q <- Q(pr, K = 2, run = which.min(cross.entropy(pr, K = 2)))
  st <- if (f == "obs20k.geno") obs_st else sim_st
  stopifnot(nrow(q) == length(st))
  res[[f]] <- c(dataset = sub("\\.geno$", "", f), round(summ(q, st), 5))
  cat(f, res[[f]], "\n")
}
if (task == 1) {   # observed, all SNPs, from the existing 163-tree K = 2 run (snmf.R)
  o <- read.delim(file.path(P, "results/revision/snmf_all163_K2.Q"))
  res[["obs_all"]] <- c(dataset = "obs_all", round(summ(as.matrix(o[, c("coastal", "andean")]), o$stand), 5))
}
write.table(do.call(rbind, res), file.path(D, sprintf("snmf_task%02d.tsv", task)),
            sep = "\t", row.names = FALSE, quote = FALSE)
