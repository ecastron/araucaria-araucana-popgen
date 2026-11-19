#!/usr/bin/env Rscript
# Distribution of 16S reads per library AFTER removing chloroplast + mitochondria ONLY.
# No prevalence filter, no abundance filter, no sample floor. All 296 libraries retained.

P <- "/home/ecastron/araucaria"
st <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org <- rep(FALSE, ncol(st))
for (rk in colnames(tx)) org <- org | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case = TRUE)
org[is.na(org)] <- FALSE
man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
sink(file.path(P, "results/Q5/depth_distribution.txt"), split = TRUE)

tot <- rowSums(st)
bac <- rowSums(st[, !org, drop = FALSE])
comp <- man$compartment[match(rownames(st), man$run_accession)]

qs <- c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1)
show <- function(x, lab) {
  cat(sprintf("\n%s  (n = %d)\n", lab, length(x)))
  q <- quantile(x, qs)
  cat("  "); for (i in seq_along(qs)) cat(sprintf("%5.0f%%", 100*qs[i])); cat("\n  ")
  for (i in seq_along(qs)) cat(sprintf("%6s", format(round(q[i]), big.mark=","))); cat("\n")
  cat(sprintf("  mean %.0f | sd %.0f | zero-read libraries %d\n", mean(x), sd(x), sum(x == 0)))
}

cat("################ 16S library depth, organelles REMOVED, no other filter ################\n")
show(bac, "ALL COMPARTMENTS")
for (cc in sort(unique(na.omit(comp)))) show(bac[which(comp == cc)], toupper(cc))

cat("\n\n################ histogram (log10 bins) ################\n")
brk <- c(-1, 0, 1, 10, 30, 100, 300, 1000, 3000, 10000, 30000, Inf)
lab <- c("0", "1", "2-10", "11-30", "31-100", "101-300", "301-1k",
         "1k-3k", "3k-10k", "10k-30k", ">30k")
cut_ <- cut(bac, brk, labels = lab)
cat(sprintf("  %-10s %6s %7s | %s\n", "reads", "n", "cum%", "needle / branch / root"))
cum <- 0
for (i in seq_along(lab)) {
  n <- sum(cut_ == lab[i], na.rm = TRUE); cum <- cum + n
  pc <- sapply(c("needle","branch","root"),
               function(cc) sum(cut_ == lab[i] & comp == cc, na.rm = TRUE))
  cat(sprintf("  %-10s %6d %6.1f%% | %4d / %4d / %4d %s\n", lab[i], n, 100*cum/length(bac),
              pc[1], pc[2], pc[3], strrep("#", round(n/3))))
}

cat("\n\n################ host fraction per library ################\n")
hf <- 1 - bac/tot
show(round(100*hf, 2), "HOST % OF READS, ALL COMPARTMENTS")
for (cc in sort(unique(na.omit(comp)))) show(round(100*hf[which(comp == cc)], 2), paste("HOST %", toupper(cc)))

cat("\n\n################ libraries retained at each floor ################\n")
cat(sprintf("  %-14s %8s %9s %9s %8s\n", "floor", "all", "needle", "branch", "root"))
for (fl in c(0, 1, 10, 30, 100, 200, 300, 500, 1000, 2000, 5000)) {
  cat(sprintf("  >=%-12d %8d %9d %9d %8d\n", fl, sum(bac >= fl),
      sum(bac >= fl & comp == "needle", na.rm=TRUE),
      sum(bac >= fl & comp == "branch", na.rm=TRUE),
      sum(bac >= fl & comp == "root",   na.rm=TRUE)))
}
sink()
