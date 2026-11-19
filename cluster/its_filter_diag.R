#!/usr/bin/env Rscript
# Why does the ITS filter drop 45% of reads?
#
# Suspicion: maxEE is an ABSOLUTE cap on expected errors per read, so it gets stricter as reads
# get longer. 16S reads are truncated to 145 bp before the cap is applied; ITS reads cannot be
# truncated (length-polymorphic marker) so the cap acts over the full ~251 bp, where the
# low-quality tail lives. If so the loss is a side effect of the no-truncation decision, not a
# property of the data, and the fix is to relax the cap for ITS rather than to truncate.
#
# Tests, on a random sample of libraries: which mate causes the rejection, and what a relaxed
# cap recovers. Nothing here changes the run; it only measures.

suppressPackageStartupMessages({library(dada2)})
set.seed(1)

PROJ <- "/home/ecastron/araucaria"
CUT  <- file.path(PROJ, "microbiome/ITS/cutadapt")
TMP  <- file.path(PROJ, "microbiome/ITS/filter_diag"); dir.create(TMP, showWarnings = FALSE)

cf <- sort(list.files(CUT, pattern = "_1.fastq.gz$", full.names = TRUE))
sel <- sample(seq_along(cf), 12)
fwd <- cf[sel]; rev <- sub("_1.fastq.gz$", "_2.fastq.gz", fwd)

# Expected-error distribution per mate, so the cause is measured rather than inferred.
ee <- function(f, n = 20000) {
  fq <- ShortRead::readFastq(f)
  if (length(fq) == 0) return(c(NA, NA, NA))
  q  <- as(Biostrings::quality(fq), "matrix")
  e  <- rowSums(10^(-q / 10), na.rm = TRUE)
  c(median(e, na.rm = TRUE), mean(e > 2, na.rm = TRUE), median(Biostrings::width(fq)))
}
cat("per-mate expected errors (median EE, frac EE>2, median length)\n")
for (i in seq_along(fwd)) {
  a <- ee(fwd[i]); b <- ee(rev[i])
  cat(sprintf("  %-14s R1 EE=%5.2f frac>2=%.2f len=%3.0f | R2 EE=%6.2f frac>2=%.2f len=%3.0f\n",
              sub("_1.fastq.gz", "", basename(fwd[i])), a[1], a[2], a[3], b[1], b[2], b[3]))
}

# What each candidate cap retains. maxEE is per mate; a pair dies if either mate fails.
grid <- list(c(2, 2), c(2, 5), c(3, 5), c(5, 5), c(2, Inf))
cat("\nretention by maxEE (pairs kept / pairs in):\n")
for (g in grid) {
  o <- filterAndTrim(fwd, file.path(TMP, paste0("f_", basename(fwd))),
                     rev, file.path(TMP, paste0("r_", basename(rev))),
                     truncLen = 0, minLen = 50, maxN = 0, maxEE = g, truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread = 8, verbose = FALSE)
  cat(sprintf("  maxEE=c(%s) -> %.1f%% retained (%s / %s)\n",
              paste(g, collapse = ","), 100 * sum(o[, 2]) / sum(o[, 1]),
              format(sum(o[, 2]), big.mark = ","), format(sum(o[, 1]), big.mark = ",")))
  unlink(list.files(TMP, full.names = TRUE))
}
