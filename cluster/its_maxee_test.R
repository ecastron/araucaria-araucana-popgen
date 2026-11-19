#!/usr/bin/env Rscript
# Does relaxing the ITS R2 error cap actually buy usable data?
#
# The filter diagnostic showed maxEE=c(2,2) rejects pairs almost entirely on R2 quality, and
# that relaxing R2 recovers a large fraction. That is not sufficient grounds to change the
# parameter: reads admitted by a looser cap could simply fail later at merging, or inflate the
# ASV count with noise. The decision-relevant comparison is reads surviving to the FINAL table
# and the ASVs they produce, so this runs the full denoise-merge-chimera path under each cap.
#
# Accept the relaxed cap only if final reads rise substantially AND the merge rate holds AND
# the ASV count does not balloon relative to the reads gained.

suppressPackageStartupMessages({library(dada2)})
set.seed(1)

PROJ <- "/home/ecastron/araucaria"
CUT  <- file.path(PROJ, "microbiome/ITS/cutadapt")
TMP  <- file.path(PROJ, "microbiome/ITS/maxee_test")

cf <- sort(list.files(CUT, pattern = "_1.fastq.gz$", full.names = TRUE))
sel <- sample(seq_along(cf), 24)
fwd <- cf[sel]; rev <- sub("_1.fastq.gz$", "_2.fastq.gz", fwd)

for (g in list(c(2, 2), c(2, 5), c(2, 7))) {
  tag <- paste0("ee", g[1], "_", g[2])
  d <- file.path(TMP, tag); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ff <- file.path(d, paste0("f_", basename(fwd))); rf <- file.path(d, paste0("r_", basename(rev)))

  o <- filterAndTrim(fwd, ff, rev, rf, truncLen = 0, minLen = 50, maxN = 0,
                     maxEE = g, truncQ = 2, rm.phix = TRUE, compress = TRUE,
                     multithread = 8, verbose = FALSE)
  alive <- file.exists(ff) & file.exists(rf) & file.info(ff)$size > 100
  ff <- ff[alive]; rf <- rf[alive]

  eF <- learnErrors(ff, multithread = 8, randomize = TRUE, verbose = FALSE)
  eR <- learnErrors(rf, multithread = 8, randomize = TRUE, verbose = FALSE)
  dF <- dada(ff, err = eF, multithread = 8, pool = FALSE, verbose = FALSE)
  dR <- dada(rf, err = eR, multithread = 8, pool = FALSE, verbose = FALSE)
  mg <- mergePairs(dF, ff, dR, rf, verbose = FALSE)
  st <- makeSequenceTable(mg)
  nc <- removeBimeraDenovo(st, method = "consensus", multithread = 8, verbose = FALSE)

  gN <- function(x) sum(getUniques(x))
  merged <- sum(vapply(mg, gN, numeric(1)))
  cat(sprintf(
    "maxEE=c(%d,%s) | filtered %7s (%.1f%% of input) | merged %7s (%.1f%% of filtered) | final %7s | ASVs %4d | libs %d | median lib %5.0f\n",
    g[1], g[2], format(sum(o[, 2]), big.mark = ","), 100 * sum(o[, 2]) / sum(o[, 1]),
    format(merged, big.mark = ","), 100 * merged / sum(o[alive, 2]),
    format(sum(nc), big.mark = ","), ncol(nc), nrow(nc), median(rowSums(nc))))
  unlink(d, recursive = TRUE)
}
