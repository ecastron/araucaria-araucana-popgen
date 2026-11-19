#!/usr/bin/env Rscript
# R1-only versus merged, part 2: what R1-only costs in taxonomy, and whether the two pipelines
# describe the same communities.
#
# Part 1 established, on 20 libraries:
#   retention  R1-only 84.7% vs merged 52.1% of input reads
#   ASVs       R1-only 805    vs merged 567
#   read-through 78.1% of R1 reads already span the whole amplicon
# and found NO pile-up of merged ASVs against the ~490 bp merge ceiling (11 ASVs within 20 bp,
# none beyond), so the merge ceiling is NOT demonstrably truncating the length distribution.
# The case for R1-only therefore rests on read retention and ASV recovery, not on dropout.
#
# Two things still decide it, and they pull in opposite directions:
#   - a 251 bp fragment may resolve to genus less often than the full spacer  (cost)
#   - 33% more reads and 42% more ASVs                                        (benefit)
# The benefit is only real if the extra ASVs are fungi rather than noise, and only matters if
# the two pipelines would lead to different ecological conclusions. Both are tested here.

suppressPackageStartupMessages({library(dada2); library(vegan)})
set.seed(1)

PROJ <- "/home/ecastron/araucaria"
ITS  <- file.path(PROJ, "microbiome/ITS")
CUT  <- file.path(ITS, "cutadapt")
TMP  <- file.path(ITS, "r1_tax"); dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
TH   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
UNITE <- file.path(PROJ, "refs/sh_general_release_dynamic_19.02.2025.fasta")

cf  <- sort(list.files(CUT, pattern = "_1.fastq.gz$", full.names = TRUE))
sel <- sample(seq_along(cf), 40)          # wider than part 1, for a stabler ordination
fwd <- cf[sel]; rev <- sub("_1.fastq.gz$", "_2.fastq.gz", fwd)

den <- function(files, err_files = NULL) {
  e <- learnErrors(files, multithread = TH, randomize = TRUE, verbose = FALSE)
  dada(files, err = e, multithread = TH, pool = FALSE, verbose = FALSE)
}

# ---- R1-only ------------------------------------------------------------------------------
ff <- file.path(TMP, paste0("f_", basename(fwd)))
o1 <- filterAndTrim(fwd, ff, truncLen = 0, minLen = 50, maxN = 0, maxEE = 2,
                    truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
a1 <- file.exists(ff) & file.info(ff)$size > 100
d1 <- den(ff[a1])
t1 <- removeBimeraDenovo(makeSequenceTable(d1), method = "consensus",
                         multithread = TH, verbose = FALSE)
rownames(t1) <- sub("^f_", "", sub("_1.fastq.gz$", "", basename(ff[a1])))

# ---- merged -------------------------------------------------------------------------------
pf <- file.path(TMP, paste0("p_", basename(fwd))); qf <- file.path(TMP, paste0("q_", basename(rev)))
o2 <- filterAndTrim(fwd, pf, rev, qf, truncLen = 0, minLen = 50, maxN = 0, maxEE = c(2, 2),
                    truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
a2 <- file.exists(pf) & file.exists(qf) & file.info(pf)$size > 100
dF <- den(pf[a2]); dR <- den(qf[a2])
mg <- mergePairs(dF, pf[a2], dR, qf[a2], verbose = FALSE)
t2 <- removeBimeraDenovo(makeSequenceTable(mg), method = "consensus",
                         multithread = TH, verbose = FALSE)
rownames(t2) <- sub("^p_", "", sub("_1.fastq.gz$", "", basename(pf[a2])))

cat(sprintf("=== retention (n = %d libraries) ===\nR1-only: %.1f%% of input, %s reads, %d ASVs\nmerged : %.1f%% of input, %s reads, %d ASVs\n",
            length(sel), 100 * sum(t1) / sum(o1[, 1]), format(sum(t1), big.mark = ","), ncol(t1),
            100 * sum(t2) / sum(o2[, 1]), format(sum(t2), big.mark = ","), ncol(t2)))

# ---- taxonomy, one table at a time to stay inside memory ----------------------------------
tax <- function(tt, lab) {
  x <- assignTaxonomy(tt, UNITE, multithread = TH, tryRC = TRUE)
  # Weight by reads as well as by ASV: an unassignable ASV carrying 3 reads matters far less
  # than one carrying 30,000, and the unweighted rate hides that.
  w <- colSums(tt)
  f <- function(r) c(asv = 100 * mean(!is.na(x[, r])),
                     reads = 100 * sum(w[!is.na(x[, r])]) / sum(w))
  for (r in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    v <- f(r); cat(sprintf("  %-8s %-8s %5.1f%% of ASVs | %5.1f%% of reads\n", lab, r, v[1], v[2]))
  }
  # Non-fungal or unassignable-at-phylum ASVs are the noise check on the extra R1 ASVs.
  cat(sprintf("  %-8s ASVs with no phylum: %d of %d\n", lab, sum(is.na(x[, "Phylum"])), nrow(x)))
  invisible(x)
}
cat("\n=== taxonomy ===\n")
x1 <- tax(t1, "R1-only"); rm(x1); gc()
x2 <- tax(t2, "merged");  rm(x2); gc()

# ---- do the two pipelines describe the same communities? ---------------------------------
common <- intersect(rownames(t1), rownames(t2))
# A library can survive filtering yet end with zero reads in one table (the merged pipeline
# empties some libraries the R1 pipeline keeps). Relative abundance is undefined for those, so
# they must be dropped before any distance is computed -- otherwise decostand returns NaN and
# vegdist fails. Dropping them is also the honest comparison: only libraries BOTH pipelines
# retain can be compared, and how many are lost is itself reported.
nz <- common[rowSums(t1[common, , drop = FALSE]) > 0 &
             rowSums(t2[common, , drop = FALSE]) > 0]
cat(sprintf("\n=== community concordance ===\nshared libraries: %d | usable in both (non-empty): %d | dropped: %s\n",
            length(common), length(nz), paste(setdiff(common, nz), collapse = ", ")))
common <- nz
if (length(common) > 5) {
  b1 <- vegdist(decostand(t1[common, , drop = FALSE], "total"), "bray")
  b2 <- vegdist(decostand(t2[common, , drop = FALSE], "total"), "bray")
  m <- mantel(b1, b2, permutations = 999)
  p <- protest(cmdscale(b1, k = 2), cmdscale(b2, k = 2), permutations = 999)
  cat(sprintf("Mantel r = %.3f (p = %.3f)\nProcrustes correlation = %.3f (p = %.3f)\n",
              m$statistic, m$signif, p$scale, p$signif))
  cat(sprintf("observed richness: R1-only median %.0f | merged median %.0f\n",
              median(rowSums(t1[common, ] > 0)), median(rowSums(t2[common, ] > 0))))
  cat("\nINTERPRETATION: a high Mantel/Procrustes correlation means the choice does not change\nthe ecology, so the pipeline retaining more reads should be preferred. A low one means the\nchoice is consequential and must be justified on its own terms, not on read counts.\n")
}
saveRDS(list(r1 = t1, merged = t2), file.path(ITS, "r1_vs_merged_tables.rds"))
unlink(TMP, recursive = TRUE)
