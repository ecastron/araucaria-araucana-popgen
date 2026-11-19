#!/usr/bin/env Rscript
# ITS: forward reads only, or merged pairs?
#
# The two options fail in OPPOSITE directions on the same length polymorphism, which is why
# this has to be measured rather than argued:
#
#   merged  : an amplicon too long for the mates to overlap cannot be merged at all, so the
#             taxon DISAPPEARS. With 2x251 and minOverlap = 12 the ceiling is ~490 bp.
#   R1 only : the same amplicon is kept but TRUNCATED at 251 bp, so the taxon survives with
#             less taxonomic resolution.
#
# Dropout is the more damaging of the two for a community study: a truncated sequence still
# counts the taxon as present, an unmergeable one does not. Against that, R1-only loses
# resolution, and loses the error-checking that agreement across the overlap provides.
#
# Four things decide it:
#   1. Is the merged ASV length distribution cut off at the merge ceiling? (dropout evidence)
#   2. How often does R1 alone already span the whole amplicon? (read-through fraction)
#   3. What does R1-only cost in genus-level assignment?
#   4. Do the two give the SAME community answer? If they do, retention wins the argument.

suppressPackageStartupMessages({library(dada2); library(vegan)})
set.seed(1)

PROJ <- "/home/ecastron/araucaria"
ITS  <- file.path(PROJ, "microbiome/ITS")
CUT  <- file.path(ITS, "cutadapt")
TMP  <- file.path(ITS, "r1_test"); dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

# ---- 1. is the merged length distribution truncated at the ceiling? ----------------------
nc <- readRDS(file.path(ITS, "seqtab_nochim.rds"))
L  <- nchar(colnames(nc))
cat("=== merged ASV lengths ===\n")
print(quantile(L, c(0, .5, .9, .95, .99, 1)))
cat(sprintf("ASVs within 20 bp of the ~490 merge ceiling: %d (%.2f%%)\n",
            sum(L >= 470), 100 * mean(L >= 470)))
cat(sprintf("ASVs 400-460: %d | 460-490: %d | >490: %d\n",
            sum(L >= 400 & L < 460), sum(L >= 460 & L <= 490), sum(L > 490)))
# A distribution ending naturally tapers; one ending at a hard ceiling stacks against it.

# ---- 2. how often does R1 alone already cover the amplicon? ------------------------------
# cutadapt trimmed the reverse primer off R1 wherever the read ran through into it. A trimmed
# R1 therefore means the full amplicon is contained in R1 and nothing is gained by merging.
cf <- sort(list.files(CUT, pattern = "_1.fastq.gz$", full.names = TRUE))
sel <- sample(seq_along(cf), 20)
tot <- 0; through <- 0
for (f in cf[sel]) {
  fq <- ShortRead::readFastq(f)
  if (!length(fq)) next
  w <- Biostrings::width(fq)
  tot <- tot + length(w); through <- through + sum(w < 251)
}
cat(sprintf("\n=== read-through ===\nR1 reads shorter than 251 bp (full amplicon in R1): %.1f%% of %s\n",
            100 * through / tot, format(tot, big.mark = ",")))

# ---- 3 & 4. run R1-only and compare against the merged table -----------------------------
fwd <- cf[sel]; rev <- sub("_1.fastq.gz$", "_2.fastq.gz", fwd)
ff  <- file.path(TMP, paste0("f_", basename(fwd)))
rf  <- file.path(TMP, paste0("r_", basename(rev)))

# R1-only: no R2 cap to satisfy, which is where the 45% loss came from.
# filterAndTrim is run single-threaded throughout: its parallel dispatch dropped a worker on
# this data ("scheduled core did not deliver a result") when a library filtered down to zero
# reads, which kills the whole call. Slower, but it cannot fail halfway.
o1 <- filterAndTrim(fwd, ff, truncLen = 0, minLen = 50, maxN = 0, maxEE = 2,
                    truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# Paired, exactly as the production run did it, on the same libraries for a fair contrast.
o2 <- filterAndTrim(fwd, file.path(TMP, paste0("p_", basename(fwd))),
                    rev, file.path(TMP, paste0("q_", basename(rev))),
                    truncLen = 0, minLen = 50, maxN = 0, maxEE = c(2, 2),
                    truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)

a1 <- file.exists(ff) & file.info(ff)$size > 100
e1 <- learnErrors(ff[a1], multithread = 8, randomize = TRUE, verbose = FALSE)
d1 <- dada(ff[a1], err = e1, multithread = 8, pool = FALSE, verbose = FALSE)
t1 <- removeBimeraDenovo(makeSequenceTable(d1), method = "consensus", multithread = 8, verbose = FALSE)

pf <- file.path(TMP, paste0("p_", basename(fwd))); qf <- file.path(TMP, paste0("q_", basename(rev)))
a2 <- file.exists(pf) & file.exists(qf) & file.info(pf)$size > 100
eF <- learnErrors(pf[a2], multithread = 8, randomize = TRUE, verbose = FALSE)
eR <- learnErrors(qf[a2], multithread = 8, randomize = TRUE, verbose = FALSE)
dF <- dada(pf[a2], err = eF, multithread = 8, pool = FALSE, verbose = FALSE)
dR <- dada(qf[a2], err = eR, multithread = 8, pool = FALSE, verbose = FALSE)
mg <- mergePairs(dF, pf[a2], dR, qf[a2], verbose = FALSE)
t2 <- removeBimeraDenovo(makeSequenceTable(mg), method = "consensus", multithread = 8, verbose = FALSE)

cat(sprintf("\n=== retention ===\nR1-only : %6.1f%% of input reads reach the final table (%s reads, %d ASVs)\n",
            100 * sum(t1) / sum(o1[, 1]), format(sum(t1), big.mark = ","), ncol(t1)))
cat(sprintf("merged  : %6.1f%% of input reads reach the final table (%s reads, %d ASVs)\n",
            100 * sum(t2) / sum(o2[, 1]), format(sum(t2), big.mark = ","), ncol(t2)))
cat(sprintf("libraries under 1000 reads -- R1-only: %d | merged: %d\n",
            sum(rowSums(t1) < 1000), sum(rowSums(t2) < 1000)))

# Taxonomy: does the shorter sequence cost genus-level assignment?
UNITE <- file.path(PROJ, "refs/sh_general_release_dynamic_19.02.2025.fasta")
gr <- function(tt) 100 * mean(!is.na(tt[, "Genus"]))
x1 <- assignTaxonomy(t1, UNITE, multithread = 8, tryRC = TRUE)
x2 <- assignTaxonomy(t2, UNITE, multithread = 8, tryRC = TRUE)
cat(sprintf("\n=== taxonomy (%% ASVs assigned to genus) ===\nR1-only: %.1f%% | merged: %.1f%%\n",
            gr(x1), gr(x2)))
cat(sprintf("phylum-level: R1-only %.1f%% | merged %.1f%%\n",
            100 * mean(!is.na(x1[, "Phylum"])), 100 * mean(!is.na(x2[, "Phylum"]))))

# The decisive test: same libraries, same ecology? If the ordinations agree, the extra reads
# are free and R1-only wins; if they disagree, the choice changes conclusions and needs care.
common <- intersect(rownames(t1), rownames(t2))
if (length(common) > 4) {
  b1 <- vegdist(decostand(t1[common, , drop = FALSE], "total"), "bray")
  b2 <- vegdist(decostand(t2[common, , drop = FALSE], "total"), "bray")
  m  <- mantel(b1, b2, permutations = 999)
  p  <- protest(cmdscale(b1, k = 2), cmdscale(b2, k = 2), permutations = 999)
  cat(sprintf("\n=== do they agree? (n = %d libraries) ===\nMantel r = %.3f (p = %.3f)\nProcrustes correlation = %.3f (p = %.3f)\n",
              length(common), m$statistic, m$signif, p$scale, p$signif))
}
unlink(TMP, recursive = TRUE)
