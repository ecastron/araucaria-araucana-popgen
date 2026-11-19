#!/usr/bin/env Rscript
# ITS, final production run: FORWARD READS ONLY.
#
# Decided 4 August 2026 on measured evidence, replacing the merged-pair run:
#
#            reads retained   ASVs   median richness   reads assigned to genus
#   R1-only        84.7%       805          30                  83.5%
#   merged         52.1%       567          20                  72.3%
#   agreement between the two: Mantel r = 0.889, Procrustes = 0.958
#
# The two pipelines describe the SAME communities, so the one retaining more data wins -- that
# was the decision rule fixed before the comparison was run. The merged run lost 45% of reads
# entirely to R2 quality: maxEE is an absolute per-read cap, so it bites far harder on an
# untruncated 251 bp ITS read (whose tail is degraded) than on a 16S read truncated to 145 bp.
# ITS cannot be truncated -- it is length-polymorphic -- so the cap could not be rescued.
#
# Dropping R2 also removes the merge ceiling. With 2x251 and minOverlap 12, amplicons beyond
# ~490 bp cannot merge and vanish entirely; R1-only truncates them instead of losing them.
# Measured here: 78% of R1 reads already span the whole amplicon (cutadapt trimmed the reverse
# primer), so merging was contributing nothing for three quarters of the data.
#
# Cost accepted: 251 bp of ITS1 rather than the full spacer, hence slightly lower per-ASV
# assignment at high ranks. Read-weighted assignment is BETTER at every rank below class,
# because the merged pipeline was discarding abundant, well-classified reads.
#
# Supported by Pauvert et al. (2019), who benchmarked 360 pipeline/parameter combinations
# against a 189-strain mock community and concluded that single forward reads with DADA2 is a
# good option for fungal community characterisation.

suppressPackageStartupMessages({library(dada2)})

PROJ    <- "/home/ecastron/araucaria"
AMP     <- "ITS"
OUT     <- file.path(PROJ, "microbiome", AMP)
CUT     <- file.path(OUT, "cutadapt")
THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))

# cutadapt output from the earlier run is reused: primer removal is identical for R1-only and
# is the expensive step. Only R1 is taken forward.
fwd <- sort(list.files(CUT, pattern = "_1.fastq.gz$", full.names = TRUE))
runs <- sub("_1.fastq.gz$", "", basename(fwd))
stopifnot(length(fwd) > 0)
cat(sprintf("[ITS-R1] libraries with cutadapt output: %d\n", length(fwd)))

filt <- file.path(OUT, "filtered_r1", paste0(runs, "_F.fq.gz"))
dir.create(dirname(filt[1]), recursive = TRUE, showWarnings = FALSE)

# multithread = FALSE deliberately: filterAndTrim's parallel dispatch dropped a worker on this
# data when a library filtered to zero reads ("scheduled core did not deliver a result"), which
# aborts the whole call. Single-threaded is slower but cannot fail halfway.
out <- filterAndTrim(fwd, filt,
                     truncLen = 0, minLen = 50,
                     maxN = 0, maxEE = 2, truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread = FALSE)
write.csv(out, file.path(OUT, "filter_stats_r1.csv"))

alive <- file.exists(filt) & file.info(filt)$size > 100
filt <- filt[alive]; runs <- runs[alive]
cat(sprintf("[ITS-R1] libraries surviving filter: %d of %d\n", length(runs), length(alive)))

err <- learnErrors(filt, multithread = THREADS, randomize = TRUE)
saveRDS(err, file.path(OUT, "errF_r1.rds"))

dd <- dada(filt, err = err, multithread = THREADS, pool = FALSE)
seqtab <- makeSequenceTable(dd)
rownames(seqtab) <- runs
cat(sprintf("[ITS-R1] ASVs before chimera removal: %d\n", ncol(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = THREADS, verbose = TRUE)
cat(sprintf("[ITS-R1] ASVs after chimera removal: %d (%.1f%% of reads retained)\n",
            ncol(seqtab.nochim), 100 * sum(seqtab.nochim) / sum(seqtab)))
saveRDS(seqtab.nochim, file.path(OUT, "seqtab_nochim_r1.rds"))

tax <- assignTaxonomy(seqtab.nochim,
                      file.path(PROJ, "refs/sh_general_release_dynamic_19.02.2025.fasta"),
                      multithread = THREADS, tryRC = TRUE)
saveRDS(tax, file.path(OUT, "taxonomy_r1.rds"))

getN <- function(x) sum(getUniques(x))
track <- cbind(out[alive, , drop = FALSE],
               denoised = sapply(dd, getN),
               nonchim  = rowSums(seqtab.nochim))
rownames(track) <- runs
write.csv(track, file.path(OUT, "read_tracking_r1.csv"))

cat(sprintf("[ITS-R1] DONE. seqtab %d x %d | total reads in %s -> out %s (%.1f%%)\n",
            nrow(seqtab.nochim), ncol(seqtab.nochim),
            format(sum(out[, 1]), big.mark = ","),
            format(sum(seqtab.nochim), big.mark = ","),
            100 * sum(seqtab.nochim) / sum(out[, 1])))
cat(sprintf("[ITS-R1] libraries under 1000 reads: %d | under 5000: %d\n",
            sum(rowSums(seqtab.nochim) < 1000), sum(rowSums(seqtab.nochim) < 5000)))
