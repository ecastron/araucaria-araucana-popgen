#!/usr/bin/env Rscript
# Step 1: microbiome ASV inference, one amplicon per invocation.
#
# 16S and ITS are processed with DIFFERENT parameters, and must be. ITS1 is length-
# polymorphic across taxa, so truncLen would silently delete every taxon whose amplicon is
# shorter than the cutoff. Measured here: 62% of ITS R1 reads contain the reverse primer
# (read-through), against ~0% for 16S, because ITS1 amplicons are frequently shorter than the
# 251 bp read while the 253 bp V4 amplicon exceeds the 151 bp 16S read.
#
#   16S : truncLen used (fixed-length amplicon), no primer removal needed
#   ITS : NO truncLen; cutadapt strips read-through first; minLen guards short reads
#
# Error models are learned per amplicon, which is the batch structure fixed in the spec:
# one flow cell per amplicon, uniform read length within amplicon, and no plate splitting.
#
# All 80 trees are processed; inference is restricted to the 54 paired trees downstream.

suppressPackageStartupMessages({library(dada2)})

args <- commandArgs(trailingOnly = TRUE)
AMP  <- args[1]                       # "16S" or "ITS"
PROJ <- "/home/ecastron/araucaria"
RAW  <- file.path(PROJ, "raw/microbiome")
OUT  <- file.path(PROJ, "microbiome", AMP)
THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

manifest <- read.csv(file.path(PROJ, "metadata/manifest microbiome runs.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
runs <- manifest$run_accession[manifest$amplicon == AMP]
cat(sprintf("[%s] runs in manifest: %d\n", AMP, length(runs)))

fwd <- file.path(RAW, paste0(runs, "_1.fastq.gz"))
rev <- file.path(RAW, paste0(runs, "_2.fastq.gz"))
keep <- file.exists(fwd) & file.exists(rev)
runs <- runs[keep]; fwd <- fwd[keep]; rev <- rev[keep]
cat(sprintf("[%s] files present: %d\n", AMP, length(runs)))

# ---- primer / read-through handling -----------------------------------------------------
if (AMP == "ITS") {
  # Primers were already stripped from the 5' ends by the submitter, but 62% of reads run
  # through into the reverse primer at the 3' end. cutadapt removes that; without it the
  # primer sequence is treated as biological and inflates ASV counts.
  CUT <- file.path(OUT, "cutadapt"); dir.create(CUT, showWarnings = FALSE)
  cutadapt <- "/opt/miniconda/miniconda-py312_24.3.0-0/envs/cutadapt/bin/cutadapt"
  FWD_RC <- "TTACTTCCTCTAAATGACCAAG"   # revcomp ITS1F
  REV_RC <- "GCATCGATGAAGAACGCAGC"     # revcomp ITS2
  cf <- file.path(CUT, basename(fwd)); cr <- file.path(CUT, basename(rev))
  for (i in seq_along(fwd)) {
    if (file.exists(cf[i]) && file.info(cf[i])$size > 0) next
    system2(cutadapt, args = c("-a", REV_RC, "-A", FWD_RC,
                               "-n", "2", "--minimum-length", "50",
                               "-j", THREADS,
                               "-o", cf[i], "-p", cr[i], fwd[i], rev[i]),
            stdout = FALSE, stderr = FALSE)
  }
  ok <- file.exists(cf) & file.info(cf)$size > 100
  fwd <- cf[ok]; rev <- cr[ok]; runs <- runs[ok]
  cat(sprintf("[ITS] after cutadapt: %d libraries\n", length(runs)))
}

# ---- filter and trim --------------------------------------------------------------------
filtF <- file.path(OUT, "filtered", paste0(runs, "_F.fq.gz"))
filtR <- file.path(OUT, "filtered", paste0(runs, "_R.fq.gz"))

if (AMP == "16S") {
  # Fixed 253 bp amplicon against 2x151 reads. 145+145 leaves ~37 bp overlap for merging.
  out <- filterAndTrim(fwd, filtF, rev, filtR,
                       truncLen = c(145, 145),
                       maxN = 0, maxEE = c(2, 2), truncQ = 2,
                       rm.phix = TRUE, compress = TRUE, multithread = THREADS)
} else {
  # No truncLen: truncating a length-polymorphic marker discards short-amplicon taxa.
  out <- filterAndTrim(fwd, filtF, rev, filtR,
                       truncLen = 0, minLen = 50,
                       maxN = 0, maxEE = c(2, 2), truncQ = 2,
                       rm.phix = TRUE, compress = TRUE, multithread = THREADS)
}
write.csv(out, file.path(OUT, "filter_stats.csv"))

# Libraries that lost everything (the deposit contains ITS runs with 1-3 reads) cannot be
# denoised and are dropped here rather than crashing the run.
alive <- file.exists(filtF) & file.exists(filtR) & file.info(filtF)$size > 100
filtF <- filtF[alive]; filtR <- filtR[alive]; runs <- runs[alive]
cat(sprintf("[%s] libraries surviving filter: %d\n", AMP, length(runs)))

# ---- error model, per amplicon ----------------------------------------------------------
errF <- learnErrors(filtF, multithread = THREADS, randomize = TRUE)
errR <- learnErrors(filtR, multithread = THREADS, randomize = TRUE)
saveRDS(errF, file.path(OUT, "errF.rds")); saveRDS(errR, file.path(OUT, "errR.rds"))

ddF <- dada(filtF, err = errF, multithread = THREADS, pool = FALSE)
ddR <- dada(filtR, err = errR, multithread = THREADS, pool = FALSE)

mrg <- mergePairs(ddF, filtF, ddR, filtR, verbose = FALSE)
seqtab <- makeSequenceTable(mrg)
rownames(seqtab) <- runs
cat(sprintf("[%s] ASVs before chimera removal: %d\n", AMP, ncol(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = THREADS, verbose = TRUE)
cat(sprintf("[%s] ASVs after chimera removal: %d (%.1f%% of reads retained)\n",
            AMP, ncol(seqtab.nochim), 100 * sum(seqtab.nochim) / sum(seqtab)))
saveRDS(seqtab.nochim, file.path(OUT, "seqtab_nochim.rds"))

# ---- taxonomy ---------------------------------------------------------------------------
REFS <- file.path(PROJ, "refs")
if (AMP == "16S") {
  tax <- assignTaxonomy(seqtab.nochim,
                        file.path(REFS, "silva_nr99_v138.2_toSpecies_trainset.fa.gz"),
                        multithread = THREADS, tryRC = FALSE)
} else {
  # tryRC because ITS reads may be in either orientation relative to the reference.
  tax <- assignTaxonomy(seqtab.nochim,
                        file.path(REFS, "sh_general_release_dynamic_19.02.2025.fasta"),
                        multithread = THREADS, tryRC = TRUE)
}
saveRDS(tax, file.path(OUT, "taxonomy.rds"))

# Track reads through every step, so loss is auditable rather than assumed.
getN <- function(x) sum(getUniques(x))
track <- cbind(out[alive, , drop = FALSE],
               denoisedF = sapply(ddF, getN), denoisedR = sapply(ddR, getN),
               merged = sapply(mrg, getN), nonchim = rowSums(seqtab.nochim))
rownames(track) <- runs
write.csv(track, file.path(OUT, "read_tracking.csv"))
cat(sprintf("[%s] DONE. seqtab %d x %d\n", AMP, nrow(seqtab.nochim), ncol(seqtab.nochim)))
