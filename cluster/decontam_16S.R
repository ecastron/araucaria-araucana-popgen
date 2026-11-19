#!/usr/bin/env Rscript
# STEP 1 of the endophyte re-analysis: remove host organellar reads from the 16S table.
#
# WHY THIS IS FIRST AND ALONE. 91.6% of 16S reads are Araucaria chloroplast (75.2%) and
# mitochondria (16.4%). Organellar haplotypes ARE host genotype, so leaving them in a matrix
# used to test "does host genetics predict community composition" can manufacture the answer.
# Nothing downstream is re-run until we know how much real bacterial signal survives, because
# if libraries collapse below usable depth the honest answer may be that 16S cannot support
# the question at all -- and that verdict must be reachable, not designed around.
#
# THRESHOLD DECLARED BEFORE LOOKING: a library needs >=1000 non-organellar reads to enter any
# downstream model, the same floor Q5 already applied. Libraries below it are dropped and
# counted, not rescued.

P <- "/home/ecastron/araucaria"
st <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")),
                    stringsAsFactors = FALSE)
stopifnot(ncol(st) == nrow(tx))
sink(file.path(P, "results/Q5/decontam_16S.txt"), split = TRUE)

cat(sprintf("input: %d libraries x %d ASVs, %s reads\n\n",
            nrow(st), ncol(st), format(sum(st), big.mark = ",")))

# ---- identify organellar ASVs, reporting WHICH rank matched -------------------------------
cat("--- organellar ASVs by taxonomic rank ---\n")
mark <- function(pat) {
  hits <- lapply(colnames(tx), function(rk) which(grepl(pat, tx[[rk]], ignore.case = TRUE)))
  names(hits) <- colnames(tx)
  for (rk in names(hits)) if (length(hits[[rk]]))
    cat(sprintf("  %-11s at %-8s : %5d ASVs\n", pat, rk, length(hits[[rk]])))
  sort(unique(unlist(hits)))
}
cp <- mark("Chloroplast")
mt <- mark("Mitochondri")
org <- sort(unique(c(cp, mt)))
cat(sprintf("\nchloroplast ASVs %d | mitochondrial ASVs %d | union %d\n",
            length(cp), length(mt), length(org)))
cat(sprintf("reads removed: %s (%.2f%%)\n",
            format(sum(st[, org, drop = FALSE]), big.mark = ","),
            100 * sum(st[, org, drop = FALSE]) / sum(st)))

# also drop non-bacterial kingdoms and kingdom-unassigned ASVs
king <- tx[[1]]
bad_k <- which(is.na(king) | !grepl("Bacteria|Archaea", king))
cat(sprintf("non-prokaryote / unassigned kingdom ASVs: %d (%.3f%% of reads)\n",
            length(bad_k), 100 * sum(st[, bad_k, drop = FALSE]) / sum(st)))

drop <- sort(unique(c(org, bad_k)))
clean <- st[, -drop, drop = FALSE]
cat(sprintf("\nsurviving: %d ASVs, %s reads (%.2f%% of input)\n",
            ncol(clean), format(sum(clean), big.mark = ","),
            100 * sum(clean) / sum(st)))

# ---- per-library survival ------------------------------------------------------------------
before <- rowSums(st); after <- rowSums(clean)
cat("\n--- per-library depth AFTER removing host organellar reads ---\n")
cat(sprintf("  median %.0f | mean %.0f | min %.0f | max %.0f\n",
            median(after), mean(after), min(after), max(after)))
qs <- quantile(after, c(.05, .25, .5, .75, .95))
cat(sprintf("  5%%=%.0f  25%%=%.0f  50%%=%.0f  75%%=%.0f  95%%=%.0f\n",
            qs[1], qs[2], qs[3], qs[4], qs[5]))
cat("\n  libraries surviving each depth floor:\n")
for (k in c(100, 500, 1000, 2000, 5000)) {
  cat(sprintf("    >=%5d reads : %3d of %3d (%.1f%%)\n",
              k, sum(after >= k), length(after), 100 * mean(after >= k)))
}
cat(sprintf("\n  host fraction per library: median %.1f%% (range %.1f-%.1f%%)\n",
            100 * median(1 - after / before), 100 * min(1 - after / before),
            100 * max(1 - after / before)))

# ---- by compartment, since that is the analysis unit ---------------------------------------
man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
i <- match(rownames(clean), man$run_accession)
comp <- man$compartment[i]
if (all(is.na(comp))) cat("\n  [manifest key did not match seqtab rownames]\n") else {
  cat("\n--- surviving libraries by compartment (>=1000 non-organellar reads) ---\n")
  cat(sprintf("  %-10s %6s %6s %10s %10s\n", "compartment", "libs", "pass", "med depth", "med host%"))
  for (cc in sort(unique(na.omit(comp)))) {
    s <- which(comp == cc)
    cat(sprintf("  %-10s %6d %6d %10.0f %9.1f%%\n", cc, length(s), sum(after[s] >= 1000),
                median(after[s]), 100 * median(1 - after[s] / before[s])))
  }
}

# ---- prevalence AFTER decontamination ------------------------------------------------------
keep_lib <- after >= 1000
cl2 <- clean[keep_lib, , drop = FALSE]
cl2 <- cl2[, colSums(cl2) > 0, drop = FALSE]
prev <- colSums(cl2 > 0)
rel  <- colSums(cl2) / sum(cl2)
cat(sprintf("\n--- after dropping low-depth libraries: %d libraries x %d ASVs ---\n",
            nrow(cl2), ncol(cl2)))
cat(sprintf("  sparsity: %.1f%%\n", 100 * mean(cl2 == 0)))
cat(sprintf("  %-34s %7s %10s\n", "filter", "ASVs", "reads kept"))
for (k in c(1, 2, 3, 5, 10)) {
  m <- prev >= k
  cat(sprintf("  prevalence >=%2d libraries          %7d %9.2f%%\n",
              k, sum(m), 100 * sum(cl2[, m, drop = FALSE]) / sum(cl2)))
}
for (a in c(1e-5, 1e-4)) {
  m <- rel >= a
  cat(sprintf("  rel. abundance >= %.0e             %7d %9.2f%%\n",
              a, sum(m), 100 * sum(cl2[, m, drop = FALSE]) / sum(cl2)))
}
m5 <- prev >= 5 & rel >= 1e-4
cat(sprintf("  prevalence>=5 AND rel>=1e-4        %7d %9.2f%%   <- proposed\n",
            sum(m5), 100 * sum(cl2[, m5, drop = FALSE]) / sum(cl2)))
cat(sprintf("  sparsity under proposed filter: %.1f%%\n", 100 * mean(cl2[, m5] == 0)))

saveRDS(list(clean = clean, kept_libs = keep_lib, drop_asv = drop,
             prev = prev, rel = rel),
        file.path(P, "microbiome/16S/seqtab_decontam.rds"))
cat(sprintf("\nwritten: %s\n", file.path(P, "microbiome/16S/seqtab_decontam.rds")))
sink()
