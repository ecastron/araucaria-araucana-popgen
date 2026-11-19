#!/usr/bin/env Rscript
# Core fungal community on the Figure 5 sample set (needle or branch libraries of genotyped trees,
# >= 1,000 reads), with the published package microbiome (Lahti & Shetty) on phyloseq objects.
#   1. detection x prevalence grid with microbiome::core_members / core_abundance on relative abundance
#   2. the same grid after phyloseq::rarefy_even_depth to the shallowest library (equal depth)
#   3. BESPOKE leak check (no package applies: the SRA project has no negative controls or DNA
#      concentrations, so decontam cannot run): for each ASV in >= 50% of libraries, prevalence and
#      read counts by range, and the share of detections carried by 1 or 2 reads.
# Usage: Rscript core_microbiome.R [r1|merged]   (default r1, the pre-registered table)
# Family level via microbiome::aggregate_taxa (unassigned ASVs pooled as one 'Unknown' family, which is
# excluded from the family core below because it is not a taxon).
# Outputs: results/core_<table>/core_grid.tsv, core_members_50.tsv, core_leak_check.tsv, core_summary.txt
suppressPackageStartupMessages({ library(phyloseq); library(microbiome) })
TAB <- {a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "r1"}; stopifnot(TAB %in% c("r1", "merged"))
SFX <- if (TAB == "r1") "_r1" else ""
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, sprintf("results/core_%s", TAB)); dir.create(OUT, showWarnings = FALSE)
man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
st  <- readRDS(file.path(P, sprintf("microbiome/ITS/seqtab_nochim%s.rds", SFX)))
tax <- as.matrix(readRDS(file.path(P, sprintf("microbiome/ITS/taxonomy%s.rds", SFX))))
m <- man[match(rownames(st), man$run_accession), ]
k <- !is.na(m$tree_id) & m$compartment %in% c("needle", "branch") & rowSums(st) >= 1000   # Figure 5 set
x <- st[k, , drop = FALSE]; m <- m[k, ]; x <- x[, colSums(x) > 0, drop = FALSE]
ids <- sprintf("ASV%04d", seq_len(ncol(x))); seqs <- colnames(x); colnames(x) <- ids
tx <- tax[seqs, , drop = FALSE]; rownames(tx) <- ids
sd <- data.frame(range = ifelse(m$cordillera == "Costa", "Coastal", "Andean"), compartment = m$compartment,
                 tree = m$tree_id, row.names = rownames(x))
ps <- phyloseq(otu_table(x, taxa_are_rows = FALSE), tax_table(tx), sample_data(sd))
depth <- min(sample_sums(ps))
sink(file.path(OUT, "core_summary.txt"), split = TRUE)
cat(sprintf("table %s | %d libraries, %d trees, %d ASVs | shallowest library %d reads\n",
            TAB, nsamples(ps), length(unique(sd$tree)), ntaxa(ps), depth))
grid <- expand.grid(detection = c(0, 1e-4, 1e-3, 1e-2), prevalence = c(0.10, 0.25, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90))
run_grid <- function(p, label) {
  rel <- microbiome::transform(p, "compositional")
  if (taxa_are_rows(rel)) rel <- phyloseq::t(rel)   # samples as rows for the abundance sum
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    d <- grid$detection[i]; pr <- grid$prevalence[i]
    il <- d > 0   # detection 0 means present (> 0 reads); otherwise relative abundance >= d
    cm <- setdiff(core_members(rel, detection = d, prevalence = pr, include.lowest = il), "Unknown")
    ca <- if (length(cm)) mean(rowSums(as(otu_table(rel), "matrix")[, cm, drop = FALSE])) else 0
    data.frame(data = label, detection = d, prevalence = pr, n_core = length(cm), mean_core_abundance = ca)
  }))
}
ps_rar <- rarefy_even_depth(ps, sample.size = depth, rngseed = 2026, replace = FALSE, verbose = FALSE)
fam <- function(p) aggregate_taxa(p, level = "Family")   # microbiome; ASVs without a Family are pooled as "Unknown"
G <- rbind(cbind(level = "ASV", rbind(run_grid(ps, "relative abundance"), run_grid(ps_rar, sprintf("rarefied to %d", depth)))),
           cbind(level = "Family", rbind(run_grid(fam(ps), "relative abundance"), run_grid(fam(ps_rar), sprintf("rarefied to %d", depth)))))
write.table(G, file.path(OUT, "core_grid.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(G, row.names = FALSE)
# members at >= 50% prevalence (any detection), with taxonomy
rel <- microbiome::transform(ps, "compositional")
pv <- prevalence(rel, detection = 0, include.lowest = FALSE)     # present = > 0 reads
mem <- names(pv)[pv >= 0.5]
cat(sprintf("\nASVs in >= 50%% of libraries: %d\n", length(mem)))
L <- do.call(rbind, lapply(mem, function(a) {
  cnt <- as.numeric(otu_table(ps)[, a]); r <- as.numeric(otu_table(rel)[, a])
  row <- data.frame(asv = a, genus = tx[a, "Genus"], family = tx[a, "Family"], order = tx[a, "Order"],
                    prevalence = mean(cnt > 0),
                    prev_detect_1e3 = mean(r >= 1e-3), prev_detect_1e2 = mean(r >= 1e-2),
                    prev_rarefied = mean(as.numeric(otu_table(ps_rar)[, a]) > 0))
  for (g in c("Coastal", "Andean")) {
    s <- sd$range == g; c2 <- cnt[s]
    row[[paste0("prev_", g)]] <- mean(c2 > 0)
    row[[paste0("median_reads_present_", g)]] <- if (any(c2 > 0)) median(c2[c2 > 0]) else NA
    row[[paste0("share_detections_le2reads_", g)]] <- if (any(c2 > 0)) mean(c2[c2 > 0] <= 2) else NA
    row[[paste0("median_relab_present_", g)]] <- if (any(c2 > 0)) median(r[s][c2 > 0]) else NA
  }
  row }))
write.table(L, file.path(OUT, "core_leak_check.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(L, row.names = FALSE)
# families in >= 40% of libraries at equal depth and 0.1% detection, with prevalence (Figure S4 text)
fr <- microbiome::transform(fam(ps_rar), "compositional")
fp <- sort(prevalence(fr, detection = 1e-3, include.lowest = TRUE), decreasing = TRUE)
fp <- fp[names(fp) != "Unknown" & fp >= 0.4]
write.table(data.frame(family = names(fp), prevalence_rarefied_0.001 = round(fp, 4)),
            file.path(OUT, "core_families.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nFamilies in >= 40% of libraries (rarefied, >= 0.1% of reads):\n"); print(round(fp, 3))
cat("\nsessionInfo: microbiome", as.character(packageVersion("microbiome")), "phyloseq", as.character(packageVersion("phyloseq")), "\n")
sink()
