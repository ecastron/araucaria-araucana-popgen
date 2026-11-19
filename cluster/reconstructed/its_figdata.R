#!/usr/bin/env Rscript
# Reconstructs the ITS descriptive figure tables in figdata/:
#   its_prevalence.csv   per-ASV prevalence (Figure S4; 2,645 ASVs, 164 samples)
#   its_richness.csv     observed ASV richness per sample (old S7 panel C)
#   its_composition.csv  class-level composition per stand (old S7 panel B)
#   its_ordination.csv   CLR PCA scores per sample (old S7 panel A)
#   its_stand_pairs.tsv  stand-centroid CLR distance vs geographic distance (36 pairs)
#
# Reconstructed 18-19 Nov 2026 to replace lost inline code (originals written 26 Aug 2026).
# Outputs go to results/reconstructed/ and are compared with results/figdata/ by
# compare_its_figdata.R; nothing in results/figdata/ is touched.
#
# INPUTS
#   microbiome/ITS/seqtab_nochim.rds   DADA2 MERGED-PAIR table (dada2_step1.R, 282 x 4,061).
#                                      NOT seqtab_nochim_r1.rds: only the merged table
#                                      reproduces the 164 samples / 2,645 ASVs of the originals.
#   microbiome/ITS/taxonomy.rds        UNITE taxonomy for that table
#   metadata/manifest microbiome runs.csv, metadata/manifest rad runs.csv,
#   metadata/chelsa_per_tree.tsv       (coordinates for the stand-pair table)
#
# SAMPLE SET (prevalence, richness, composition, ordination): library level; needle or
# branch; tree_id present in the RAD manifest (that is where the cordillera label comes
# from); >= 250 reads. Gives 164 libraries from 53 trees in 9 stands, 2,645 non-zero ASVs.
# This is the "164 with a cordillera assignment / 53 trees" subset of ANALYSIS_SPEC.md.
#
# STAND-PAIR SET: the logic of q9_microbiome_ibd.R (ITS block), which wrote
# results/Q5/microbiome_ibd.txt: counts pooled per tree over needle and branch, trees with
# >= 250 reads (79), ASVs in >= 5 trees (333; q9 itself also applied rel >= 1e-4), CLR, stand
# centroids, Euclidean distance; haversine km between stand-mean tree coordinates;
# between = 1 when one stand is west of -72.5 (Coastal) and the other east of it.
#
# VERIFICATION (18 Nov 2026, R 4.4.0, vegan 2.6.6.1; compare_its_figdata.R, log in
# results/reconstructed/its_compare.txt):
#   its_prevalence.csv, its_richness.csv, its_composition.csv: byte-identical (md5).
#   its_ordination.csv: all labels identical; PC1, PC2 max |diff| 1.1e-13 (last printed
#     digit only, floating-point summation order). PC1 10.14%, PC2 6.41% (old S7 axis labels
#     10.1%, 6.4%).
#   its_stand_pairs.tsv: labels and `between` identical; km max |diff| 1.0e-12, dissim max
#     |diff| 1.0e-13. r(ln km, dissim) = +0.138 (the log's +0.136 is the q9 run with the
#     extra abundance filter, not this file).
#
# ITS TABLE SELECTION (added 19 Nov 2026). "merged" (default) reproduces everything above
# unchanged. "r1" reads the forward-read table (dada2_its_r1.R), the pipeline chosen on
# 4 Aug 2026 (ANALYSIS_SPEC.md). Choose with the first command-line argument or ITS_TABLE;
# OUT_DIR redirects output. The 164 / 2,645 / 53 and 248 assertions describe the merged table
# only, so they are checked only when ITS_TABLE is "merged".
# The final section (also added 19 Nov 2026) recomputes Figure S4 on the Figure 5 sample set
# (fungal_fig.R: needle or branch, tree_id present, >= 1,000 reads) so that S4 and Figure 5
# describe one data set, and reruns the core sweep (counts only) on both sample definitions.
suppressPackageStartupMessages({library(vegan)})
P   <- "/home/ecastron/araucaria"
ITS_TABLE <- {a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else Sys.getenv("ITS_TABLE", "merged")}
stopifnot(ITS_TABLE %in% c("merged", "r1"))
SFX <- if (ITS_TABLE == "r1") "_r1" else ""
cat(sprintf("ITS table: %s (seqtab_nochim%s.rds, taxonomy%s.rds)\n", ITS_TABLE, SFX, SFX))
OUT <- Sys.getenv("OUT_DIR", file.path(P, "results/reconstructed")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

man    <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
radman <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
ch     <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)
st     <- readRDS(file.path(P, sprintf("microbiome/ITS/seqtab_nochim%s.rds", SFX)))
tax    <- as.data.frame(readRDS(file.path(P, sprintf("microbiome/ITS/taxonomy%s.rds", SFX))), stringsAsFactors = FALSE)
cord_of <- setNames(radman$cordillera, radman$tree_id)

# Same zero replacement as every other CLR in this project (q5/q6/q9 scripts).
clr <- function(m) { m <- as.matrix(m); if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

# ---------------- library-level sample set ----------------
m <- man[match(rownames(st), man$run_accession), ]
k <- !is.na(m$tree_id) & m$compartment %in% c("needle", "branch") &
     !is.na(cord_of[m$tree_id]) & rowSums(st) >= 250
x <- st[k, , drop = FALSE]; m <- m[k, ]
x <- x[, colSums(x) > 0, drop = FALSE]
m$range <- ifelse(cord_of[m$tree_id] == "Costa", "Coastal", "Andean")
if (ITS_TABLE == "merged") stopifnot(nrow(x) == 164, ncol(x) == 2645, length(unique(m$tree_id)) == 53)
cat(sprintf("old S4 set (RAD-manifest tree, >= 250 reads): %d libraries, %d trees, %d ASVs\n",
    nrow(x), length(unique(m$tree_id)), ncol(x)))

# prevalence: share of the 164 samples in which an ASV has >= 1 read (seqtab column order)
write.csv(data.frame(prevalence = unname(colMeans(x > 0))),
          file.path(OUT, "its_prevalence.csv"), row.names = FALSE)

# richness: ASVs observed per sample, no rarefaction
write.csv(data.frame(richness = unname(rowSums(x > 0)), stand = m$site_code,
                     cordillera = m$range, compartment = m$compartment),
          file.path(OUT, "its_richness.csv"), row.names = FALSE)

# composition: per-sample relative abundance summed to UNITE class, then averaged over the
# samples of each stand. The eight classes kept are ranked by their mean over all samples.
cls <- tax[colnames(x), "Class"]; cls[is.na(cls)] <- "unassigned"
A   <- rowsum(t(x / rowSums(x)), cls)                       # class x sample
byS <- sapply(split(seq_len(nrow(x)), m$site_code), function(i) rowMeans(A[, i, drop = FALSE]))
top <- names(sort(rowMeans(A), decreasing = TRUE))[1:8]
write.csv(round(byS[top, , drop = FALSE], 5), file.path(OUT, "its_composition.csv"))

# ordination: PCA (prcomp, centered, unscaled) of CLR. ASVs are those present in >= 3 of
# the 248 above-ground libraries with >= 250 reads (the ag_250 library set, with or without a
# RAD-manifest tree), 637 ASVs, kept even where they are absent from all 164 plotted samples;
# CLR and PCA are then computed on the 164 plotted samples. This filter is recovered from the
# output (it is the only one tried that reproduces it), not from surviving code.
lib248 <- st[!is.na(man$tree_id[match(rownames(st), man$run_accession)]) &
             man$compartment[match(rownames(st), man$run_accession)] %in% c("needle", "branch") &
             rowSums(st) >= 250, , drop = FALSE]
if (ITS_TABLE == "merged") stopifnot(nrow(lib248) == 248)
keep <- colnames(lib248)[colSums(lib248 > 0) >= 3]
X  <- clr(st[rownames(x), keep, drop = FALSE])
pc <- prcomp(X)
write.csv(data.frame(PC1 = unname(pc$x[, 1]), PC2 = unname(pc$x[, 2]), stand = m$site_code,
                     cordillera = m$range, compartment = m$compartment, tree = m$tree_id),
          file.path(OUT, "its_ordination.csv"), row.names = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
cat(sprintf("ordination: %d ASVs, PC1 %.2f%%, PC2 %.2f%%\n", ncol(X), ve[1], ve[2]))

# ---------------- stand-pair table (q9_microbiome_ibd.R logic) ----------------
hav <- function(la1, lo1, la2, lo2) { R <- 6371000; p1 <- la1*pi/180; p2 <- la2*pi/180
  2*R*asin(sqrt(sin((p2-p1)/2)^2 + cos(p1)*cos(p2)*sin((lo2-lo1)*pi/180/2)^2)) }
i   <- match(rownames(st), man$run_accession)
ag  <- man$compartment[i] %in% c("branch", "needle") & !is.na(man$tree_id[i])
s2  <- st[ag, , drop = FALSE]; tid <- man$tree_id[i][ag]; tr <- sort(unique(tid))
cnt <- t(sapply(tr, function(z) colSums(s2[tid == z, , drop = FALSE])))
sel <- rowSums(cnt) >= 250; cnt <- cnt[sel, , drop = FALSE]; tr <- tr[sel]
# NOTE: prevalence >= 5 trees ONLY. q9 also required relative abundance >= 1e-4 (302 ASVs);
# the surviving its_stand_pairs.tsv is reproduced only without that second condition (333 ASVs).
cnt <- cnt[, colSums(cnt > 0) >= 5, drop = FALSE]
cat(sprintf("stand pairs: %d trees, %d ASVs\n", nrow(cnt), ncol(cnt)))
mg <- man[match(tr, man$tree_id), ]; j <- match(tr, ch$tree_id)
lat <- ch$lat[j]; lon <- ch$lon[j]
stand <- droplevels(factor(mg$site_code))
cord <- ifelse(lon < -72.5, "Costa", "Andes")       # as in q9: the Central Valley gap
ns <- tapply(cord, stand, function(z) z[1])
Xt <- clr(cnt)
cent <- t(sapply(levels(stand), function(s) colMeans(Xt[stand == s, , drop = FALSE])))
Mc <- as.matrix(dist(cent))
slat <- tapply(lat, stand, mean); slon <- tapply(lon, stand, mean)
pr <- t(combn(levels(stand), 2))
sp <- data.frame(stand_a = pr[, 1], stand_b = pr[, 2],
                 km = hav(slat[pr[, 1]], slon[pr[, 1]], slat[pr[, 2]], slon[pr[, 2]]) / 1000,
                 dissim = Mc[pr], between = as.integer(ns[pr[, 1]] != ns[pr[, 2]]))
write.table(sp, file.path(OUT, "its_stand_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("stand pairs: r(ln km, dissim) = %+.3f over %d pairs\n", cor(log(sp$km), sp$dissim), nrow(sp)))

# ---------------- Figure S4 on the Figure 5 sample set (added 19 Nov 2026) ----------------
# Panel A: per-ASV prevalence. Panel B: core sweep, the number of taxa reaching >= detection
# relative abundance in >= prevalence of samples, at ASV level and at the rank-propagated Family
# level of q6_family2.R (ASVs with no Family inherit their finest assigned rank, so no reads are
# dropped). Counts only: the permutation tests of q6_core_sweep.R / q6_family2.R are not part of
# Figure S4. Three sample definitions are reported so old and new numbers sit side by side:
#   fig5   needle or branch, tree_id present, >= 1,000 reads (fungal_fig.R; Figure 5)
#   old164 the panel A set above (RAD-manifest tree, >= 250 reads)
#   paired the panel B set of q6_core_sweep.R (Q5 manifests, ngsrelate trees, both
#          compartments >= 250 reads)
fam  <- tax[colnames(st), "Family"]; lab <- as.character(fam)
for (r in c("Order", "Class", "Phylum", "Kingdom")) {
  need <- is.na(lab); if (!any(need)) break
  src <- as.character(tax[colnames(st), r]); ok <- need & !is.na(src)
  lab[ok] <- paste0(tolower(substr(r, 1, 1)), "__", src[ok], "_unclassified")
}
lab[is.na(lab)] <- "k__Fungi_unclassified"
names(lab) <- colnames(st)   # indexed by name: the old164 set drops all-zero columns

m5 <- man[match(rownames(st), man$run_accession), ]
k5 <- !is.na(m5$tree_id) & m5$compartment %in% c("needle", "branch") & rowSums(st) >= 1000
x5 <- st[k5, , drop = FALSE]

# paired set, selection copied from q6_core_sweep.R lines 8-39
rel <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
tree_of <- setNames(radman$tree_id, radman$run_accession)
ta <- tree_of[bare(as.character(rel$ida))]; tb <- tree_of[bare(as.character(rel$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta != tb & rel$nSites >= 2000 & is.finite(rel$theta)
ktrees <- sort(unique(c(ta[gd], tb[gd])))
parts <- list()
for (comp in c("needle", "branch")) {
  o <- readRDS(file.path(P, sprintf("results/Q5/q5_ITS_%s_objects.rds", comp)))
  mm <- o$manifest; mm <- mm[!is.na(mm$tree_id), ]
  ids <- intersect(rownames(st), mm$run_accession); mm <- mm[match(ids, mm$run_accession), ]
  xx <- st[ids, , drop = FALSE]; rownames(xx) <- mm$tree_id
  parts[[comp]] <- xx[rownames(xx) %in% ktrees, , drop = FALSE]
}
sh <- intersect(rownames(parts$needle), rownames(parts$branch))
xp <- rbind(parts$needle[sh, , drop = FALSE], parts$branch[sh, , drop = FALSE])
tp <- rep(sh, 2); kp <- rowSums(xp) >= 250; xp <- xp[kp, , drop = FALSE]; tp <- tp[kp]
okp <- names(which(table(tp) == 2)); xp <- xp[tp %in% okp, , drop = FALSE]

sets <- list(fig5 = x5, old164 = x, paired = xp)
sumrows <- list(); core <- list()
for (nm in names(sets)) {
  X <- sets[[nm]]; Xa <- X[, colSums(X) > 0, drop = FALSE]
  pv <- colMeans(Xa > 0)
  sumrows[[nm]] <- data.frame(set = nm, table = ITS_TABLE, n_samples = nrow(X),
      n_asv = ncol(Xa), pct_single_sample = 100 * mean(colSums(Xa > 0) == 1),
      max_prevalence = 100 * max(pv), n_asv_ge50pct = sum(pv >= 0.5),
      median_reads = median(rowSums(X)))
  if (nm == "fig5") write.csv(data.frame(prevalence = unname(pv)),
                              file.path(OUT, "its_prevalence_fig5set.csv"), row.names = FALSE)
  for (lev in c("ASV", "Family")) {
    M <- if (lev == "ASV") X else t(rowsum(t(X), lab[colnames(X)]))
    R <- M / rowSums(M)
    for (det in c(1e-4, 1e-3)) for (prv in c(0.25, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)) {
      inc <- colSums(R >= det) >= prv * nrow(R)
      core[[length(core) + 1]] <- data.frame(set = nm, table = ITS_TABLE, level = lev,
          detection = det, prevalence = prv, n_taxa = sum(inc),
          pct_reads = 100 * sum(M[, inc, drop = FALSE]) / sum(M))
    }
  }
}
S4 <- do.call(rbind, sumrows); CORE <- do.call(rbind, core)
print(S4, row.names = FALSE)
print(CORE[CORE$n_taxa > 0, ], row.names = FALSE)
write.csv(S4, file.path(OUT, "its_s4_summary.csv"), row.names = FALSE)
write.csv(CORE, file.path(OUT, "its_core_sweep.csv"), row.names = FALSE)
cat("DONE\n")
