#!/usr/bin/env Rscript
# Figure S2E, polyRAD route: every polyRAD-route row and its own unfiltered baseline, same loci
# universe (34,674 polyRAD loci), same 7 depth bins, same F estimator.
#
# PRODUCES  results/derived/s2e/polyrad_route.tsv     one row per intervention
#           results/derived/s2e/locus_filter_F.tsv    figdata/locus_filter_F.tsv
#           results/derived/s2e/ipyrad_binomial_F.tsv the "ipyrad binomial" series used in
#                                                     figures/build_polyrad_F.py
#           stdout: all of the above, human readable
#
# ESTIMATORS  Two within-stand F estimators exist in the 2026 scripts on this object:
#   "valid" (hobs_paralog_filter.R, fb_F_by_depth.R): stand allele frequency from genotypes
#           with >= 6 reads only; this gives the drawn S2E baseline 1.207 and FreeBayes 1.231.
#   "fixed" (polyrad_F_fixed.R): numerator counts genotypes at ANY depth, denominator only those
#           with >= 6 reads; this gives the 1.175 drawn as "beta-binomial overdispersion".
# Both are computed here for every row, so each row can be compared on a single estimator.
#
# SPREAD RULE  bins with > 2000 genotypes (hobs_paralog_filter.R rule) for the main columns.
#
# GENOTYPE MODELS  polyRAD posterior (d = 5, IteratePopStructLD; GetProbableGenotypes) and, on the
#   SAME loci and read counts, the binomial call rule of scripts/ref/*.py (het if L(0.5) >
#   L(err 0.001)), to reproduce the ipyrad-binomial series and give the polyRAD model row its
#   own baseline.
# Inputs: results/polyrad/rad_genotyped.rds, metadata/popmap_177.tsv
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
suppressPackageStartupMessages(library(polyRAD))
P   <- "/home/ecastron/araucaria"
O   <- file.path(P, "results/derived/s2e"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
r   <- readRDS(file.path(P, "results/polyrad/rad_genotyped.rds"))

G   <- GetProbableGenotypes(r, omit1allelePerLocus = TRUE)$genotypes          # taxa x alleles
gl  <- r$alleles2loc[match(colnames(G), colnames(r$alleleDepth))]
D   <- matrix(0, nrow(G), ncol(G), dimnames = dimnames(G))
for (L in unique(gl)) D[, which(gl == L)] <- rowSums(r$alleleDepth[, which(r$alleles2loc == L), drop = FALSE])
A   <- r$alleleDepth[, colnames(G)]                                             # reads for that allele

# binomial calls from the same read counts: dosage of THIS allele (0, 1, 2)
# p is always 0.5, E or 1 - E here, so the p = 0 / p = 1 branches of the Python lnb() never apply
lnb <- function(k, n, p) k * log(p) + (n - k) * log(1 - p)
E <- 0.001
Bn <- ifelse(lnb(A, D, 0.5) > pmax(lnb(A, D, E), lnb(A, D, 1 - E)), 1, ifelse(A > D - A, 2, 0))
Bn[D == 0] <- NA

pop <- read.table(file.path(P, "metadata/popmap_177.tsv"), stringsAsFactors = FALSE)
stv <- setNames(pop[[2]], pop[[1]])[sub("\\.[12]$", "", rownames(G))]
EDG <- c(6, 10, 15, 20, 30, 50, 100, Inf); NM <- c("6-10", "10-15", "15-20", "20-30", "30-50", "50-100", "100+")

Fdepth <- function(Gm, keepLoci, estimator = c("valid", "fixed")) {
  estimator <- match.arg(estimator)
  cols <- which(gl %in% keepLoci)
  ho <- he <- ng <- numeric(length(NM))
  for (s in unique(stv[!is.na(stv)])) {
    rows <- which(stv == s); if (length(rows) < 6) next
    Gs <- Gm[rows, cols, drop = FALSE]; Ds <- D[rows, cols, drop = FALSE]
    v <- !is.na(Gs) & Ds >= 6; n_ok <- colSums(v)
    num <- if (estimator == "valid") colSums(ifelse(v, Gs, 0), na.rm = TRUE) else colSums(Gs, na.rm = TRUE)
    q <- ifelse(n_ok > 0, num / (2 * n_ok), NA)
    ok <- which(is.finite(q) & q > 0 & q < 1 & n_ok >= 6)
    for (b in seq_along(NM)) {
      inb <- Ds[, ok, drop = FALSE] >= EDG[b] & Ds[, ok, drop = FALSE] < EDG[b + 1] & !is.na(Gs[, ok, drop = FALSE])
      ho[b] <- ho[b] + sum((Gs[, ok, drop = FALSE] == 1) & inb, na.rm = TRUE)
      he[b] <- he[b] + sum(sweep(inb, 2, 2 * q[ok] * (1 - q[ok]), "*")); ng[b] <- ng[b] + sum(inb)
    }
  }
  list(F = 1 - ho / he, n = ng)
}
spr <- function(z, minn = 2000) { g <- z$n > minn; c(spread = max(z$F[g]) - min(z$F[g]), bins = sum(g)) }

all_loci <- sort(unique(gl))
cat(sprintf("polyRAD loci: %d, alleles scored: %d, taxa: %d\n", length(all_loci), ncol(G), nrow(G)))

# ---- locus metrics (hobs_paralog_filter.R and locus_structural_filter.R, unchanged logic) ----
W <- GetWeightedMeanGenotypes(r)
a2w <- r$alleles2loc[match(colnames(W), colnames(r$alleleDepth))]
locus_hobs <- tapply(colMeans(W > 0.2 & W < 0.8, na.rm = TRUE), a2w, max)
ld <- t(r$locDepth); pres <- ld; pres[pres == 0] <- NA
lmean <- rowMeans(pres, na.rm = TRUE); lcv <- apply(pres, 1, sd, na.rm = TRUE) / lmean
lmiss <- rowSums(ld == 0, na.rm = TRUE) / ncol(ld)
loc_idx <- seq_len(nrow(ld))      # locDepth column order = locus index used in alleles2loc

# ---- 1. figdata/locus_filter_F.tsv -------------------------------------------
SF <- list(c(Inf, Inf, Inf), c(0.80, 95, 0.70), c(0.80, 64, 0.65), c(0.60, 64, 0.65), c(0.55, 43, 0.60))
lf <- NULL
for (t in SF) {
  keep <- if (is.infinite(t[1])) all_loci else loc_idx[which(lcv <= t[1] & lmean <= t[2] & lmiss <= t[3])]
  z <- Fdepth(G, keep, "valid")
  lab <- if (is.infinite(t[1])) "no locus filter" else sprintf("cv<=%.2f depth<=%g miss<=%.2f", t[1], t[2], t[3])
  lf <- rbind(lf, data.frame(filter = lab, loci = length(keep), pct = round(100 * length(keep) / length(all_loci)),
              t(setNames(round(z$F, 3), NM)), spread = round(spr(z)[["spread"]], 3),
              t(setNames(z$n, paste0("n_", NM))), check.names = FALSE))
}
write.table(lf, file.path(O, "locus_filter_F.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== locus_filter_F (valid estimator, bins > 2000 genotypes) ==\n"); print(lf, row.names = FALSE)

# ---- 2. ipyrad-binomial series on the polyRAD loci ------------------------------------------------
zb_valid <- Fdepth(Bn, all_loci, "valid"); zb_fixed <- Fdepth(Bn, all_loci, "fixed")
ib <- data.frame(depth = NM, n_valid = zb_valid$n, F_valid = round(zb_valid$F, 4),
                 n_fixed = zb_fixed$n, F_fixed = round(zb_fixed$F, 4),
                 F_hardcoded = c(0.8694, 0.8153, 0.7847, 0.6512, 0.3061, -0.1960, -0.3363))
write.table(ib, file.path(O, "ipyrad_binomial_F.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== binomial calls on the 34,674 polyRAD loci vs the series in build_polyrad_F.py ==\n"); print(ib, row.names = FALSE)

# ---- 3. S2E rows on this route ----------------------------------------------------------------------
rows <- list(
  list("none - unfiltered baseline (polyRAD posterior)", G, all_loci),
  list("H_obs ceiling 0.45", G, as.integer(names(locus_hobs))[locus_hobs <= 0.45]),
  list("H_obs ceiling 0.20", G, as.integer(names(locus_hobs))[locus_hobs <= 0.20]),
  list("structural filter cv<=0.80 depth<=95 miss<=0.70", G, loc_idx[which(lcv <= 0.80 & lmean <= 95 & lmiss <= 0.70)]),
  list("structural filter cv<=0.80 depth<=64 miss<=0.65", G, loc_idx[which(lcv <= 0.80 & lmean <= 64 & lmiss <= 0.65)]),
  list("binomial calls, same loci (genotype-model baseline)", Bn, all_loci))
out <- NULL
for (est in c("valid", "fixed")) {
  zb <- Fdepth(G, all_loci, est); sb <- spr(zb)
  zbin <- Fdepth(Bn, all_loci, est); sbin <- spr(zbin)
  for (x in rows) {
    z <- Fdepth(x[[2]], x[[3]], est); s <- spr(z)
    out <- rbind(out, data.frame(mechanism = x[[1]], estimator = est, n_loci_kept = length(x[[3]]),
      pct_loci_kept = round(100 * length(x[[3]]) / length(all_loci), 1), n_bins = s[["bins"]],
      baseline_spread = round(sb[["spread"]], 4), filtered_spread = round(s[["spread"]], 4),
      reduction = round(sb[["spread"]] - s[["spread"]], 4), binomial_same_loci_spread = round(sbin[["spread"]], 4),
      genotypes_used = sum(z$n), t(setNames(round(z$F, 4), paste0("F_", NM))), t(setNames(z$n, paste0("n_", NM))),
      check.names = FALSE))
  }
}
write.table(out, file.path(O, "polyrad_route.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== polyRAD-route rows ==\n"); print(out[, 1:10], row.names = FALSE)
cat("written ->", O, "\n")
