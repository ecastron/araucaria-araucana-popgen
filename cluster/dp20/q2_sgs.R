#!/usr/bin/env Rscript
# Q2: fine-scale spatial genetic structure. Regress kinship on ln(geographic distance), summarise
# as Sp, and derive neighbourhood size Nb.
#
# POOLED ACROSS STANDS, with stand as a stratum -- never per stand. Each stand holds 10-20 trees,
# an order of magnitude below the published requirement: Cavers et al. (2005) put the lower limit
# near 100 individuals, Jump et al. (2006) at 150-200 in Fagus. Thousands of SNPs do not rescue
# this, because what limits the estimate is the number of spatially distributed PAIRS. Per-stand
# Sp values are therefore not reported (spec, Section 7).
#
# THE GATE IS ONE-DIRECTIONAL. A significant Sp makes the spatial term in Q6 load-bearing. A
# NON-significant result at this sample size is NOT evidence of no structure and must not be used
# to claim the microsite confound has no route. Read the other way, the gate converts low power
# into false confidence.
#
#   Sp = -b_F / (1 - F1)     b_F = regression slope of kinship on ln(distance)
#                            F1  = mean kinship in the first distance class
#   Nb = 1 / Sp
#
# sigma is NOT derived: that needs effective density, which was never measured.
#
# Permutation null: shuffle tree POSITIONS WITHIN stands, which preserves both the among-stand
# structure (real, and not what Q2 is about) and each stand's spatial extent.

suppressPackageStartupMessages({library(vegan)})
set.seed(1)
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/dp20/Q2")

rel <- read.delim(file.path(OUT, "ngsrelate.tsv"), check.names = FALSE)
cat(sprintf("ngsRelate output: %d pairs, columns: %s\n", nrow(rel),
            paste(head(names(rel), 12), collapse = ", ")))

# Use theta, the coancestry coefficient -- the direct analogue of the kinship coefficient SGS
# analysis is built on. NOT KING: on this data KING is numerically broken, with 13,232 of 15,576
# pairs outside |1.5| and a minimum of -8.5e10, because the low per-pair site count and high
# inferred inbreeding break its ratio estimator. theta is bounded (0 to 0.533 here) and behaves.
#
# The chosen column is range-checked before use: a kinship coefficient outside [-1, 1] is a
# broken estimator, not a finding, and must never reach a regression.
prefer <- c("theta", "rab", "R1")
kincol <- NA
for (cc in prefer) {
  if (cc %in% names(rel)) {
    v <- rel[[cc]]
    if (mean(abs(v) <= 1, na.rm = TRUE) > 0.99) { kincol <- cc; break }
    cat(sprintf("  rejecting %s: %.1f%% of values outside [-1,1]\n",
                cc, 100 * mean(abs(v) > 1, na.rm = TRUE)))
  }
}
if (is.na(kincol)) stop("no relatedness column passes the range check")
cat(sprintf("using relatedness column: %s\n", kincol))

# Per-pair site count varies 40-fold; pairs resting on very few shared sites carry almost no
# information and would add noise to the regression at full weight.
if ("nSites" %in% names(rel)) {
  cat(sprintf("shared sites per pair: median %d, %.1f%% below 2000\n",
              median(rel$nSites), 100 * mean(rel$nSites < 2000)))
}
if (all(c("Fa", "Fb") %in% names(rel))) {
  mf <- median(c(rel$Fa, rel$Fb), na.rm = TRUE)
  cat(sprintf("median inferred inbreeding F = %.3f\n", mf))
  if (mf > 0.3) {
    cat("  WARNING: F this high is implausible for an outcrossing conifer. It reflects\n")
    cat("  heterozygote under-calling at low depth, and it biases every relatedness estimate\n")
    cat("  derived from the same likelihoods. Treat Sp below as provisional.\n")
  }
}
cat("\n")

idcols <- intersect(c("ida", "idb", "a", "b"), names(rel))
if (length(idcols) < 2) stop("no sample identifier columns in ngsRelate output")
A <- as.character(rel[[idcols[1]]]); B <- as.character(rel[[idcols[2]]])
K <- rel[[kincol]]

man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
pop <- read.delim(file.path(P, "metadata/popmap_177.tsv"), header = FALSE,
                  col.names = c("lib", "stand"))
ch  <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)
radman <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)

bare <- function(x) sub("\\.[12]$", "", x)
tree_of <- setNames(radman$tree_id, radman$run_accession)
lat <- setNames(ch$lat, ch$tree_id); lon <- setNames(ch$lon, ch$tree_id)
stand_of <- setNames(pop$stand, pop$lib)

ta <- tree_of[bare(A)]; tb <- tree_of[bare(B)]
sa <- stand_of[bare(A)]; sb <- stand_of[bare(B)]
ok <- !is.na(ta) & !is.na(tb) & ta != tb & !is.na(sa) & sa == sb &
      !is.na(lat[ta]) & !is.na(lat[tb])
cat(sprintf("pairs usable (same stand, different trees, coordinates present): %d\n", sum(ok)))
if (sum(ok) < 50) stop("too few within-stand pairs for SGS")

hav <- function(la1, lo1, la2, lo2) {
  R <- 6371000; p1 <- la1*pi/180; p2 <- la2*pi/180
  dp <- p2-p1; dl <- (lo2-lo1)*pi/180
  2*R*asin(sqrt(sin(dp/2)^2 + cos(p1)*cos(p2)*sin(dl/2)^2))
}
d <- hav(lat[ta[ok]], lon[ta[ok]], lat[tb[ok]], lon[tb[ok]])
k <- K[ok]; st <- sa[ok]
nS <- if ("nSites" %in% names(rel)) rel$nSites[ok] else rep(Inf, sum(ok))
keep <- d > 0 & is.finite(d) & is.finite(k) & nS >= 2000
d <- d[keep]; k <- k[keep]; st <- st[keep]
cat(sprintf("after removing zero distances and pairs with <2000 shared sites: %d pairs, %d stands\n", length(d),
            length(unique(st))))
cat(sprintf("distance range: %.1f - %.1f m (median %.1f)\n\n", min(d), max(d), median(d)))

ld <- log(d)
# STAND MUST BE IN THE MODEL. Pooling raw across stands conflates the WITHIN-stand
# kinship-distance relationship (what Q2 asks) with the AMONG-stand one (real -- 15.8% of
# variance, isolation by distance -- but Q1's business). Stands differ in both mean kinship and
# mean pairwise distance, so without a stand term the among-stand component dominates and can
# even flip the sign of the slope. Observed here: pooled slope +0.000503, which would have been
# reported as SGS of the wrong sign.
fit <- lm(k ~ ld + factor(st))
b <- unname(coef(fit)["ld"])
fit_pooled <- lm(k ~ ld)
cat(sprintf("pooled slope, NO stand term (biased, shown for contrast): %+.6f\n",
            coef(fit_pooled)[2]))
q <- quantile(d, 0.2)
F1 <- mean(k[d <= q])
Sp <- -b / (1 - F1)
Nb <- if (Sp > 0) 1/Sp else NA

cat("--- kinship-distance regression (within stands, stand as stratum) ---\n")
print(summary(fit)$coefficients[1:2, , drop = FALSE])
cat(sprintf("\nF1 (mean kinship, first distance class <= %.0f m) = %.4f\n", q, F1))
cat(sprintf("b_F (slope on ln distance)                        = %.6f\n", b))
cat(sprintf("Sp = -b_F / (1 - F1)                              = %.6f\n", Sp))
cat(sprintf("Nb = 1 / Sp                                       = %s\n",
            if (is.na(Nb)) "undefined (Sp <= 0)" else sprintf("%.1f", Nb)))
cat("\nsigma is NOT reported: it requires effective density, which was not measured.\n")

# Permutation null: positions shuffled WITHIN stands.
nperm <- 999; bs <- numeric(nperm)
idx <- split(seq_along(d), st)
fst <- factor(st)
for (i in seq_len(nperm)) {
  kk <- k
  for (g in idx) kk[g] <- k[sample(g)]
  bs[i] <- unname(coef(lm(kk ~ ld + fst))["ld"])
}
pval <- (sum(bs <= b) + 1) / (nperm + 1)
cat(sprintf("\npermutation p (slope more negative than chance, positions shuffled within stands) = %.4f\n",
            pval))
if (pval < 0.05 && b < 0) {
  cat("SIGNIFICANT spatial genetic structure: kinship declines with distance within stands,\n")
  cat("so the spatial term in Q6 is load-bearing.\n")
} else if (b >= 0) {
  cat("NO spatial genetic structure detected: the within-stand slope is POSITIVE or zero, which\n")
  cat("is not SGS regardless of the permutation p-value. A significant p with a positive slope\n")
  cat("means the observed slope is unusual, not that kinship declines with distance.\n")
} else {
  cat("NOT significant. READ THIS IN ONE DIRECTION ONLY: at 10-20 trees per stand this is NOT\n")
  cat("evidence of absent structure, and must not be used to argue the microsite confound in Q6\n")
  cat("has no route.\n")
}

cat("\n--- kinship by distance class (descriptive) ---\n")
br <- unique(quantile(d, seq(0, 1, 0.1)))
cl <- cut(d, br, include.lowest = TRUE)
tab <- data.frame(class = levels(cl),
                  n = as.integer(table(cl)),
                  median_m = tapply(d, cl, median),
                  mean_kinship = tapply(k, cl, mean))
print(tab, row.names = FALSE, digits = 4)

write.csv(tab, file.path(OUT, "q2_kinship_by_distance.csv"), row.names = FALSE)
saveRDS(list(d = d, k = k, stand = st, Sp = Sp, Nb = Nb, b = b, F1 = F1, p = pval),
        file.path(OUT, "q2_sgs.rds"))
cat(sprintf("\nwritten: %s\n", file.path(OUT, "q2_sgs.rds")))
