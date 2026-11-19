#!/usr/bin/env Rscript
# WHERE DOES THE NEEDLE 16S EFFECT DIE?
#
# Three numbers now exist for the same question:
#   published, organellar included, prevalence>=2 : R2 = 0.101, p = 0.026  (0.007 first run)
#   organellar removed, same loose filter, 178 ASVs: R2 = 0.113, p = 0.013
#   organellar removed, prevalence>=5 & rel>=1e-4, 16 ASVs: R2 = 0.073, p = 0.271
#
# So the effect lives in ASVs that a routine prevalence/abundance filter discards, in libraries
# holding a median of 43 bacterial reads. Two explanations are indistinguishable so far:
#   (A) real signal in genuinely rare endophytes, thrown away by an over-aggressive filter
#   (B) structured noise -- at 43 reads the CLR is mostly zero-replacement, so distances are
#       driven by WHICH ASVs happen to be observed, not by their abundances
#
# TWO DIAGNOSTICS THAT SEPARATE THEM, declared before running:
#   1. FILTER GRADIENT. (A) predicts a gradual decline as rare ASVs are removed. (B) predicts
#      the effect tracks ASV COUNT itself -- more sparse columns, more apparent structure.
#   2. PRESENCE/ABSENCE SUBSTITUTION. Replace counts with 0/1 and refit. Under (A) discarding
#      abundance should WEAKEN the effect. Under (B) it should be unchanged or stronger,
#      because the signal was never in the abundances.
# Diagnostic 2 is the decisive one.

suppressPackageStartupMessages({library(vegan)})
set.seed(7)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q6/filter_gradient.txt"), split = TRUE)

st16 <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx16 <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org <- rep(FALSE, ncol(st16))
for (rk in colnames(tx16)) org <- org | grepl("Chloroplast|Mitochondri", tx16[[rk]], ignore.case = TRUE)
org[is.na(org)] <- FALSE; names(org) <- colnames(st16)

rel <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
radman <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
tree_of <- setNames(radman$tree_id, radman$run_accession)
ta <- tree_of[bare(as.character(rel$ida))]; tb <- tree_of[bare(as.character(rel$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta != tb & rel$nSites >= 2000 & is.finite(rel$theta)
trees <- sort(unique(c(ta[gd], tb[gd])))
K <- matrix(NA_real_, length(trees), length(trees), dimnames = list(trees, trees))
for (i in which(gd)) { K[ta[i], tb[i]] <- rel$theta[i]; K[tb[i], ta[i]] <- rel$theta[i] }
diag(K) <- max(rel$theta[gd], na.rm = TRUE)

clr <- function(m) { m <- as.matrix(m)
  if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

o <- readRDS(file.path(P, "results/Q5/q5_16S_needle_objects.rds"))
m <- o$manifest
keep <- !is.na(m$tree_id) & m$tree_id %in% rownames(K); m <- m[keep, ]
dup <- !duplicated(m$tree_id); m <- m[dup, ]
ids <- intersect(rownames(st16), m$run_accession); m <- m[match(ids, m$run_accession), ]
bac <- st16[ids, !org[colnames(st16)], drop = FALSE]
cat(sprintf("needle: %d trees | non-organellar reads median %.0f\n\n", nrow(m), median(rowSums(bac))))

Ks <- K[m$tree_id, m$tree_id]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm = TRUE)
diag(Ks) <- max(Ks, na.rm = TRUE)
pc <- cmdscale(as.dist(max(Ks) - Ks), k = 4); colnames(pc) <- paste0("G", 1:4)
mem <- pcnm(dist(cbind(m$longitude, m$latitude)))$vectors[, 1:4, drop = FALSE]
colnames(mem) <- paste0("S", 1:4)

fit <- function(X, nperm = 9999) {
  dat <- data.frame(pc, mem, stand = factor(m$site_code), age = factor(m$age_class),
                    logdepth = log(rowSums(exp(X))))
  cap <- try(dbrda(dist(X) ~ G1+G2+G3+G4 + Condition(stand + logdepth + S1+S2+S3+S4 + age),
                   data = dat), silent = TRUE)
  if (inherits(cap, "try-error") || is.null(cap$CCA)) return(NULL)
  a <- anova(cap, permutations = nperm)
  list(r2 = cap$CCA$tot.chi / cap$tot.chi, F = a$F[1], p = a$`Pr(>F)`[1])
}

cat("--- DIAGNOSTIC 1: filter gradient (abundance-weighted CLR) ---\n")
cat(sprintf("  %-12s %-12s %6s %10s %8s %9s\n", "prevalence", "rel.abund", "ASVs", "partial R2", "F", "p"))
grid <- expand.grid(pv = c(2, 3, 5, 10), ab = c(0, 1e-5, 1e-4))
for (i in seq_len(nrow(grid))) {
  B <- bac[, colSums(bac > 0) >= grid$pv[i] & colSums(bac) / sum(bac) >= grid$ab[i], drop = FALSE]
  if (ncol(B) < 5) { cat(sprintf("  >=%-10d %-12s %6d   (too few ASVs)\n", grid$pv[i],
                        format(grid$ab[i], scientific = TRUE), ncol(B))); next }
  r <- fit(clr(B))
  cat(sprintf("  >=%-10d %-12s %6d %10.4f %8.3f %9.4f %s\n", grid$pv[i],
              format(grid$ab[i], scientific = TRUE), ncol(B),
              r$r2, r$F, r$p, if (r$p < 0.05) "<<" else ""))
}

cat("\n--- DIAGNOSTIC 2 (DECISIVE): does the effect need the abundances? ---\n")
cat("  If presence/absence reproduces or beats CLR, the signal was in WHICH ASVs were seen,\n")
cat("  not how much of them -- i.e. sampling noise at 43 reads, not community structure.\n\n")
cat(sprintf("  %-12s %6s   %-22s %-22s\n", "prevalence", "ASVs", "CLR (abundance)", "presence/absence"))
for (pv in c(2, 3, 5)) {
  B <- bac[, colSums(bac > 0) >= pv, drop = FALSE]
  if (ncol(B) < 5) next
  rc <- fit(clr(B))
  PA <- (B > 0) * 1
  PA <- PA[, colSums(PA) > 0 & colSums(PA) < nrow(PA), drop = FALSE]
  rp <- if (ncol(PA) >= 5) fit(clr(PA + 0.5)) else NULL
  cat(sprintf("  >=%-10d %6d   R2=%.4f p=%.4f%-5s  %s\n", pv, ncol(B),
      rc$r2, rc$p, if (rc$p < 0.05) " <<" else "",
      if (is.null(rp)) "n/a" else sprintf("R2=%.4f p=%.4f%s", rp$r2, rp$p,
      if (rp$p < 0.05) " <<" else "")))
}
sink()
