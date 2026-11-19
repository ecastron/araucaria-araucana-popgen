#!/usr/bin/env Rscript
# Q6, RESPECIFIED: per-ASV heritability with the FULL relatedness matrix.
#
# WHY THE PREVIOUS SPECIFICATION WAS WEAK, and this is not a fishing expedition after a null:
#   (1) COMPRESSION. dbRDA needs sample-level predictors, so 78,802 SNPs were reduced to pairwise
#       theta and then to FOUR PCoA axes. Here the kinship matrix enters WHOLE, as the covariance
#       of a random effect -- no compression at any step.
#   (2) DILUTION. A community-level test asks whether OVERALL composition tracks relatedness. An
#       effect concentrated in a handful of taxa is buried by the hundreds that are not. The
#       Pinus radiata root study found 83 bacterial ASVs with h2 = 0.2-0.46 in data where the
#       whole-community effect was ABSENT. That is precisely the pattern the old Q6 could not see.
#
# MODEL, per ASV:  y = X beta + u + e,  u ~ N(0, K sigma2_g),  e ~ N(0, I sigma2_e)
#   y  CLR abundance across trees
#   X  intercept + stand + age class + 4 PCNM spatial vectors   (the Q6 covariates, unchanged)
#   K  ngsRelate coancestry, used in full
#   h2 = sigma2_g / (sigma2_g + sigma2_e)
#
# INFERENCE BY PERMUTATION, not by a likelihood-ratio boundary approximation. Testing a variance
# component at zero puts the parameter on the boundary of its space, where the usual chi-square
# reference is wrong. Tree labels are permuted against K, so the covariate structure and the
# between-ASV correlation are preserved. ONE permutation is shared across all ASVs per replicate,
# which keeps the joint null honest rather than treating ASVs as independent.
#
# TWO TESTS, BOTH DECLARED BEFORE RUNNING:
#   GLOBAL  mean h2 across ASVs vs its permutation null. This is the primary test -- it asks
#           whether there is ANY host genetic control, and it does not pay a multiple-testing cost.
#   PER-ASV empirical p per ASV, then Benjamini-Hochberg. Secondary, and expected to be
#           underpowered at n = 33-48; a null here is uninformative and will be reported as such.

suppressPackageStartupMessages({library(rrBLUP); library(vegan); library(parallel)})
set.seed(2026)
P <- "/home/ecastron/araucaria"; NPERM <- 999
CORES <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
dir.create(file.path(P, "results/Q6"), showWarnings = FALSE, recursive = TRUE)
sink(file.path(P, "results/Q6/heritability.txt"), split = TRUE)
cat(sprintf("permutations %d | cores %d\n", NPERM, CORES))

man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
rad <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
ch  <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)
rf  <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
to <- setNames(rad$tree_id, rad$run_accession)
ta <- to[bare(as.character(rf$ida))]; tb <- to[bare(as.character(rf$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta != tb & rf$nSites >= 2000 & is.finite(rf$theta)
kt <- sort(unique(c(ta[gd], tb[gd])))
K0 <- matrix(NA_real_, length(kt), length(kt), dimnames = list(kt, kt))
for (z in which(gd)) { K0[ta[z], tb[z]] <- rf$theta[z]; K0[tb[z], ta[z]] <- rf$theta[z] }
diag(K0) <- max(rf$theta[gd], na.rm = TRUE)

clr <- function(m) { m <- as.matrix(m); if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

pool <- function(amp) {
  st <- readRDS(file.path(P, sprintf("microbiome/%s/seqtab_nochim.rds", amp)))
  tx <- as.data.frame(readRDS(file.path(P, sprintf("microbiome/%s/taxonomy.rds", amp))), stringsAsFactors = FALSE)
  if (amp == "16S") {
    o <- rep(FALSE, ncol(st))
    for (rk in colnames(tx)) o <- o | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case = TRUE)
    o[is.na(o)] <- FALSE; st <- st[, !o, drop = FALSE]; tx <- tx[!o, , drop = FALSE]
  }
  i <- match(rownames(st), man$run_accession)
  k <- man$compartment[i] %in% c("branch", "needle") & !is.na(man$tree_id[i])
  st <- st[k, , drop = FALSE]; tid <- man$tree_id[i][k]; tr <- sort(unique(tid))
  list(cnt = t(sapply(tr, function(x) colSums(st[tid == x, , drop = FALSE]))), trees = tr, tax = tx)
}

h2_of <- function(y, X, K) {
  r <- try(mixed.solve(y = y, X = X, K = K, method = "REML"), silent = TRUE)
  if (inherits(r, "try-error") || !is.finite(r$Vu) || !is.finite(r$Ve)) return(NA_real_)
  if (r$Vu + r$Ve <= 0) return(NA_real_)
  r$Vu / (r$Vu + r$Ve)
}

for (amp in c("16S", "ITS")) {
  b <- pool(amp)
  sel <- rowSums(b$cnt) >= 250
  cnt <- b$cnt[sel, , drop = FALSE]; tr <- b$trees[sel]
  g <- tr %in% rownames(K0); cnt <- cnt[g, , drop = FALSE]; tg <- tr[g]
  cnt <- cnt[, colSums(cnt > 0) >= 5 & colSums(cnt)/sum(cnt) >= 1e-4, drop = FALSE]
  n <- nrow(cnt); m <- ncol(cnt)
  cat(sprintf("\n################ %s above-ground: %d trees x %d ASVs ################\n", amp, n, m))
  if (n < 25 || m < 10) { cat("insufficient; skipped\n"); next }

  Y <- clr(cnt)
  K <- K0[tg, tg]; K[is.na(K)] <- mean(K[upper.tri(K)], na.rm = TRUE); diag(K) <- max(K, na.rm = TRUE)
  K <- (K + t(K)) / 2                     # enforce symmetry after imputation
  mg <- man[match(tg, man$tree_id), ]; j <- match(tg, ch$tree_id)
  mem <- pcnm(dist(cbind(ch$lon[j], ch$lat[j])))$vectors[, 1:4, drop = FALSE]
  colnames(mem) <- paste0("S", 1:4)
  X <- model.matrix(~ factor(mg$site_code) + factor(mg$age_class) + mem)
  qr_ <- qr(X); X <- X[, qr_$pivot[seq_len(qr_$rank)], drop = FALSE]   # drop aliased columns
  cat(sprintf("fixed-effect columns retained: %d (rank)\n", ncol(X)))

  obs <- vapply(seq_len(m), function(k) h2_of(Y[, k], X, K), numeric(1))
  ok <- is.finite(obs)
  cat(sprintf("ASVs fitted: %d of %d | mean h2 = %.4f | median %.4f | max %.4f\n",
              sum(ok), m, mean(obs[ok]), median(obs[ok]), max(obs[ok])))

  perm <- simplify2array(mclapply(seq_len(NPERM), function(bi) {
    o <- sample(n)
    vapply(seq_len(m), function(k) h2_of(Y[o, k], X, K), numeric(1))
  }, mc.cores = CORES))                                   # m x NPERM

  # ---- GLOBAL test ----
  nullmean <- colMeans(perm, na.rm = TRUE)
  gp <- (1 + sum(nullmean >= mean(obs[ok]), na.rm = TRUE)) / (1 + sum(is.finite(nullmean)))
  cat(sprintf("\n[GLOBAL] observed mean h2 = %.4f | permutation null mean = %.4f (sd %.4f)\n",
              mean(obs[ok]), mean(nullmean, na.rm = TRUE), sd(nullmean, na.rm = TRUE)))
  cat(sprintf("[GLOBAL] p = %.4f  -> %s\n", gp,
              if (gp < 0.05) "host genetic control DETECTED across the community" else
                             "no community-wide host genetic control"))

  # ---- per-ASV ----
  pv <- vapply(seq_len(m), function(k) {
    if (!ok[k]) return(NA_real_)
    (1 + sum(perm[k, ] >= obs[k], na.rm = TRUE)) / (1 + sum(is.finite(perm[k, ])))
  }, numeric(1))
  q <- p.adjust(pv, "BH")
  cat(sprintf("[PER-ASV] p < 0.05: %d | BH q < 0.05: %d | q < 0.10: %d\n",
              sum(pv < 0.05, na.rm = TRUE), sum(q < 0.05, na.rm = TRUE), sum(q < 0.10, na.rm = TRUE)))
  hit <- which(q < 0.10)
  if (length(hit)) {
    tx <- b$tax[match(colnames(cnt)[hit], rownames(b$tax)), , drop = FALSE]
    cat("\n  taxa with q < 0.10:\n")
    for (z in seq_along(hit)) cat(sprintf("   h2=%.3f p=%.4f q=%.3f  %s\n", obs[hit[z]], pv[hit[z]], q[hit[z]],
        paste(na.omit(unlist(tx[z, ])), collapse = "; ")))
  }
  saveRDS(list(h2 = obs, p = pv, q = q, perm = perm, asv = colnames(cnt), n = n),
          file.path(P, sprintf("results/Q6/herit_%s.rds", amp)))
}
cat("\nA null GLOBAL test at n = 33-48 is weak evidence, not proof of absence; the estimator is\n")
cat("attenuated by kinship measurement error, which was declared before any of this was run.\n")
sink()
