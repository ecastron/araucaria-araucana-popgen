#!/usr/bin/env Rscript
# Q6 RE-RUN on decontaminated, properly filtered matrices, with a fragility test.
#
# The organellar control showed the needle 16S effect is NOT a host-DNA artifact: partitioning
# the published matrix gave organellar p = 0.082 and NON-organellar p = 0.013, R2 = 0.113. The
# host reads were diluting the signal, not creating it. So the effect survives -- but it now
# rests on a median of 43 non-organellar reads per needle library, which is very thin.
#
# THE FRAGILITY TEST, declared before running: re-fit at increasing floors on NON-ORGANELLAR
# depth. If the effect is carried by the near-empty libraries it will vanish as the floor rises.
# If it holds or strengthens on the better-sequenced subset, it is real. Losing significance at
# a higher floor is confounded with losing sample size, so n is reported at every step and a
# drop is NOT read as refutation on its own.
#
# 9999 permutations throughout: at 999 the same model returned p = 0.007, 0.026 and 0.013 on
# repeat runs, which is Monte Carlo noise large enough to change how a reader sees the result.

suppressPackageStartupMessages({library(vegan)})
set.seed(42)
MINTREES <- 8   # 25-tree floor removed on request 24 Aug 2026
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/nofloor/Q6_baseline/q6_rerun_clean.txt"), split = TRUE)

st16 <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx16 <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org <- rep(FALSE, ncol(st16))
for (rk in colnames(tx16)) org <- org | grepl("Chloroplast|Mitochondri", tx16[[rk]], ignore.case = TRUE)
org[is.na(org)] <- FALSE; names(org) <- colnames(st16)
stI <- readRDS(file.path(P, "microbiome/ITS/seqtab_nochim.rds"))

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

fit <- function(X, m, nperm = 9999) {
  tr <- m$tree_id
  Ks <- K[tr, tr]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm = TRUE)
  diag(Ks) <- max(Ks, na.rm = TRUE)
  pc <- cmdscale(as.dist(max(Ks) - Ks), k = min(4, length(tr) - 2))
  colnames(pc) <- paste0("G", seq_len(ncol(pc)))
  mem <- tryCatch({ p <- pcnm(dist(cbind(m$longitude, m$latitude)))
                    p$vectors[, seq_len(min(4, ncol(p$vectors))), drop = FALSE] },
                  error = function(e) NULL)
  dat <- data.frame(pc, stand = factor(m$site_code), age = factor(m$age_class),
                    logdepth = log(rowSums(exp(X))))
  s_ <- character(0)
  if (!is.null(mem)) { colnames(mem) <- paste0("S", seq_len(ncol(mem))); dat <- cbind(dat, mem); s_ <- colnames(mem) }
  if (nlevels(dat$stand) < 2) return(NULL)
  fo <- as.formula(paste("dist(X) ~", paste(colnames(pc), collapse = " + "),
        "+ Condition(stand + logdepth", if (length(s_)) paste("+", paste(s_, collapse = " + ")) else "", "+ age)"))
  cap <- try(dbrda(fo, data = dat), silent = TRUE)
  if (inherits(cap, "try-error") || is.null(cap$CCA)) return(NULL)
  a <- anova(cap, permutations = nperm)
  list(r2 = cap$CCA$tot.chi / cap$tot.chi, F = a$F[1], p = a$`Pr(>F)`[1])
}

# align a raw count table to the tree set Q6 uses
align <- function(raw, comp) {
  o <- readRDS(file.path(P, "results/Q5", sprintf("q5_%s_%s_objects.rds",
       if (identical(raw, st16)) "16S" else "ITS", comp)))
  m <- o$manifest
  keep <- !is.na(m$tree_id) & m$tree_id %in% rownames(K)
  m <- m[keep, ]; dup <- !duplicated(m$tree_id); m <- m[dup, ]
  ids <- intersect(rownames(raw), m$run_accession)
  m <- m[match(ids, m$run_accession), ]
  list(cnt = raw[ids, , drop = FALSE], m = m)
}

cat("=========== 16S, host organellar reads REMOVED ===========\n")
for (comp in c("needle", "branch", "root")) {
  a <- tryCatch(align(st16, comp), error = function(e) NULL)
  if (is.null(a) || nrow(a$m) < MINTREES) { cat(sprintf("\n-- %s: <%d genotyped trees, skipped\n", comp, MINTREES)); next }
  bac <- a$cnt[, !org[colnames(a$cnt)], drop = FALSE]
  d <- rowSums(bac)
  cat(sprintf("\n-- %s: %d trees | non-organellar reads median %.0f (range %.0f-%.0f)\n",
              comp, nrow(a$m), median(d), min(d), max(d)))
  cat(sprintf("   %-22s %5s %6s %9s %8s %9s\n", "depth floor", "n", "ASVs", "partial R2", "F", "p"))
  for (fl in c(0, 50, 100, 250, 500)) {
    sel <- d >= fl
    if (sum(sel) < MINTREES) { cat(sprintf("   >=%-20d %5d   (below %d trees, not fitted)\n", fl, sum(sel), MINTREES)); next }
    B <- bac[sel, , drop = FALSE]
    B <- B[, colSums(B > 0) >= 5 & colSums(B) / sum(B) >= 1e-4, drop = FALSE]
    if (ncol(B) < 10) { cat(sprintf("   >=%-20d %5d   (<10 ASVs survive filter)\n", fl, sum(sel))); next }
    r <- fit(clr(B), a$m[sel, ])
    if (is.null(r)) { cat(sprintf("   >=%-20d %5d   model failed\n", fl, sum(sel))); next }
    cat(sprintf("   >=%-20d %5d %6d %9.4f %8.3f %9.4f %s\n", fl, sum(sel), ncol(B),
                r$r2, r$F, r$p, if (!is.na(r$p) && r$p < 0.05) "<<" else ""))
  }
}

cat("\n\n=========== ITS (no contamination; prevalence>=5 & rel>=1e-4) ===========\n")
for (comp in c("needle", "branch", "root")) {
  a <- tryCatch(align(stI, comp), error = function(e) NULL)
  if (is.null(a) || nrow(a$m) < MINTREES) { cat(sprintf("\n-- %s: <%d genotyped trees, skipped\n", comp, MINTREES)); next }
  B <- a$cnt
  B <- B[, colSums(B > 0) >= 5 & colSums(B) / sum(B) >= 1e-4, drop = FALSE]
  r <- fit(clr(B), a$m)
  cat(sprintf("\n-- %-8s n=%d  ASVs=%4d (from %d)  ", comp, nrow(a$m), ncol(B), ncol(a$cnt)))
  if (is.null(r)) cat("model failed\n") else
    cat(sprintf("partial R2 = %.4f  F = %.3f  p = %.4f %s\n", r$r2, r$F, r$p, if (!is.na(r$p) && r$p < 0.05) "<<" else ""))
}
sink()
