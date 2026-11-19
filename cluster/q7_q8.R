#!/usr/bin/env Rscript
# Q7: is the host genetic effect DIFFERENTIAL between fungi and bacteria?
# Q8: does that effect operate at the spatial scale of gene dispersal estimated in Q2?
#
# Q7 IS REPORTED AS EFFECT SIZES WITH INTERVALS, NOT A SIGNIFICANCE VERDICT. A
# difference-of-effects test needs roughly four times the power of a single effect, and this
# design does not have it: 45-49 trees per stratum. A null difference here would be
# uninformative, so the honest output is the two effect sizes, their uncertainty, and the
# uncertainty on their difference -- stated before running so a null cannot be spun as evidence
# of no difference.
#
# The prior runs the OTHER WAY. Published work maps host control onto FUNGAL communities: immune
# receptor genes in switchgrass, and fungal-but-not-bacterial effects in Pinus radiata roots.
# Q6 found the opposite -- bacterial, in needles, with fungi null. The design says test the
# direction, do not fit to it, so the comparison is made symmetrically.
#
# Q8 IS EXPLORATORY and declared as such. It asks whether endophyte communities decay over the
# same distances as kinship. Q3 already showed the design has almost no pairs below 25 m, so the
# fine-scale end of this comparison is unresolvable; what can be compared is the decay slope
# across the range that IS sampled.

suppressPackageStartupMessages({library(vegan)})
set.seed(1)
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/Q7Q8")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "q7q8_summary.txt"), split = TRUE)

# ---- kinship matrix, as built for Q2/Q6 ----------------------------------------------------
rel <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
radman <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
tree_of <- setNames(radman$tree_id, radman$run_accession)
ta <- tree_of[bare(as.character(rel$ida))]; tb <- tree_of[bare(as.character(rel$idb))]
good <- !is.na(ta) & !is.na(tb) & ta != tb & rel$nSites >= 2000 & is.finite(rel$theta)
trees <- sort(unique(c(ta[good], tb[good])))
K <- matrix(NA_real_, length(trees), length(trees), dimnames = list(trees, trees))
for (i in which(good)) { K[ta[i], tb[i]] <- rel$theta[i]; K[tb[i], ta[i]] <- rel$theta[i] }
diag(K) <- max(rel$theta[good], na.rm = TRUE)

hav <- function(la1, lo1, la2, lo2) {
  R <- 6371000; p1 <- la1*pi/180; p2 <- la2*pi/180
  dp <- p2-p1; dl <- (lo2-lo1)*pi/180
  2*R*asin(sqrt(sin(dp/2)^2 + cos(p1)*cos(p2)*sin(dl/2)^2))
}

prep <- function(amp, comp) {
  f <- file.path(P, "results/Q5", sprintf("q5_%s_%s_objects.rds", amp, comp))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f); m <- o$manifest
  keep <- m$tree_id %in% rownames(K)
  X <- o$clr[keep, , drop = FALSE]; m <- m[keep, ]
  dup <- !duplicated(m$tree_id); X <- X[dup, , drop = FALSE]; m <- m[dup, ]
  if (nrow(m) < 25) return(NULL)
  Ks <- K[m$tree_id, m$tree_id]
  Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm = TRUE); diag(Ks) <- max(Ks, na.rm = TRUE)
  list(X = X, m = m, K = Ks, D = dist(X))
}

# partial R2 for host genetics, everything else removed -- the Q6 statistic
partR2 <- function(d, idx = seq_len(nrow(d$m))) {
  m <- d$m[idx, ]; X <- d$X[idx, , drop = FALSE]; Ks <- d$K[idx, idx]
  if (length(unique(m$site_code)) < 2) return(NA_real_)
  pc <- cmdscale(as.dist(max(Ks) - Ks), k = min(4, length(idx) - 2))
  colnames(pc) <- paste0("G", seq_len(ncol(pc)))
  mem <- tryCatch({ p <- pcnm(dist(cbind(m$longitude, m$latitude)))
                    p$vectors[, seq_len(min(4, ncol(p$vectors))), drop = FALSE] },
                  error = function(e) NULL)
  dat <- data.frame(pc, stand = factor(m$site_code), age = factor(m$age_class),
                    logdepth = log(rowSums(exp(X))))
  if (!is.null(mem)) { colnames(mem) <- paste0("S", seq_len(ncol(mem))); dat <- cbind(dat, mem) }
  cond <- paste("Condition(stand + logdepth + age",
                if (!is.null(mem)) paste("+", paste(colnames(mem), collapse = " + ")) else "", ")")
  fo <- as.formula(paste("dist(X) ~", paste(colnames(pc), collapse = " + "), "+", cond))
  cap <- try(dbrda(fo, data = dat), silent = TRUE)
  if (inherits(cap, "try-error") || is.null(cap$CCA)) return(NA_real_)
  cap$CCA$tot.chi / cap$tot.chi
}

cat("Q7: is the host genetic effect differential between fungi and bacteria?\n\n")
NB <- 300
res <- list()
for (comp in c("branch", "needle")) {
  cat(sprintf("---------- %s ----------\n", comp))
  for (amp in c("16S", "ITS")) {
    d <- prep(amp, comp)
    if (is.null(d)) { cat(sprintf("  %s: insufficient data\n", amp)); next }
    obs <- partR2(d)
    bs <- replicate(NB, {
      i <- sample(nrow(d$m), replace = TRUE)
      if (length(unique(d$m$site_code[i])) < 2) return(NA_real_)
      partR2(d, i) })
    bs <- bs[is.finite(bs)]
    ci <- quantile(bs, c(.025, .975), na.rm = TRUE)
    res[[paste(comp, amp)]] <- list(obs = obs, bs = bs)
    cat(sprintf("  %-4s n=%d  partial R2 = %.4f   95%% CI [%.4f, %.4f]  (%d bootstraps)\n",
                amp, nrow(d$m), obs, ci[1], ci[2], length(bs)))
  }
  a <- res[[paste(comp, "16S")]]; b <- res[[paste(comp, "ITS")]]
  if (!is.null(a) && !is.null(b)) {
    n <- min(length(a$bs), length(b$bs))
    diff <- a$bs[seq_len(n)] - b$bs[seq_len(n)]
    dci <- quantile(diff, c(.025, .975), na.rm = TRUE)
    cat(sprintf("\n  16S minus ITS: %+.4f   95%% CI [%+.4f, %+.4f]\n",
                a$obs - b$obs, dci[1], dci[2]))
    cat(sprintf("  interval %s zero -> %s\n",
                if (dci[1] < 0 && dci[2] > 0) "INCLUDES" else "excludes",
                if (dci[1] < 0 && dci[2] > 0)
                  "no detectable difference between fungi and bacteria" else
                  "a differential effect is detectable"))
  }
  cat("\n")
}
cat("Report these as effect sizes. A difference-of-effects test at n = 45-49 has low power, so\n")
cat("an interval spanning zero is uninformative rather than evidence of no difference.\n")

# ---- Q8 -------------------------------------------------------------------------------------
cat("\n\nQ8: does the community pattern operate at the gene-dispersal scale? (EXPLORATORY)\n\n")
cat(sprintf("kinship decay from Q2: b_F = -0.006742 per ln(m), Sp = 0.0071, Nb = 141\n\n"))
cat(sprintf("  %-14s %6s %12s %12s %10s\n", "stratum", "pairs", "slope/ln(m)", "p", "range (m)"))
for (comp in c("branch", "needle")) for (amp in c("16S", "ITS")) {
  d <- prep(amp, comp); if (is.null(d)) next
  m <- d$m; Dm <- as.matrix(d$D)
  iu <- which(upper.tri(Dm), arr.ind = TRUE)
  same <- m$site_code[iu[,1]] == m$site_code[iu[,2]]
  iu <- iu[same, , drop = FALSE]
  if (nrow(iu) < 40) { cat(sprintf("  %-14s %6d   (too few within-stand pairs)\n",
                                   paste(amp, comp), nrow(iu))); next }
  gd <- hav(m$latitude[iu[,1]], m$longitude[iu[,1]], m$latitude[iu[,2]], m$longitude[iu[,2]])
  cd <- Dm[iu]; st <- m$site_code[iu[,1]]
  ok <- gd > 0 & is.finite(gd) & is.finite(cd)
  gd <- gd[ok]; cd <- cd[ok]; st <- st[ok]
  fit <- lm(cd ~ log(gd) + factor(st))
  co <- summary(fit)$coefficients["log(gd)", ]
  cat(sprintf("  %-14s %6d %12.4f %12.4f %5.0f-%.0f\n", paste(amp, comp), length(gd),
              co[1], co[4], min(gd), max(gd)))
}
cat("\nA positive slope means communities become MORE dissimilar with distance inside a stand,\n")
cat("which is the community analogue of the kinship decline. Comparing the two scales directly\n")
cat("is not possible here: Q3 showed the design has 9 pairs below 25 m, so the distances where\n")
cat("kin structure is strongest are barely sampled. Q8 therefore remains exploratory.\n")
sink()
