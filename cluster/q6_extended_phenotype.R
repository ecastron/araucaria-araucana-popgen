#!/usr/bin/env Rscript
# Q6: does host relatedness predict endophyte community composition WITHIN stands, independent
# of geographic distance and tree size? This is the extended-phenotype test.
#
# METHOD IS dbRDA, NOT MANTEL. The Mantel test has poor power for this class of question and
# inflated type I error when both matrices are spatially autocorrelated -- which is exactly our
# situation, since Q2 established that kinship declines with distance. Kinship is inherently
# pairwise, so it enters the ordination as principal coordinates of the kinship matrix rather
# than as a distance matrix to be correlated.
#
# SIZE IS FITTED JOINTLY WITH RELATEDNESS, not sequentially. age_class and height are the same
# variable (Spearman rho = 0.968) so only one is used, and age_class is complete where height is
# not. Seedlings are spatially aggregated, so size and relatedness may be collinear; sequential
# testing would hand shared variance to whichever term went first.
#
# DIRECTION OF BIAS, DECLARED BEFORE RUNNING. Kinship carries measurement error from heterozygote
# under-calling (inferred F = 0.70, impossible in a dioecious species). Depth was verified not to
# track within-stand distance (p = 0.47-0.99), so that error is noise rather than trend. Noise in
# a PREDICTOR causes regression dilution: it attenuates the estimate toward zero. Therefore a
# POSITIVE result here is trustworthy, and a NULL result is ambiguous -- either no effect, or a
# real effect flattened by noisy kinship. This asymmetry must be stated with whatever comes out.

suppressPackageStartupMessages({library(vegan)})
set.seed(1)
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/Q6")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "q6_summary.txt"), split = TRUE)
cat("Q6 extended phenotype: host relatedness and endophyte community\n\n")

# ---- kinship, tree x tree --------------------------------------------------------------
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
cat(sprintf("kinship matrix: %d trees, %.1f%% of off-diagonal pairs estimated\n",
            nrow(K), 100 * mean(!is.na(K[upper.tri(K)]))))

man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
ch  <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)

run_one <- function(amp, comp) {
  cat(sprintf("\n---------- %s / %s ----------\n", amp, comp))
  f <- file.path(P, "results/Q5", sprintf("q5_%s_%s_objects.rds", amp, comp))
  if (!file.exists(f)) { cat("  Q5 objects absent\n"); return(invisible(NULL)) }
  o <- readRDS(f)
  m <- o$manifest
  tr <- m$tree_id
  keep <- tr %in% rownames(K)
  cat(sprintf("  microbiome libraries: %d | with host genotype: %d\n", length(tr), sum(keep)))
  if (sum(keep) < 25) { cat("  fewer than 25 paired trees; not fitted\n"); return(invisible(NULL)) }

  X <- o$clr[keep, , drop = FALSE]; m <- m[keep, ]; tr <- tr[keep]
  # one library per tree: duplicates would be pseudoreplication in a tree-level model
  dup <- !duplicated(tr); X <- X[dup, , drop = FALSE]; m <- m[dup, ]; tr <- tr[dup]
  cat(sprintf("  unique trees: %d across %d stands\n", length(tr), length(unique(m$site_code))))
  if (length(tr) < 25) { cat("  too few unique trees\n"); return(invisible(NULL)) }

  Ksub <- K[tr, tr]
  # a few pairs lack an estimate; fill with the observed mean so PCoA is defined, and report it
  miss <- mean(is.na(Ksub[upper.tri(Ksub)]))
  Ksub[is.na(Ksub)] <- mean(Ksub[upper.tri(Ksub)], na.rm = TRUE)
  diag(Ksub) <- max(Ksub, na.rm = TRUE)
  cat(sprintf("  kinship cells imputed with the mean: %.1f%%\n", 100 * miss))

  # kinship -> distance -> principal coordinates. These are sample-level predictors, which is
  # what dbRDA needs, and they avoid the Mantel framework entirely.
  Kd <- as.dist(max(Ksub) - Ksub)
  pc <- cmdscale(Kd, k = min(4, length(tr) - 2), eig = TRUE)
  gax <- pc$points
  colnames(gax) <- paste0("G", seq_len(ncol(gax)))
  posev <- pc$eig[pc$eig > 0]
  cat(sprintf("  genetic PCoA axes retained: %d (%.1f%% of positive eigenvalue mass)\n",
              ncol(gax), 100 * sum(posev[seq_len(ncol(gax))]) / sum(posev)))

  D <- dist(X)
  xy <- cbind(m$longitude, m$latitude)
  mem <- tryCatch({ p <- pcnm(dist(xy)); p$vectors[, seq_len(min(4, ncol(p$vectors))), drop = FALSE] },
                  error = function(e) NULL)
  age <- factor(m$age_class)
  logdepth <- log(rowSums(exp(X)))       # CLR back-scale is monotone in library size
  stand <- factor(m$site_code)
  cat(sprintf("  age classes present: %s\n", paste(levels(age), collapse = ", ")))

  dat <- data.frame(gax, age = age, stand = stand, logdepth = logdepth)
  if (!is.null(mem)) { colnames(mem) <- paste0("S", seq_len(ncol(mem))); dat <- cbind(dat, mem) }
  gterms <- colnames(gax); sterms <- if (is.null(mem)) character(0) else colnames(mem)

  # WITHIN stands: stand is conditioned out, so nothing here rides on among-stand differences.
  cond <- "Condition(stand + logdepth)"
  form <- as.formula(paste("D ~", paste(c(gterms, sterms, "age"), collapse = " + "), "+", cond))
  cap <- try(dbrda(form, data = dat), silent = TRUE)
  if (inherits(cap, "try-error")) { cat("  dbRDA failed:", conditionMessage(attr(cap, "condition")), "\n"); return(invisible(NULL)) }

  cat("\n  [dbRDA] genetics + space + size, stand and depth conditioned out\n")
  a <- anova(cap, by = "terms", permutations = 999)
  print(a)
  cat(sprintf("\n  constrained inertia: %.1f%%\n", 100 * cap$CCA$tot.chi / cap$tot.chi))

  # the genetic term on its own, with space and size partialled out
  form_g <- as.formula(paste("D ~", paste(gterms, collapse = " + "), "+ Condition(stand + logdepth",
                             if (length(sterms)) paste("+", paste(sterms, collapse = " + ")) else "",
                             "+ age)"))
  capg <- try(dbrda(form_g, data = dat), silent = TRUE)
  if (!inherits(capg, "try-error")) {
    ag <- anova(capg, permutations = 999)
    r2 <- capg$CCA$tot.chi / capg$tot.chi
    cat(sprintf("\n  HOST GENETICS alone, space + size + stand + depth removed:\n"))
    cat(sprintf("    partial R2 = %.4f   F = %.3f   p = %.4f\n", r2, ag$F[1], ag$`Pr(>F)`[1]))
    verdict <- if (isTRUE(ag$`Pr(>F)`[1] < 0.05)) "DETECTED" else "not detected"
    cat(sprintf("    host genetic effect: %s\n", verdict))
  }

  if (length(sterms)) {
    vp <- try(varpart(D, gax, mem, ~ age), silent = TRUE)
    if (!inherits(vp, "try-error")) {
      cat("\n  [variance partitioning] genetics | space | size\n")
      print(vp$part$indfract[, "Adj.R.square", drop = FALSE], digits = 3)
    }
  }
  invisible(NULL)
}

for (amp in c("16S", "ITS")) for (comp in c("branch", "needle")) {
  tryCatch(run_one(amp, comp),
           error = function(e) cat("\n  ERROR", amp, comp, ":", conditionMessage(e), "\n"))
}

cat("\n================ reading this result ================\n")
cat("Kinship carries measurement error that is NOISE, not trend (depth does not track\n")
cat("within-stand distance: p = 0.47-0.99). Noise in a predictor attenuates toward zero, so:\n")
cat("  a POSITIVE result is trustworthy -- noise makes it harder, not easier, to obtain\n")
cat("  a NULL result is AMBIGUOUS -- either no effect, or one flattened by noisy kinship\n")
cat("This asymmetry was declared before the analysis ran.\n")
sink()
