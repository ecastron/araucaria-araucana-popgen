#!/usr/bin/env Rscript
# DECISIVE CONTROL for the Q6 headline result (16S needle, partial R2 = 0.101, p = 0.007).
#
# The 16S needle libraries are 99.9% Araucaria chloroplast + mitochondrial reads (median 27
# non-organellar reads per library). Organellar haplotypes ARE host genotype. So the Q6 model
# may have been regressing host genotype on host genotype.
#
# THE TEST, on the EXACT library set Q6 used, splitting the same counts two ways:
#   (a) ORGANELLAR ASVs only   -- if this reproduces or exceeds the effect, the association is
#                                 host DNA correlating with itself, not a microbial community.
#   (b) NON-ORGANELLAR only    -- the actual bacteria. If this is null, the published effect
#                                 came from (a).
# Both outcomes are informative and neither is the "hoped for" result; written down before running.

suppressPackageStartupMessages({library(vegan)})
set.seed(1)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q6/organelle_control.txt"), split = TRUE)

st <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
isorg <- rep(FALSE, ncol(st))
for (rk in colnames(tx)) isorg <- isorg | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case = TRUE)
isorg[is.na(isorg)] <- FALSE
cat(sprintf("organellar ASVs: %d of %d\n", sum(isorg), length(isorg)))

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

clr <- function(m) { m <- as.matrix(m)
  if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

fit_one <- function(X, m, tag) {
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
  st_ <- character(0)
  if (!is.null(mem)) { colnames(mem) <- paste0("S", seq_len(ncol(mem)))
                       dat <- cbind(dat, mem); st_ <- colnames(mem) }
  D <- dist(X)
  fo <- as.formula(paste("D ~", paste(colnames(pc), collapse = " + "),
        "+ Condition(stand + logdepth", if (length(st_)) paste("+", paste(st_, collapse = " + ")) else "", "+ age)"))
  cap <- try(dbrda(fo, data = dat), silent = TRUE)
  if (inherits(cap, "try-error") || is.null(cap$CCA)) { cat(sprintf("  %-28s model failed\n", tag)); return(invisible()) }
  a <- anova(cap, permutations = 999)
  cat(sprintf("  %-28s n=%d  ASVs=%4d  partial R2 = %.4f  F = %.3f  p = %.4f  %s\n",
              tag, nrow(m), ncol(X), cap$CCA$tot.chi / cap$tot.chi, a$F[1], a$`Pr(>F)`[1],
              if (isTRUE(a$`Pr(>F)`[1] < 0.05)) "<< DETECTED" else ""))
}

for (comp in c("needle", "branch")) {
  f <- file.path(P, "results/Q5", sprintf("q5_16S_%s_objects.rds", comp))
  if (!file.exists(f)) next
  o <- readRDS(f); m <- o$manifest
  keep <- m$tree_id %in% rownames(K) & !is.na(m$tree_id)
  m <- m[keep, ]; dup <- !duplicated(m$tree_id); m <- m[dup, ]
  libs <- m$run_accession
  if (is.null(libs) || !all(libs %in% rownames(st))) {
    hit <- intersect(rownames(st), unlist(m[sapply(m, is.character)]))
    cat(sprintf("\n[%s] could not key manifest to seqtab rownames; skipped\n", comp)); next
  }
  sub <- st[libs, , drop = FALSE]
  cat(sprintf("\n===== 16S / %s : %d trees =====\n", comp, nrow(m)))
  cat(sprintf("  median non-organellar reads per library: %.0f (host fraction %.2f%%)\n",
      median(rowSums(sub[, !isorg, drop = FALSE])),
      100 * sum(sub[, isorg, drop = FALSE]) / sum(sub)))
  A <- sub[, isorg, drop = FALSE]; A <- A[, colSums(A > 0) >= 2, drop = FALSE]
  B <- sub[, !isorg, drop = FALSE]; B <- B[, colSums(B > 0) >= 2, drop = FALSE]
  fit_one(clr(sub[, colSums(sub > 0) >= 2, drop = FALSE]), m, "ALL ASVs (as published)")
  fit_one(clr(A), m, "ORGANELLAR only")
  fit_one(clr(B), m, "NON-organellar (bacteria)")
}
cat("\nIf ORGANELLAR reproduces the effect and NON-organellar does not, the published Q6\n")
cat("result is host DNA predicting host DNA and must be withdrawn.\n")
sink()
