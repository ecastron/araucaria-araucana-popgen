#!/usr/bin/env Rscript
# Q5 RE-RUN. The stand-vs-climate comparison was computed on a 16S matrix that was 91.6% host
# organellar sequence, so it needs redoing on bacteria alone. ITS was never contaminated and is
# re-run only with the prevalence filter, to keep both amplicons on one footing.
#
# A depth floor of 100 NON-ORGANELLAR reads is applied to 16S and reported, because a library
# holding 43 bacterial reads cannot describe a community regardless of what the model returns.

suppressPackageStartupMessages({library(vegan)})
set.seed(11)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q5/q5_rerun_clean.txt"), split = TRUE)

st16 <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx16 <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org <- rep(FALSE, ncol(st16))
for (rk in colnames(tx16)) org <- org | grepl("Chloroplast|Mitochondri", tx16[[rk]], ignore.case = TRUE)
org[is.na(org)] <- FALSE; names(org) <- colnames(st16)
stI <- readRDS(file.path(P, "microbiome/ITS/seqtab_nochim.rds"))

clr <- function(m) { m <- as.matrix(m)
  if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

go <- function(amp, comp, raw, decon, floor_) {
  f <- file.path(P, "results/Q5", sprintf("q5_%s_%s_objects.rds", amp, comp))
  if (!file.exists(f)) return(invisible())
  o <- readRDS(f); m <- o$manifest; env <- o$env
  ids <- intersect(rownames(raw), m$run_accession)
  k <- match(ids, m$run_accession); m <- m[k, ]; env <- env[k, ]
  cnt <- raw[ids, , drop = FALSE]
  if (decon) cnt <- cnt[, !org[colnames(cnt)], drop = FALSE]
  d <- rowSums(cnt)
  sel <- d >= floor_
  cat(sprintf("\n-- %s / %-7s libraries %3d -> %3d pass depth floor %d  (median depth %.0f)\n",
              amp, comp, length(d), sum(sel), floor_, median(d)))
  if (sum(sel) < 20) { cat("   <20 libraries; not fitted\n"); return(invisible()) }
  cnt <- cnt[sel, , drop = FALSE]; m <- m[sel, ]; env <- env[sel, ]
  cnt <- cnt[, colSums(cnt > 0) >= 5 & colSums(cnt) / sum(cnt) >= 1e-4, drop = FALSE]
  if (ncol(cnt) < 10) { cat("   <10 ASVs survive filter; not fitted\n"); return(invisible()) }
  X <- clr(cnt); D <- dist(X)
  bio <- env[, grep("^bio[0-9]+$", names(env)), drop = FALSE]
  bio <- bio[, apply(bio, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  ok <- complete.cases(bio)
  stand <- droplevels(factor(m$site_code))
  if (nlevels(stand) < 2 || sum(ok) < 20) { cat("   insufficient strata\n"); return(invisible()) }
  a1 <- adonis2(D ~ stand, permutations = 9999)
  pcx <- prcomp(scale(bio[ok, , drop = FALSE]))$x[, 1:3]
  a2 <- adonis2(as.dist(as.matrix(D)[ok, ok]) ~ pcx, permutations = 9999)
  cat(sprintf("   ASVs %4d | stands %d\n", ncol(cnt), nlevels(stand)))
  cat(sprintf("   stand identity : R2 = %.4f  p = %.4f\n", a1$R2[1], a1$`Pr(>F)`[1]))
  cat(sprintf("   climate PC1-3  : R2 = %.4f  p = %.4f\n", a2$R2[1], a2$`Pr(>F)`[1]))
  cat(sprintf("   ratio stand/climate = %.2f\n", a1$R2[1] / a2$R2[1]))
}

cat("=========== 16S, ORGANELLAR REMOVED, depth floor 100 bacterial reads ===========\n")
for (cc in c("needle", "branch", "root")) go("16S", cc, st16, TRUE, 100)
cat("\n\n=========== ITS, prevalence>=5 & rel>=1e-4 ===========\n")
for (cc in c("needle", "branch", "root")) go("ITS", cc, stI, FALSE, 1000)
sink()
