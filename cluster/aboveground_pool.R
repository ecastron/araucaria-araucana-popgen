#!/usr/bin/env Rscript
# POOL branch + needle per tree into an "above-ground" compartment.
#
# Counts are summed across all branch and needle libraries belonging to one tree. Summing raw
# counts is the correct merge for compositional count data; averaging proportions would give a
# 27-read needle library the same weight as a 27,000-read branch library.
#
# THE TEST THAT DECIDES WHETHER THIS IS WORTH DOING, declared before running: does pooling
# recover trees that BRANCH ALONE does not? Needle holds a median of 27 bacterial reads against
# branch's 130 (and branch's 90th percentile is 14,556), so pooling may simply be branch plus
# rounding error. If usable-tree counts are the same either way, the "above-ground" framing adds
# a compartment label without adding data, and should be dropped rather than reported.

suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q5/aboveground_pool.txt"), split = TRUE)

st <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org <- rep(FALSE, ncol(st))
for (rk in colnames(tx)) org <- org | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case = TRUE)
org[is.na(org)] <- FALSE
bac <- st[, !org, drop = FALSE]
man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
i <- match(rownames(bac), man$run_accession)
comp <- man$compartment[i]; tid <- man$tree_id[i]
ag <- which(comp %in% c("branch", "needle") & !is.na(tid))

cat(sprintf("above-ground libraries: %d (branch %d, needle %d) across %d trees\n",
            length(ag), sum(comp[ag] == "branch"), sum(comp[ag] == "needle"),
            length(unique(tid[ag]))))
tt <- table(tid[ag], comp[ag])
cat(sprintf("trees with BOTH tissues: %d | branch only: %d | needle only: %d\n",
            sum(tt[, "branch"] > 0 & tt[, "needle"] > 0),
            sum(tt[, "branch"] > 0 & tt[, "needle"] == 0),
            sum(tt[, "branch"] == 0 & tt[, "needle"] > 0)))

trees <- sort(unique(tid[ag]))
pool <- t(sapply(trees, function(x) colSums(bac[ag[tid[ag] == x], , drop = FALSE])))
bo   <- t(sapply(trees, function(x) { k <- ag[tid[ag] == x & comp[ag] == "branch"]
                 if (!length(k)) rep(0, ncol(bac)) else colSums(bac[k, , drop = FALSE]) }))
dp <- rowSums(pool); db <- rowSums(bo)

cat("\n--- pooled above-ground depth per tree ---\n")
q <- quantile(dp, c(0, .1, .25, .5, .75, .9, 1))
cat(sprintf("  min %.0f | 10%% %.0f | 25%% %.0f | median %.0f | 75%% %.0f | 90%% %.0f | max %.0f\n",
            q[1], q[2], q[3], q[4], q[5], q[6], q[7]))

rad <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
geno <- unique(rad$tree_id)
cat("\n--- DOES POOLING RECOVER TREES THAT BRANCH ALONE DOES NOT? ---\n")
cat(sprintf("  %-10s %-22s %-22s %s\n", "floor", "branch alone", "pooled above-ground", "gain"))
for (fl in c(100, 250, 500, 1000, 2000, 5000)) {
  a <- sum(db >= fl); b <- sum(dp >= fl)
  ag_ <- sum(db >= fl & trees %in% geno); bg <- sum(dp >= fl & trees %in% geno)
  cat(sprintf("  >=%-8d %3d trees (%2d geno)      %3d trees (%2d geno)      %+d trees (%+d geno)\n",
              fl, a, ag_, b, bg, b - a, bg - ag_))
}
cat(sprintf("\n  needle's median contribution to a pooled tree: %.1f%% of its reads\n",
            100 * median((dp - db) / pmax(dp, 1))))

# ---- if it clears 25 genotyped trees, actually fit Q5 and Q6 -------------------------------
FL <- 500
sel <- dp >= FL
cat(sprintf("\n\n================ above-ground at >=%d bacterial reads ================\n", FL))
cat(sprintf("trees: %d | genotyped: %d\n", sum(sel), sum(sel & trees %in% geno)))
X <- pool[sel, , drop = FALSE]; tr <- trees[sel]
X <- X[, colSums(X > 0) >= 5 & colSums(X) / sum(X) >= 1e-4, drop = FALSE]
cat(sprintf("ASVs after prevalence>=5 & rel>=1e-4: %d\n", ncol(X)))
ch <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)
j <- match(tr, ch$tree_id)
meta <- man[match(tr, man$tree_id), ]
clr <- function(m) { m <- as.matrix(m); if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }
Y <- clr(X); D <- dist(Y)
stand <- droplevels(factor(meta$site_code))
bio <- ch[j, grep("^bio[0-9]+$", names(ch)), drop = FALSE]
bio <- bio[, apply(bio, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
ok <- complete.cases(bio)
if (nlevels(stand) >= 2 && sum(ok) >= 20) {
  a1 <- adonis2(D ~ stand, permutations = 9999)
  a2 <- adonis2(as.dist(as.matrix(D)[ok, ok]) ~ prcomp(scale(bio[ok, ]))$x[, 1:3], permutations = 9999)
  cat(sprintf("[Q5] stand R2 = %.4f (p=%.4f) | climate R2 = %.4f (p=%.4f) | ratio %.2f\n",
              a1$R2[1], a1$`Pr(>F)`[1], a2$R2[1], a2$`Pr(>F)`[1], a1$R2[1]/a2$R2[1]))
}
rel_ <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
to <- setNames(rad$tree_id, rad$run_accession)
ta <- to[bare(as.character(rel_$ida))]; tb <- to[bare(as.character(rel_$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta != tb & rel_$nSites >= 2000 & is.finite(rel_$theta)
kt <- sort(unique(c(ta[gd], tb[gd])))
K <- matrix(NA_real_, length(kt), length(kt), dimnames = list(kt, kt))
for (z in which(gd)) { K[ta[z], tb[z]] <- rel_$theta[z]; K[tb[z], ta[z]] <- rel_$theta[z] }
diag(K) <- max(rel_$theta[gd], na.rm = TRUE)
g <- tr %in% rownames(K)
cat(sprintf("[Q6] genotyped trees available: %d\n", sum(g)))
if (sum(g) >= 25) {
  Yg <- Y[g, , drop = FALSE]; mg <- meta[g, ]; jg <- j[g]
  Ks <- K[tr[g], tr[g]]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm = TRUE)
  diag(Ks) <- max(Ks, na.rm = TRUE)
  pc <- cmdscale(as.dist(max(Ks) - Ks), k = 4); colnames(pc) <- paste0("G", 1:4)
  mem <- pcnm(dist(cbind(ch$lon[jg], ch$lat[jg])))$vectors[, 1:4, drop = FALSE]
  colnames(mem) <- paste0("S", 1:4)
  dat <- data.frame(pc, mem, stand = droplevels(factor(mg$site_code)),
                    age = factor(mg$age_class), logdepth = log(rowSums(exp(Yg))))
  cap <- try(dbrda(dist(Yg) ~ G1+G2+G3+G4 + Condition(stand + logdepth + S1+S2+S3+S4 + age),
                   data = dat), silent = TRUE)
  if (!inherits(cap, "try-error") && !is.null(cap$CCA)) {
    a <- anova(cap, permutations = 9999)
    cat(sprintf("[Q6] n=%d  host genetics partial R2 = %.4f  F = %.3f  p = %.4f  %s\n",
                nrow(mg), cap$CCA$tot.chi/cap$tot.chi, a$F[1], a$`Pr(>F)`[1],
                if (a$`Pr(>F)`[1] < 0.05) "<< DETECTED" else "not detected"))
  } else cat("[Q6] model failed\n")
} else cat("[Q6] below 25 genotyped trees, NOT FITTED\n")
sink()
