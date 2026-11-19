#!/usr/bin/env Rscript
# Above-ground pool (branch + needle per tree), floor 250 reads, NO prevalence/abundance filter.
# Only all-zero columns are dropped, which is mandatory for CLR and is not a prevalence filter.
#
# THE PRESENCE/ABSENCE CONTROL RUNS ALONGSIDE, NOT AFTERWARDS. Removing the prevalence filter
# restores exactly the condition that produced the earlier needle artifact: at prevalence>=2 the
# old 16S needle model gave p = 0.013, and 83% of it was reproducible from presence/absence
# alone, i.e. from WHICH ASVs happened to be seen at low depth rather than their abundances.
# So any significant result here is reported next to its presence/absence twin:
#   CLR >> presence/absence  -> the signal is in abundances, credible
#   CLR ~= presence/absence  -> the signal is in the zero pattern, the artifact again
# A prevalence gradient is printed too, so the threshold's effect is visible rather than chosen.

suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q5/ag_250_nofilter.txt"), split = TRUE)

man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
rad <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
ch  <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)
relf <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
to <- setNames(rad$tree_id, rad$run_accession)
ta <- to[bare(as.character(relf$ida))]; tb <- to[bare(as.character(relf$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta != tb & relf$nSites >= 2000 & is.finite(relf$theta)
kt <- sort(unique(c(ta[gd], tb[gd])))
K <- matrix(NA_real_, length(kt), length(kt), dimnames = list(kt, kt))
for (z in which(gd)) { K[ta[z], tb[z]] <- relf$theta[z]; K[tb[z], ta[z]] <- relf$theta[z] }
diag(K) <- max(relf$theta[gd], na.rm = TRUE)

clr <- function(m) { m <- as.matrix(m); if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

build <- function(amp) {
  st <- readRDS(file.path(P, sprintf("microbiome/%s/seqtab_nochim.rds", amp)))
  if (amp == "16S") {
    tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
    o <- rep(FALSE, ncol(st))
    for (rk in colnames(tx)) o <- o | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case = TRUE)
    o[is.na(o)] <- FALSE; st <- st[, !o, drop = FALSE]
  }
  i <- match(rownames(st), man$run_accession)
  keep <- man$compartment[i] %in% c("branch", "needle") & !is.na(man$tree_id[i])
  st <- st[keep, , drop = FALSE]; tid <- man$tree_id[i][keep]
  trees <- sort(unique(tid))
  list(cnt = t(sapply(trees, function(x) colSums(st[tid == x, , drop = FALSE]))), trees = trees)
}

fitQ6 <- function(Y, tr, tag, nperm = 9999) {
  g <- tr %in% rownames(K)
  if (sum(g) < 25) { cat(sprintf("   %-26s %d genotyped trees - not fitted\n", tag, sum(g))); return(invisible()) }
  Yg <- Y[g, , drop = FALSE]; tg <- tr[g]
  mg <- man[match(tg, man$tree_id), ]; j <- match(tg, ch$tree_id)
  Ks <- K[tg, tg]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm = TRUE); diag(Ks) <- max(Ks, na.rm = TRUE)
  pc <- cmdscale(as.dist(max(Ks) - Ks), k = 4); colnames(pc) <- paste0("G", 1:4)
  mem <- pcnm(dist(cbind(ch$lon[j], ch$lat[j])))$vectors[, 1:4, drop = FALSE]; colnames(mem) <- paste0("S", 1:4)
  dat <- data.frame(pc, mem, stand = droplevels(factor(mg$site_code)),
                    age = factor(mg$age_class), logdepth = log(rowSums(exp(Yg))))
  cap <- try(dbrda(dist(Yg) ~ G1+G2+G3+G4 + Condition(stand + logdepth + S1+S2+S3+S4 + age),
                   data = dat), silent = TRUE)
  if (inherits(cap, "try-error") || is.null(cap$CCA)) { cat(sprintf("   %-26s model failed\n", tag)); return(invisible()) }
  a <- anova(cap, permutations = nperm)
  cat(sprintf("   %-26s n=%2d ASVs=%4d  R2=%.4f  F=%.3f  p=%.4f %s\n", tag, nrow(Yg), ncol(Yg),
              cap$CCA$tot.chi/cap$tot.chi, a$F[1], a$`Pr(>F)`[1],
              if (a$`Pr(>F)`[1] < 0.05) "<< DETECTED" else ""))
  invisible(a$`Pr(>F)`[1])
}

for (amp in c("16S", "ITS")) {
  b <- build(amp); d <- rowSums(b$cnt)
  sel <- d >= 250
  cnt <- b$cnt[sel, , drop = FALSE]; tr <- b$trees[sel]
  cnt <- cnt[, colSums(cnt) > 0, drop = FALSE]
  cat(sprintf("\n############ %s above-ground, >=250 reads, NO prevalence filter ############\n", amp))
  cat(sprintf("trees %d (of %d) | ASVs %d | median depth %.0f | sparsity %.1f%%\n",
              nrow(cnt), length(d), ncol(cnt), median(rowSums(cnt)), 100*mean(cnt == 0)))
  Y <- clr(cnt); D <- dist(Y)
  mg <- man[match(tr, man$tree_id), ]; j <- match(tr, ch$tree_id)
  stand <- droplevels(factor(mg$site_code))
  bio <- ch[j, grep("^bio[0-9]+$", names(ch)), drop = FALSE]
  bio <- bio[, apply(bio, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  ok <- complete.cases(bio)
  if (nlevels(stand) >= 2 && sum(ok) >= 20) {
    a1 <- adonis2(D ~ stand, permutations = 9999)
    a2 <- adonis2(as.dist(as.matrix(D)[ok, ok]) ~ prcomp(scale(bio[ok, ]))$x[, 1:3], permutations = 9999)
    cat(sprintf("[Q5] stand R2=%.4f (p=%.4f) | climate R2=%.4f (p=%.4f) | ratio %.2f\n",
                a1$R2[1], a1$`Pr(>F)`[1], a2$R2[1], a2$`Pr(>F)`[1], a1$R2[1]/a2$R2[1]))
  }
  cat("[Q6] host genetics:\n")
  fitQ6(Y, tr, "CLR (abundance)")
  PA <- (cnt > 0) * 1; PA <- PA[, colSums(PA) > 0 & colSums(PA) < nrow(PA), drop = FALSE]
  fitQ6(clr(PA + 0.5), tr, "presence/absence CONTROL")
  cat("[Q6] prevalence gradient:\n")
  for (pv in c(2, 3, 5, 10)) {
    C2 <- cnt[, colSums(cnt > 0) >= pv, drop = FALSE]
    if (ncol(C2) >= 10) fitQ6(clr(C2), tr, sprintf("prevalence >=%d", pv))
  }
}
sink()
