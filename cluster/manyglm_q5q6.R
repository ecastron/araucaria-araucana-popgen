#!/usr/bin/env Rscript
# Q5 and Q6 by MODEL-BASED multivariate GLM (mvabund::manyglm), as a check on the
# distance-based results.
#
# WHY. PERMANOVA tests centroids on a distance matrix and is known to CONFOUND LOCATION WITH
# DISPERSION: groups differing only in within-group spread can come out significant. That is not
# hypothetical here -- beta dispersion is unequal in this dataset (p = 0.031). manyglm instead
# fits a negative-binomial GLM per taxon with a library-size offset and combines them into one
# multivariate statistic by resampling, so the mean-variance relationship is modelled rather than
# absorbed into a distance. Counts enter raw; no CLR, no zero replacement.
#
# betadisper is run alongside to say explicitly how much of the stand effect is dispersion.
#
# HOW THE TWO ARE COMPARED. manyglm returns a summed likelihood-ratio, not an R2, so there is no
# direct equivalent of the "stand explains 2.1x climate" ratio. Deviance explained is reported
# instead. The conclusion is only claimed where BOTH methods agree; disagreement is reported as
# disagreement rather than resolved in favour of whichever is more convenient.

suppressPackageStartupMessages({library(mvabund); library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"; NB <- 499
sink(file.path(P, "results/Q5/manyglm.txt"), split = TRUE)
cat(sprintf("mvabund %s | resampling: PIT-trap, %d iterations\n", packageVersion("mvabund"), NB))

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
pool <- function(amp) {
  st <- readRDS(file.path(P, sprintf("microbiome/%s/seqtab_nochim.rds", amp)))
  if (amp == "16S") { tx <- as.data.frame(readRDS(file.path(P,"microbiome/16S/taxonomy.rds")), stringsAsFactors=FALSE)
    o <- rep(FALSE, ncol(st)); for (rk in colnames(tx)) o <- o | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case=TRUE)
    o[is.na(o)] <- FALSE; st <- st[, !o, drop=FALSE] }
  i <- match(rownames(st), man$run_accession)
  k <- man$compartment[i] %in% c("branch","needle") & !is.na(man$tree_id[i])
  st <- st[k,,drop=FALSE]; tid <- man$tree_id[i][k]; tr <- sort(unique(tid))
  list(cnt = t(sapply(tr, function(x) colSums(st[tid==x,,drop=FALSE]))), trees = tr)
}
devexp <- function(m) { d <- sum(m$deviance); n0 <- sum(m$null.deviance); 100*(n0-d)/n0 }

for (amp in c("16S","ITS")) {
  b <- pool(amp); sel <- rowSums(b$cnt) >= 250
  cnt <- b$cnt[sel,,drop=FALSE]; tr <- b$trees[sel]
  cnt <- cnt[, colSums(cnt>0) >= 5 & colSums(cnt)/sum(cnt) >= 1e-4, drop=FALSE]
  mg <- man[match(tr, man$tree_id),]; j <- match(tr, ch$tree_id)
  off <- log(rowSums(cnt))
  stand <- droplevels(factor(mg$site_code))
  bio <- ch[j, grep("^bio[0-9]+$", names(ch)), drop=FALSE]
  bio <- bio[, apply(bio,2,function(x) sd(x,na.rm=TRUE)>0), drop=FALSE]
  ok <- complete.cases(bio)
  pcx <- matrix(NA, nrow(cnt), 3); pcx[ok,] <- prcomp(scale(bio[ok,]))$x[,1:3]
  cat(sprintf("\n################ %s above-ground: %d trees x %d ASVs ################\n",
              amp, nrow(cnt), ncol(cnt)))
  mv <- mvabund(cnt)

  cat("\n--- Q5: stand vs climate, negative-binomial manyglm ---\n")
  m1 <- manyglm(mv ~ stand, family="negative.binomial", offset=off)
  a1 <- anova(m1, nBoot=NB, resamp="pit.trap", show.time="none")
  cat(sprintf("  stand   : deviance explained %.1f%%  LR = %.1f  p = %.4f\n",
              devexp(m1), a1$table[2,3], a1$table[2,4]))
  m2 <- manyglm(mv[ok,] ~ pcx[ok,], family="negative.binomial", offset=off[ok])
  a2 <- anova(m2, nBoot=NB, resamp="pit.trap", show.time="none")
  cat(sprintf("  climate : deviance explained %.1f%%  LR = %.1f  p = %.4f\n",
              devexp(m2), a2$table[2,3], a2$table[2,4]))
  cat(sprintf("  deviance ratio stand/climate = %.2f\n", devexp(m1)/devexp(m2)))

  cat("\n--- how much of the stand effect is DISPERSION, not location? ---\n")
  Dd <- vegdist(decostand(cnt, "hellinger"), "euclidean")
  bd <- betadisper(Dd, stand)
  pb <- permutest(bd, permutations=999)$tab[1,6]
  cat(sprintf("  betadisper across stands: F = %.3f  p = %.4f  -> %s\n",
      permutest(bd, permutations=999)$tab[1,4], pb,
      if (pb < 0.05) "spread DIFFERS between stands; PERMANOVA partly reflects dispersion" else
                     "spread is homogeneous; the PERMANOVA stand effect is a location effect"))

  cat("\n--- Q6: host genetics, sequential manyglm (genetics entered LAST) ---\n")
  g <- tr %in% rownames(K0)
  if (sum(g) < 25) { cat(sprintf("  %d genotyped trees - not fitted\n", sum(g))); next }
  cg <- cnt[g,,drop=FALSE]; cg <- cg[, colSums(cg>0) >= 5, drop=FALSE]
  tg <- tr[g]; mgg <- man[match(tg, man$tree_id),]; jg <- match(tg, ch$tree_id)
  K <- K0[tg,tg]; K[is.na(K)] <- mean(K[upper.tri(K)],na.rm=TRUE); diag(K) <- max(K,na.rm=TRUE)
  K <- (K+t(K))/2
  G <- cmdscale(as.dist(max(K)-K), k=4); colnames(G) <- paste0("G",1:4)
  S <- pcnm(dist(cbind(ch$lon[jg], ch$lat[jg])))$vectors[,1:4,drop=FALSE]; colnames(S) <- paste0("S",1:4)
  st2 <- droplevels(factor(mgg$site_code)); ag2 <- factor(mgg$age_class)
  mv2 <- mvabund(cg); off2 <- log(rowSums(cg))
  m3 <- manyglm(mv2 ~ st2 + ag2 + S + G, family="negative.binomial", offset=off2)
  a3 <- anova(m3, nBoot=NB, resamp="pit.trap", show.time="none")
  print(a3$table)
  gi <- grep("^G$|G$", rownames(a3$table))
  gi <- gi[length(gi)]
  cat(sprintf("\n  HOST GENETICS (entered last): LR = %.1f  p = %.4f  -> %s\n",
      a3$table[gi,3], a3$table[gi,4],
      if (!is.na(a3$table[gi,4]) && a3$table[gi,4] < 0.05) "DETECTED" else "not detected"))
}
sink()
