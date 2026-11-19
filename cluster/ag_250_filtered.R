#!/usr/bin/env Rscript
# FINAL SETTING: above-ground pool (branch + needle per tree), floor 250 reads,
# prevalence >=5 trees AND relative abundance >=1e-4.
#
# Chosen because at >=250 with NO prevalence filter the 16S POSITIVE CONTROL failed
# (stand p = 0.302): 1,462 columns at 96.7% sparsity swamped a stand effect that is
# real. The filter is retained for that reason, not for tidiness.
#
# Every model is reported beside its PRESENCE/ABSENCE twin. If discarding abundances
# does not weaken a result, that result is a zero-pattern artifact rather than
# community structure, and must not be reported as the latter.
# 16S root is carried along as an uncontaminated reference (68% host, median 9,618 reads).

suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q5/ag_250_filtered.txt"), split = TRUE)

man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
rad <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
ch  <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)
rf  <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
to <- setNames(rad$tree_id, rad$run_accession)
ta <- to[bare(as.character(rf$ida))]; tb <- to[bare(as.character(rf$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta != tb & rf$nSites >= 2000 & is.finite(rf$theta)
kt <- sort(unique(c(ta[gd], tb[gd])))
K <- matrix(NA_real_, length(kt), length(kt), dimnames = list(kt, kt))
for (z in which(gd)) { K[ta[z], tb[z]] <- rf$theta[z]; K[tb[z], ta[z]] <- rf$theta[z] }
diag(K) <- max(rf$theta[gd], na.rm = TRUE)

clr <- function(m) { m <- as.matrix(m); if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

pool <- function(amp, comps) {
  st <- readRDS(file.path(P, sprintf("microbiome/%s/seqtab_nochim.rds", amp)))
  if (amp == "16S") {
    tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
    o <- rep(FALSE, ncol(st))
    for (rk in colnames(tx)) o <- o | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case = TRUE)
    o[is.na(o)] <- FALSE; st <- st[, !o, drop = FALSE]
  }
  i <- match(rownames(st), man$run_accession)
  k <- man$compartment[i] %in% comps & !is.na(man$tree_id[i])
  st <- st[k, , drop = FALSE]; tid <- man$tree_id[i][k]
  tr <- sort(unique(tid))
  list(cnt = t(sapply(tr, function(x) colSums(st[tid == x, , drop = FALSE]))), trees = tr)
}

q6 <- function(Y, tr, tag) {
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
  if (inherits(cap, "try-error") || is.null(cap$CCA)) { cat(sprintf("   %-26s failed\n", tag)); return(invisible()) }
  a <- anova(cap, permutations = 9999)
  cat(sprintf("   %-26s n=%2d ASVs=%4d  R2=%.4f  F=%.3f  p=%.4f %s\n", tag, nrow(Yg), ncol(Yg),
              cap$CCA$tot.chi/cap$tot.chi, a$F[1], a$`Pr(>F)`[1], if (a$`Pr(>F)`[1] < 0.05) "<<" else ""))
}

go <- function(amp, comps, lab, floor_ = 250) {
  b <- pool(amp, comps); d <- rowSums(b$cnt); sel <- d >= floor_
  cnt <- b$cnt[sel, , drop = FALSE]; tr <- b$trees[sel]
  cnt <- cnt[, colSums(cnt > 0) >= 5 & colSums(cnt)/sum(cnt) >= 1e-4, drop = FALSE]
  cat(sprintf("\n############ %s ############\n", lab))
  cat(sprintf("trees %d of %d | ASVs %d | median depth %.0f | sparsity %.1f%%\n",
              nrow(cnt), length(d), ncol(cnt), median(rowSums(cnt)), 100*mean(cnt == 0)))
  if (nrow(cnt) < 20 || ncol(cnt) < 10) { cat("EXCLUDED\n"); return(invisible()) }
  Y <- clr(cnt); D <- dist(Y)
  mg <- man[match(tr, man$tree_id), ]; j <- match(tr, ch$tree_id)
  stand <- droplevels(factor(mg$site_code))
  bio <- ch[j, grep("^bio[0-9]+$", names(ch)), drop = FALSE]
  bio <- bio[, apply(bio, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  ok <- complete.cases(bio)
  if (nlevels(stand) >= 2 && sum(ok) >= 20) {
    a1 <- adonis2(D ~ stand, permutations = 9999)
    a2 <- adonis2(as.dist(as.matrix(D)[ok, ok]) ~ prcomp(scale(bio[ok, ]))$x[, 1:3], permutations = 9999)
    cat(sprintf("[Q5] stand R2=%.4f (p=%.4f) | climate R2=%.4f (p=%.4f) | ratio %.2f  %s\n",
                a1$R2[1], a1$`Pr(>F)`[1], a2$R2[1], a2$`Pr(>F)`[1], a1$R2[1]/a2$R2[1],
                if (a1$`Pr(>F)`[1] < 0.05) "CONTROL OK" else "CONTROL FAILED"))
  }
  cat("[Q6]\n")
  q6(Y, tr, "CLR (abundance)")
  PA <- (cnt > 0)*1; PA <- PA[, colSums(PA) > 0 & colSums(PA) < nrow(PA), drop = FALSE]
  if (ncol(PA) >= 10) q6(clr(PA + 0.5), tr, "presence/absence CONTROL")
}

go("16S", c("branch","needle"), "16S above-ground, >=250 reads, prev>=5 & rel>=1e-4")
go("ITS", c("branch","needle"), "ITS above-ground, >=250 reads, prev>=5 & rel>=1e-4")
go("16S", "root", "16S root (reference), >=250 reads, prev>=5 & rel>=1e-4")
go("ITS", "root", "ITS root (reference), >=250 reads, prev>=5 & rel>=1e-4")
sink()
