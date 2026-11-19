#!/usr/bin/env Rscript
# Complete prevalence gradient with the presence/absence twin at EVERY threshold, so the
# published figure shows the full decay rather than two endpoints. Above-ground pool, floor 250.
suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
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
clr <- function(m) { m <- as.matrix(m); if (any(m == 0)) m[m == 0] <- 0.65*min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x)))) }
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
fit <- function(Y, tg, mg, j) {
  Ks <- K[tg,tg]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm=TRUE); diag(Ks) <- max(Ks, na.rm=TRUE)
  pc <- cmdscale(as.dist(max(Ks)-Ks), k=4); colnames(pc) <- paste0("G",1:4)
  mem <- pcnm(dist(cbind(ch$lon[j], ch$lat[j])))$vectors[,1:4,drop=FALSE]; colnames(mem) <- paste0("S",1:4)
  dat <- data.frame(pc, mem, stand=droplevels(factor(mg$site_code)), age=factor(mg$age_class),
                    logdepth=log(rowSums(exp(Y))))
  cap <- try(dbrda(dist(Y) ~ G1+G2+G3+G4 + Condition(stand+logdepth+S1+S2+S3+S4+age), data=dat), silent=TRUE)
  if (inherits(cap,"try-error") || is.null(cap$CCA)) return(NULL)
  a <- anova(cap, permutations=9999)
  list(r2=cap$CCA$tot.chi/cap$tot.chi, F=a$F[1], p=a$`Pr(>F)`[1])
}
cat("amp\tprev\tasv\tn\tclr_r2\tclr_F\tclr_p\tpa_p\n")
for (amp in c("16S","ITS")) {
  b <- pool(amp); sel <- rowSums(b$cnt) >= 250
  cnt <- b$cnt[sel,,drop=FALSE]; tr <- b$trees[sel]
  g <- tr %in% rownames(K); cnt <- cnt[g,,drop=FALSE]; tg <- tr[g]
  mg <- man[match(tg, man$tree_id),]; j <- match(tg, ch$tree_id)
  for (pv in c(1,2,3,5,10,20)) {
    C <- cnt[, colSums(cnt > 0) >= pv, drop=FALSE]
    if (ncol(C) < 10) next
    rc <- fit(clr(C), tg, mg, j)
    PA <- (C > 0)*1; PA <- PA[, colSums(PA) > 0 & colSums(PA) < nrow(PA), drop=FALSE]
    rp <- if (ncol(PA) >= 10) fit(clr(PA + 0.5), tg, mg, j) else NULL
    cat(sprintf("%s\t%d\t%d\t%d\t%.4f\t%.3f\t%.4f\t%s\n", amp, pv, ncol(C), nrow(C),
        rc$r2, rc$F, rc$p, if (is.null(rp)) "NA" else sprintf("%.4f", rp$p)))
  }
}
