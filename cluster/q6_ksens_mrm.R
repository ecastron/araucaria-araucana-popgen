#!/usr/bin/env Rscript
# TWO CHECKS ON THE Q6 NULL, both aimed at the objection that the result is an artefact of
# compressing 78,802 SNPs into four kinship PCoA axes.
#
# 1. AXIS SENSITIVITY. Refit the same dbRDA at k = 2, 4, 8, 16, 24 kinship axes. If the null is
#    caused by compression it should weaken as k rises. If it holds at every k, compression is
#    not the explanation. Note that raising k also spends degrees of freedom, so R2 will drift
#    upward mechanically -- F and p are the statistics to read, not R2.
#
# 2. MRM as a SECONDARY check only. Multiple regression on distance matrices is Mantel-family and
#    inherits Mantel's inflated type I error under spatial autocorrelation, which Q2 established
#    we have. It is therefore NOT the primary test and is not reported as one. It is run because
#    agreement between two methods with different failure modes is worth more than either alone;
#    a DISagreement would favour dbRDA, whose assumptions fit this design better.

suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q6/ksens_mrm.txt"), split = TRUE)

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
hav <- function(la1,lo1,la2,lo2){R<-6371000;p1<-la1*pi/180;p2<-la2*pi/180
  2*R*asin(sqrt(sin((p2-p1)/2)^2+cos(p1)*cos(p2)*sin((lo2-lo1)*pi/180/2)^2))}

for (amp in c("16S","ITS")) {
  b <- pool(amp); sel <- rowSums(b$cnt) >= 250
  cnt <- b$cnt[sel,,drop=FALSE]; tr <- b$trees[sel]
  g <- tr %in% rownames(K0); cnt <- cnt[g,,drop=FALSE]; tg <- tr[g]
  cnt <- cnt[, colSums(cnt>0) >= 5 & colSums(cnt)/sum(cnt) >= 1e-4, drop=FALSE]
  n <- nrow(cnt)
  cat(sprintf("\n############ %s above-ground: %d trees x %d ASVs ############\n", amp, n, ncol(cnt)))
  Y <- clr(cnt); D <- dist(Y)
  K <- K0[tg,tg]; K[is.na(K)] <- mean(K[upper.tri(K)], na.rm=TRUE); diag(K) <- max(K, na.rm=TRUE)
  K <- (K+t(K))/2
  mg <- man[match(tg, man$tree_id),]; j <- match(tg, ch$tree_id)
  mem <- pcnm(dist(cbind(ch$lon[j], ch$lat[j])))$vectors[,1:4,drop=FALSE]; colnames(mem) <- paste0("S",1:4)
  Kd <- as.dist(max(K)-K)
  full <- cmdscale(Kd, k = n-2, eig = TRUE)
  pos <- full$eig[full$eig > 0]
  cat(sprintf("kinship PCoA: %d positive eigenvalues\n", length(pos)))
  cat(sprintf("  %-6s %-8s %9s %8s %9s %s\n","k","%eigen","partialR2","F","p",""))
  for (k in c(2,4,8,16,24)) {
    if (k > n-3) next
    pcs <- full$points[,1:k,drop=FALSE]; colnames(pcs) <- paste0("G",1:k)
    dat <- data.frame(pcs, mem, stand=droplevels(factor(mg$site_code)),
                      age=factor(mg$age_class), logdepth=log(rowSums(exp(Y))))
    fo <- as.formula(paste("D ~", paste(colnames(pcs), collapse="+"),
          "+ Condition(stand+logdepth+S1+S2+S3+S4+age)"))
    cap <- try(dbrda(fo, data=dat), silent=TRUE)
    if (inherits(cap,"try-error") || is.null(cap$CCA)) { cat(sprintf("  %-6d model failed\n", k)); next }
    rdf <- cap$CA$rank
    if (is.null(rdf) || rdf < 1) {
      cat(sprintf("  %-6d %-8.1f %9s %8s %9s  SATURATED: no residual df at n=%d\n", k,
          100*sum(pos[1:min(k,length(pos))])/sum(pos), "-", "-", "-", n)); next }
    a <- anova(cap, permutations=9999)
    pp <- a$`Pr(>F)`[1]
    cat(sprintf("  %-6d %-8.1f %9.4f %8.3f %9s %s (resid df %d)\n", k,
        100*sum(pos[1:min(k,length(pos))])/sum(pos), cap$CCA$tot.chi/cap$tot.chi, a$F[1],
        if (is.na(pp)) "NA" else sprintf("%.4f", pp),
        if (!is.na(pp) && pp<0.05) "<<" else "", rdf))
  }
  # ---- MRM, secondary ----
  dm <- as.vector(D)
  dk <- as.vector(as.dist(max(K)-K))
  dg <- as.vector(as.dist(outer(seq_len(n), seq_len(n), Vectorize(function(a,b)
        hav(ch$lat[j][a], ch$lon[j][a], ch$lat[j][b], ch$lon[j][b])))))
  ds <- as.vector(dist(model.matrix(~0+factor(mg$site_code))))
  fit <- lm(dm ~ dk + log1p(dg) + ds)
  co <- summary(fit)$coefficients
  B <- 9999; obs <- abs(co["dk","t value"]); cnt_ <- 0
  M <- as.matrix(D)
  for (z in seq_len(B)) {
    o <- sample(n); dmp <- as.vector(as.dist(M[o,o]))
    tt <- summary(lm(dmp ~ dk + log1p(dg) + ds))$coefficients["dk","t value"]
    if (abs(tt) >= obs) cnt_ <- cnt_ + 1
  }
  cat(sprintf("\n  [MRM secondary] kinship coef = %+.4f  t = %+.3f  permutation p = %.4f  %s\n",
      co["dk","Estimate"], co["dk","t value"], (1+cnt_)/(1+B),
      if ((1+cnt_)/(1+B) < 0.05) "<< significant" else "not significant"))
  cat(sprintf("  (geography and stand held in the same model; MRM is Mantel-family and secondary)\n"))
}
sink()
