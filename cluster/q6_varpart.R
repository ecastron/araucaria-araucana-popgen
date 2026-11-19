#!/usr/bin/env Rscript
# Extended-phenotype test by variance partitioning: CLIMATE | GEOGRAPHY | KINSHIP
#
# Community prep is byte-identical to q6_rerun_clean.R (decontaminated, above-ground pooled,
# 250-read floor, prevalence >=5 trees & rel >=1e-4). What differs: instead of conditioning
# stand out, the three explanatory tables compete directly.
#
# Climate = bio7 (temperature annual range), bio13 (precip wettest month), elevation.
# bio14 dropped: r = 0.875 with bio13 across the 12 stands, VIF 17-19.
# NOTE the four CHELSA/terrain variables are effectively STAND-LEVEL (ICC > 0.95 for all 19
# bioclim), so the climate fraction is bounded by the among-stand fraction by construction.
suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
sink(file.path(P,"results/Q6/varpart/q6_varpart.txt"), split=TRUE)

st16 <- readRDS(file.path(P,"microbiome/16S/seqtab_nochim.rds"))
tx16 <- as.data.frame(readRDS(file.path(P,"microbiome/16S/taxonomy.rds")), stringsAsFactors=FALSE)
org <- rep(FALSE, ncol(st16))
for (rk in colnames(tx16)) org <- org | grepl("Chloroplast|Mitochondri", tx16[[rk]], ignore.case=TRUE)
org[is.na(org)] <- FALSE; names(org) <- colnames(st16)
stI <- readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))

rel <- read.delim(file.path(P,"results/Q2/ngsrelate.tsv"), check.names=FALSE)
radman <- read.csv(file.path(P,"metadata/manifest rad runs.csv"), stringsAsFactors=FALSE)
bare <- function(v) sub("\\.[12]$","",v)
tree_of <- setNames(radman$tree_id, radman$run_accession)
ta <- tree_of[bare(as.character(rel$ida))]; tb <- tree_of[bare(as.character(rel$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta!=tb & rel$nSites>=2000 & is.finite(rel$theta)
trees <- sort(unique(c(ta[gd], tb[gd])))
K <- matrix(NA_real_, length(trees), length(trees), dimnames=list(trees,trees))
for (i in which(gd)) { K[ta[i],tb[i]] <- rel$theta[i]; K[tb[i],ta[i]] <- rel$theta[i] }
diag(K) <- max(rel$theta[gd], na.rm=TRUE)

ch <- read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"), check.names=FALSE)
te <- read.delim(file.path(P,"metadata/terrain_per_tree.tsv"), check.names=FALSE)
clim <- merge(ch[,c("tree_id","site_code","lon","lat","bio7","bio13")],
              te[,c("tree_id","dem_elev_m")], by="tree_id", all.x=TRUE)

clr <- function(m){ m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0])
  t(apply(m,1,function(x) log(x)-mean(log(x)))) }

align <- function(raw, comp) {
  o <- readRDS(file.path(P,"results/Q5",sprintf("q5_%s_%s_objects.rds",
       if (identical(raw,st16)) "16S" else "ITS", comp)))
  m <- o$manifest
  keep <- !is.na(m$tree_id) & m$tree_id %in% rownames(K)
  m <- m[keep,]; m <- m[!duplicated(m$tree_id),]
  ids <- intersect(rownames(raw), m$run_accession)
  m <- m[match(ids,m$run_accession),]
  list(cnt=raw[ids,,drop=FALSE], m=m)
}

run <- function(raw, comp, amp, drop_org=FALSE) {
  a <- tryCatch(align(raw,comp), error=function(e) NULL)
  if (is.null(a) || nrow(a$m) < 12) { cat(sprintf("\n-- %s %s: too few trees\n",amp,comp)); return(invisible()) }
  cnt <- a$cnt; if (drop_org) cnt <- cnt[, !org[colnames(cnt)], drop=FALSE]
  sel <- rowSums(cnt) >= 250
  cnt <- cnt[sel,,drop=FALSE]; m <- a$m[sel,,drop=FALSE]
  cnt <- cnt[, colSums(cnt>0)>=5 & colSums(cnt)/sum(cnt)>=1e-4, drop=FALSE]
  if (ncol(cnt)<10 || nrow(cnt)<12) { cat(sprintf("\n-- %s %s: n=%d ASVs=%d, not fitted\n",amp,comp,nrow(cnt),ncol(cnt))); return(invisible()) }
  cv <- clim[match(m$tree_id, clim$tree_id),]
  ok <- complete.cases(cv[,c("bio7","bio13","dem_elev_m","lon","lat")])
  cnt<-cnt[ok,,drop=FALSE]; m<-m[ok,,drop=FALSE]; cv<-cv[ok,,drop=FALSE]
  if (nrow(cnt)<12) { cat(sprintf("\n-- %s %s: n=%d after env join\n",amp,comp,nrow(cnt))); return(invisible()) }
  D <- dist(clr(cnt))                                       # Aitchison
  CLIM <- scale(cv[,c("bio7","bio13","dem_elev_m")])
  pcnm_ <- tryCatch(pcnm(dist(cbind(cv$lon,cv$lat)))$vectors, error=function(e) NULL)
  GEO <- if (is.null(pcnm_)) NULL else scale(pcnm_[,seq_len(min(3,ncol(pcnm_))),drop=FALSE])
  Ks <- K[m$tree_id,m$tree_id]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)],na.rm=TRUE)
  diag(Ks) <- max(Ks,na.rm=TRUE)
  KIN <- cmdscale(as.dist(max(Ks)-Ks), k=min(3,nrow(Ks)-2))
  cat(sprintf("\n=========== %s %s: n=%d trees, %d ASVs ===========\n",amp,comp,nrow(cnt),ncol(cnt)))
  vp <- varpart(D, CLIM, GEO, KIN)
  print(vp$part$indfract[,c("Adj.R.square")], digits=4)
  cat("\n  fractions [a]=climate [b]=geography [c]=kinship, adj R2:\n")
  fr <- vp$part$indfract$Adj.R.square
  nm <- rownames(vp$part$indfract)
  for (i in seq_along(fr)) cat(sprintf("   %-14s %+8.4f\n", nm[i], fr[i]))
  # significance of each unique fraction
  for (lab in c("climate","geography","kinship")) {
    X <- switch(lab, climate=CLIM, geography=GEO, kinship=KIN)
    Z <- switch(lab, climate=cbind(GEO,KIN), geography=cbind(CLIM,KIN), kinship=cbind(CLIM,GEO))
    if (is.null(X)) next
    r <- tryCatch({ cap <- dbrda(D ~ X + Condition(Z)); an <- anova(cap, permutations=9999)
                    sprintf("F = %.3f, p = %.4f, partial R2 = %.4f", an$F[1], an$`Pr(>F)`[1], cap$CCA$tot.chi/cap$tot.chi) },
                  error=function(e) paste("failed:",conditionMessage(e)))
    cat(sprintf("   %-10s | others : %s\n", lab, r))
  }
}
cat("############ 16S above-ground (organellar removed) ############\n")
for (cp in c("needle","branch")) run(st16, cp, "16S", drop_org=TRUE)
cat("\n############ ITS above-ground ############\n")
for (cp in c("needle","branch")) run(stI, cp, "ITS")
cat("\nDONE\n")
