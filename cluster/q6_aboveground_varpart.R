#!/usr/bin/env Rscript
# ABOVE-GROUND endophyte microbiome: needle + branch POOLED per tree, for power.
# Climate set chosen from the A. araucana literature, all VIF < 6, max |r| = 0.375:
#   bio7  temperature annual range      (continentality; top variable in Varas-Myrik 2022, McIntosh 2022)
#   bio18 precipitation warmest quarter (growing-season water; Mundo 2011, Munoz 2014)
#   elev  elevation                     (biophysical gradient; also carries winter cold, r=-0.91 with bio6)
#   bio3  isothermality                 (fourth independent axis)
# Tests: cordillera (Coastal vs Andean) as a factor, and variance partitioning
#        CLIMATE | GEOGRAPHY | KINSHIP.
suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
dir.create(file.path(P,"results/Q6/varpart"), showWarnings=FALSE, recursive=TRUE)
sink(file.path(P,"results/Q6/varpart/q6_aboveground.txt"), split=TRUE)

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
cord_of <- setNames(radman$cordillera, radman$tree_id)
ta <- tree_of[bare(as.character(rel$ida))]; tb <- tree_of[bare(as.character(rel$idb))]
gd <- !is.na(ta)&!is.na(tb)&ta!=tb&rel$nSites>=2000&is.finite(rel$theta)
trees <- sort(unique(c(ta[gd],tb[gd])))
K <- matrix(NA_real_,length(trees),length(trees),dimnames=list(trees,trees))
for (i in which(gd)) { K[ta[i],tb[i]]<-rel$theta[i]; K[tb[i],ta[i]]<-rel$theta[i] }
diag(K) <- max(rel$theta[gd],na.rm=TRUE)

ch <- read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"), check.names=FALSE)
te <- read.delim(file.path(P,"metadata/terrain_per_tree.tsv"), check.names=FALSE)
env <- merge(ch[,c("tree_id","site_code","lon","lat","bio7","bio18","bio3")],
             te[,c("tree_id","dem_elev_m")], by="tree_id", all.x=TRUE)

clr <- function(m){ m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0])
  t(apply(m,1,function(x) log(x)-mean(log(x)))) }

# pool needle + branch per TREE
pool <- function(raw, amp) {
  acc <- NULL
  for (comp in c("needle","branch")) {
    o <- readRDS(file.path(P,"results/Q5",sprintf("q5_%s_%s_objects.rds",amp,comp)))
    m <- o$manifest; m <- m[!is.na(m$tree_id),]
    ids <- intersect(rownames(raw), m$run_accession)
    if (!length(ids)) next
    m <- m[match(ids,m$run_accession),]
    sub <- raw[ids,,drop=FALSE]; rownames(sub) <- m$tree_id
    acc <- if (is.null(acc)) sub else {
      allt <- union(rownames(acc),rownames(sub)); allc <- union(colnames(acc),colnames(sub))
      A <- matrix(0,length(allt),length(allc),dimnames=list(allt,allc))
      A[rownames(acc),colnames(acc)] <- acc
      A[rownames(sub),colnames(sub)] <- A[rownames(sub),colnames(sub)] + sub
      A }
  }
  acc
}

run <- function(raw, amp, drop_org=FALSE) {
  cnt <- pool(raw, amp)
  if (drop_org) cnt <- cnt[, !org[colnames(cnt)], drop=FALSE]
  keep <- rownames(cnt) %in% rownames(K)
  cnt <- cnt[keep,,drop=FALSE]
  cnt <- cnt[rowSums(cnt)>=250,,drop=FALSE]
  cnt <- cnt[, colSums(cnt>0)>=5 & colSums(cnt)/sum(cnt)>=1e-4, drop=FALSE]
  tr <- rownames(cnt)
  ev <- env[match(tr, env$tree_id),]
  ok <- complete.cases(ev[,c("bio7","bio18","bio3","dem_elev_m","lon","lat")])
  cnt<-cnt[ok,,drop=FALSE]; ev<-ev[ok,,drop=FALSE]; tr<-tr[ok]
  cd <- factor(ifelse(cord_of[tr]=="Costa","Coastal","Andean"))
  cat(sprintf("\n=========== %s ABOVE-GROUND (needle + branch pooled) ===========\n", amp))
  cat(sprintf("  %d trees, %d ASVs | Coastal %d, Andean %d | %d stands\n",
      nrow(cnt), ncol(cnt), sum(cd=="Coastal"), sum(cd=="Andean"), length(unique(ev$site_code))))
  if (nrow(cnt) < 15) { cat("  too few trees\n"); return(invisible()) }
  D <- dist(clr(cnt))
  CLIM <- scale(ev[,c("bio7","bio18","bio3","dem_elev_m")])
  colnames(CLIM) <- c("bio7","bio18","bio3","elev")
  pc <- tryCatch(pcnm(dist(cbind(ev$lon,ev$lat)))$vectors, error=function(e) NULL)
  GEO <- if(is.null(pc)) NULL else scale(pc[,seq_len(min(3,ncol(pc))),drop=FALSE])
  Ks <- K[tr,tr]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
  KIN <- cmdscale(as.dist(max(Ks)-Ks), k=min(3,nrow(Ks)-2))

  cat("\n  --- CORDILLERA (Coastal vs Andean) ---\n")
  a <- adonis2(D ~ cd, permutations=9999)
  cat(sprintf("   PERMANOVA : R2 = %.4f  F = %.3f  p = %.4f\n", a$R2[1], a$F[1], a$`Pr(>F)`[1]))
  bd <- betadisper(D, cd); pb <- permutest(bd, permutations=9999)$tab
  cat(sprintf("   betadisper: F = %.3f  p = %.4f  -> %s\n", pb[1,4], pb[1,6],
      if (pb[1,6] < 0.05) "spread differs; PERMANOVA partly dispersion" else "spread homogeneous; a location effect"))
  ac <- tryCatch({ cp <- dbrda(D ~ cd + Condition(CLIM)); an <- anova(cp, permutations=9999)
      sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$`Pr(>F)`[1]) },
      error=function(e) paste("failed:", conditionMessage(e)))
  cat(sprintf("   cordillera | climate : %s\n", ac))
  ag <- tryCatch({ cp <- dbrda(D ~ CLIM + Condition(cd)); an <- anova(cp, permutations=9999)
      sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$`Pr(>F)`[1]) },
      error=function(e) paste("failed:", conditionMessage(e)))
  cat(sprintf("   climate | cordillera : %s\n", ag))

  cat("\n  --- VARIANCE PARTITION: climate | geography | kinship ---\n")
  vp <- varpart(D, CLIM, GEO, KIN)
  fr <- vp$part$indfract$Adj.R.square; nm <- rownames(vp$part$indfract)
  for (i in seq_along(fr)) cat(sprintf("   %-16s %+8.4f\n", nm[i], fr[i]))
  for (lab in c("climate","geography","kinship")) {
    X <- switch(lab, climate=CLIM, geography=GEO, kinship=KIN)
    Z <- switch(lab, climate=cbind(GEO,KIN), geography=cbind(CLIM,KIN), kinship=cbind(CLIM,GEO))
    if (is.null(X)) next
    r <- tryCatch({ cp<-dbrda(D ~ X + Condition(Z)); an<-anova(cp,permutations=9999)
      sprintf("F = %.3f, p = %.4f, partial R2 = %.4f", an$F[1], an$`Pr(>F)`[1], cp$CCA$tot.chi/cp$tot.chi) },
      error=function(e) paste("failed:",conditionMessage(e)))
    cat(sprintf("   %-10s | others : %s\n", lab, r))
  }
  cat("\n  --- individual climate terms, each conditioned on the rest ---\n")
  for (v in colnames(CLIM)) {
    Z <- CLIM[,setdiff(colnames(CLIM),v),drop=FALSE]; X <- CLIM[,v,drop=FALSE]
    r <- tryCatch({ cp<-dbrda(D ~ X + Condition(Z)); an<-anova(cp,permutations=9999)
      sprintf("F = %.3f, p = %.4f", an$F[1], an$`Pr(>F)`[1]) }, error=function(e) "failed")
    cat(sprintf("   %-6s : %s\n", v, r))
  }
}
run(st16,"16S",drop_org=TRUE)
run(stI,"ITS")
cat("\nDONE\n")
