#!/usr/bin/env Rscript
# Cordillera split using needle AND branch as PAIRED observations of the same tree.
# Tree identity is the permutation stratum, so the non-independence of a tree's two
# compartments is respected rather than averaged away.
#   - cordillera is a BETWEEN-tree factor  -> permute whole trees (freely across trees)
#   - compartment is a WITHIN-tree factor  -> permute within tree
suppressPackageStartupMessages({library(vegan); library(permute)})
set.seed(2026); P <- "/home/ecastron/araucaria"
sink(file.path(P,"results/Q6/varpart/q6_paired_cordillera.txt"), split=TRUE)

st16 <- readRDS(file.path(P,"microbiome/16S/seqtab_nochim.rds"))
tx16 <- as.data.frame(readRDS(file.path(P,"microbiome/16S/taxonomy.rds")),stringsAsFactors=FALSE)
org <- rep(FALSE, ncol(st16))
for (rk in colnames(tx16)) org <- org | grepl("Chloroplast|Mitochondri", tx16[[rk]], ignore.case=TRUE)
org[is.na(org)] <- FALSE; names(org) <- colnames(st16)
stI <- readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))

rel <- read.delim(file.path(P,"results/Q2/ngsrelate.tsv"),check.names=FALSE)
radman <- read.csv(file.path(P,"metadata/manifest rad runs.csv"),stringsAsFactors=FALSE)
bare<-function(v) sub("\\.[12]$","",v)
tree_of<-setNames(radman$tree_id,radman$run_accession)
cord_of<-setNames(radman$cordillera,radman$tree_id)
ta<-tree_of[bare(as.character(rel$ida))]; tb<-tree_of[bare(as.character(rel$idb))]
gd<-!is.na(ta)&!is.na(tb)&ta!=tb&rel$nSites>=2000&is.finite(rel$theta)
trees<-sort(unique(c(ta[gd],tb[gd])))
K<-matrix(NA_real_,length(trees),length(trees),dimnames=list(trees,trees))
for (i in which(gd)){K[ta[i],tb[i]]<-rel$theta[i]; K[tb[i],ta[i]]<-rel$theta[i]}
diag(K)<-max(rel$theta[gd],na.rm=TRUE)
ch<-read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"),check.names=FALSE)
te<-read.delim(file.path(P,"metadata/terrain_per_tree.tsv"),check.names=FALSE)
env<-merge(ch[,c("tree_id","site_code","lon","lat","bio7","bio18","bio3")],
           te[,c("tree_id","dem_elev_m")],by="tree_id",all.x=TRUE)
# SoilGrids is per stand, so it joins on site_code. Soil was never in this model,
# yet clay is the second-strongest Coastal/Andean discriminator of the 39 variables
# tested (d = +1.87, exact permutation p = 0.0045) -- stronger than three of the
# four terms that WERE included, one of which (elevation) does not discriminate at
# all (p = 0.132).
so<-read.delim(file.path(P,"metadata/soil_per_stand.tsv"),check.names=FALSE)
env<-merge(env, so[,c("site_code","clay_0-5cm","cec_0-5cm")], by="site_code", all.x=TRUE)
names(env)[names(env)=="clay_0-5cm"]<-"clay"; names(env)[names(env)=="cec_0-5cm"]<-"cec"
clr<-function(m){m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0]); t(apply(m,1,function(x) log(x)-mean(log(x))))}

build <- function(raw, amp, drop_org=FALSE){
  parts<-list()
  for (comp in c("needle","branch")){
    o<-readRDS(file.path(P,sprintf("results/Q5/q5_%s_%s_objects.rds",amp,comp)))
    m<-o$manifest; m<-m[!is.na(m$tree_id),]
    ids<-intersect(rownames(raw),m$run_accession); if(!length(ids)) next
    m<-m[match(ids,m$run_accession),]
    x<-raw[ids,,drop=FALSE]; rownames(x)<-m$tree_id
    x<-x[rownames(x)%in%trees,,drop=FALSE]
    parts[[comp]]<-x
  }
  sh<-Reduce(intersect,lapply(parts,rownames))          # trees measured in BOTH
  cols<-Reduce(union,lapply(parts,colnames))
  mk<-function(M){A<-matrix(0,length(sh),length(cols),dimnames=list(sh,cols)); A[sh,colnames(M)]<-M[sh,]; A}
  X<-rbind(mk(parts$needle),mk(parts$branch))
  if (drop_org) X<-X[,!org[colnames(X)],drop=FALSE]
  meta<-data.frame(tree=rep(sh,2), comp=factor(rep(c("needle","branch"),each=length(sh))),
                   stringsAsFactors=FALSE)
  keep<-rowSums(X)>=250
  X<-X[keep,,drop=FALSE]; meta<-meta[keep,]
  # keep only trees still represented by BOTH compartments after the depth floor
  ok<-names(which(table(meta$tree)==2))
  sel<-meta$tree %in% ok; X<-X[sel,,drop=FALSE]; meta<-meta[sel,]
  X<-X[,colSums(X>0)>=5 & colSums(X)/sum(X)>=1e-4,drop=FALSE]
  meta$cord<-factor(ifelse(cord_of[meta$tree]=="Costa","Coastal","Andean"))
  ev<-env[match(meta$tree,env$tree_id),]
  list(X=X, meta=meta, ev=ev)
}

run <- function(raw, amp, drop_org=FALSE){
  b<-build(raw,amp,drop_org); X<-b$X; meta<-b$meta; ev<-b$ev
  ntree<-length(unique(meta$tree))
  cat(sprintf("\n=========== %s : needle + branch as PAIRED observations ===========\n",amp))
  cat(sprintf("  %d trees x 2 compartments = %d observations | %d ASVs\n", ntree, nrow(X), ncol(X)))
  cat(sprintf("  Coastal %d trees, Andean %d trees\n",
      length(unique(meta$tree[meta$cord=="Coastal"])), length(unique(meta$tree[meta$cord=="Andean"]))))
  if (ntree < 12){cat("  too few trees\n"); return(invisible())}
  D<-dist(clr(X))
  # BETWEEN-tree factor: permute whole trees
  hB <- how(blocks = factor(meta$tree), nperm = 9999)
  hB <- how(plots = Plots(strata=factor(meta$tree), type="free"),
            within = Within(type="none"), nperm = 9999)
  cat("\n  --- CORDILLERA (between-tree; whole trees permuted) ---\n")
  aC <- adonis2(D ~ cord, data=meta, permutations=hB)
  cat(sprintf("   R2 = %.4f  F = %.3f  p = %.4f\n", aC$R2[1], aC$F[1], aC$`Pr(>F)`[1]))
  bd<-betadisper(D, meta$cord); pb<-permutest(bd,permutations=9999)$tab
  cat(sprintf("   betadisper F = %.3f p = %.4f -> %s\n", pb[1,4], pb[1,6],
      if(pb[1,6]<0.05) "spread differs" else "spread homogeneous; a location effect"))
  cat("\n  --- COMPARTMENT (within-tree; permuted inside each tree) ---\n")
  hW <- how(plots = Plots(strata=factor(meta$tree), type="none"),
            within = Within(type="free"), nperm = 9999)
  aP <- adonis2(D ~ comp, data=meta, permutations=hW)
  cat(sprintf("   R2 = %.4f  F = %.3f  p = %.4f\n", aP$R2[1], aP$F[1], aP$`Pr(>F)`[1]))
  cat("\n  --- cordillera + compartment together ---\n")
  aB <- adonis2(D ~ cord + comp, data=meta, permutations=hB, by="terms")
  print(aB)
  cat("\n  --- cordillera conditioned on climate (dbRDA, trees permuted) ---\n")
  CLIM<-scale(ev[,c("bio7","bio18","bio3","dem_elev_m")]); colnames(CLIM)<-c("bio7","bio18","bio3","elev")
  SOIL<-scale(ev[,c("clay","cec")]); colnames(SOIL)<-c("clay","cec")
  ENV2<-cbind(CLIM,SOIL)
  # With 3 Coastal and 9 Andean stands, a variable that separates them completely is
  # near-collinear with the cordillera factor, so absorbing the cordillera term would
  # be partly STRUCTURAL rather than mechanistic. Report it so it cannot hide.
  cat("\n  --- how far each env term duplicates the cordillera factor ---\n")
  for (v in colnames(ENV2)) {
    r2<-summary(lm(ENV2[,v] ~ meta$cord))$r.squared
    cat(sprintf("   %-6s R2 vs cordillera = %.3f\n", v, r2))
  }
  cat("\n  --- cordillera conditioned on climate + SOIL ---\n")
  # Report the constrained rank explicitly. If the conditioning set spans the
  # cordillera contrast, dbrda leaves 0 constrained axes, anova returns an empty
  # table, sprintf() yields character(0) and cat() prints NOTHING -- a silent
  # non-result that looks identical to the line never running.
  cpd<-dbrda(D~meta$cord+Condition(ENV2))
  cat(sprintf("   constrained rank with climate+soil conditioned: %d  (0 = the cordillera term is not identifiable)\n",
      ifelse(is.null(cpd$CCA), 0L, cpd$CCA$rank)))
  cpc<-dbrda(D~meta$cord+Condition(CLIM))
  cat(sprintf("   constrained rank with climate only conditioned : %d\n",
      ifelse(is.null(cpc$CCA), 0L, cpc$CCA$rank)))
  cat(sprintf("   variance inflation, clay in the full env set    : %.1f\n",
      tryCatch(max(vif.cca(dbrda(D~ENV2))), error=function(e) NA_real_)))
  rs<-tryCatch({cp<-dbrda(D~meta$cord+Condition(ENV2)); an<-anova(cp,permutations=hB)
        sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$"Pr(>F)"[1])},
        error=function(e) paste("failed:",conditionMessage(e)))
  cat(sprintf("   cordillera | climate + soil : %s\n", rs))
  rs2<-tryCatch({CD<-model.matrix(~meta$cord)[,-1,drop=FALSE]
        cp<-dbrda(D~SOIL+Condition(cbind(CLIM,CD))); an<-anova(cp,permutations=hB)
        sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$"Pr(>F)"[1])},
        error=function(e) paste("failed:",conditionMessage(e)))
  cat(sprintf("   soil | climate + cordillera : %s\n", rs2))
  rs3<-tryCatch({cp<-dbrda(D~SOIL+Condition(CLIM)); an<-anova(cp,permutations=hB)
        sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$"Pr(>F)"[1])},
        error=function(e) paste("failed:",conditionMessage(e)))
  cat(sprintf("   soil | climate (no cordillera) : %s\n", rs3))
  r1<-tryCatch({cp<-dbrda(D~meta$cord+Condition(CLIM)); an<-anova(cp,permutations=hB)
      sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$`Pr(>F)`[1])},
      error=function(e) paste("failed:",conditionMessage(e)))
  cat(sprintf("   cordillera | climate : %s\n", r1))
  r2_<-tryCatch({cp<-dbrda(D~CLIM+Condition(meta$cord)); an<-anova(cp,permutations=hB)
      sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$`Pr(>F)`[1])},
      error=function(e) paste("failed:",conditionMessage(e)))
  cat(sprintf("   climate | cordillera : %s\n", r2_))
  cat("\n  --- host genetics, same paired design ---\n")
  Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
  KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
  rk<-tryCatch({cp<-dbrda(D~KIN+Condition(CLIM+meta$comp)); an<-anova(cp,permutations=hB)
      sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$`Pr(>F)`[1])},
      error=function(e) paste("failed:",conditionMessage(e)))
  cat(sprintf("   kinship | climate + compartment : %s\n", rk))
}
run(stI,"ITS")
run(st16,"16S",drop_org=TRUE)
cat("\nDONE\n")
