#!/usr/bin/env Rscript
# ITS at FAMILY level vs ASV level, same paired needle+branch design throughout.
# Rationale: host genotype may select on broad functional/taxonomic groups rather than on exact
# sequence variants, in which case ASV resolution scatters a real signal across many rare columns.
# 16S is excluded entirely (91.6% chloroplast+mitochondria).
#
# Agglomeration sums counts across ASVs sharing a Family; ASVs with unassigned Family are dropped
# (28.9% of ASVs but only 15.9% of reads). Every downstream step is identical to the ASV pipeline
# so any difference is attributable to resolution, not to a change in design.
suppressPackageStartupMessages({library(vegan); library(permute)})
set.seed(2026); P <- "/home/ecastron/araucaria"
sink(file.path(P,"results/Q6/family/q6_family_level.txt"), split=TRUE)

stI <- readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))
tax <- readRDS(file.path(P,"microbiome/ITS/taxonomy.rds"))
stopifnot(identical(rownames(tax), colnames(stI)))

fam <- tax[,"Family"]
cat("=================== FAMILY AGGLOMERATION ===================\n")
cat(sprintf("ITS ASV table : %d samples x %d ASVs, %s reads\n", nrow(stI), ncol(stI), format(sum(stI),big.mark=",")))
keepf <- !is.na(fam)
cat(sprintf("Family assigned: %d/%d ASVs (%.1f%%), %.1f%% of reads\n",
    sum(keepf), length(fam), 100*mean(keepf), 100*sum(stI[,keepf])/sum(stI)))
cat(sprintf("DROPPED (unassigned Family): %d ASVs, %.1f%% of reads\n",
    sum(!keepf), 100*sum(stI[,!keepf])/sum(stI)))
FAM <- t(rowsum(t(stI[,keepf,drop=FALSE]), fam[keepf]))
cat(sprintf("Family table  : %d samples x %d families\n\n", nrow(FAM), ncol(FAM)))

rel <- read.delim(file.path(P,"results/Q2/ngsrelate.tsv"),check.names=FALSE)
radman <- read.csv(file.path(P,"metadata/manifest rad runs.csv"),stringsAsFactors=FALSE)
bare<-function(v) sub("\\.[12]$","",v)
tree_of<-setNames(radman$tree_id,radman$run_accession); cord_of<-setNames(radman$cordillera,radman$tree_id)
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
clr<-function(m){m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0]); t(apply(m,1,function(x) log(x)-mean(log(x))))}

build <- function(raw){
  parts<-list()
  for (comp in c("needle","branch")){
    o<-readRDS(file.path(P,sprintf("results/Q5/q5_ITS_%s_objects.rds",comp)))
    m<-o$manifest; m<-m[!is.na(m$tree_id),]
    ids<-intersect(rownames(raw),m$run_accession); if(!length(ids)) next
    m<-m[match(ids,m$run_accession),]
    x<-raw[ids,,drop=FALSE]; rownames(x)<-m$tree_id
    parts[[comp]]<-x[rownames(x)%in%trees,,drop=FALSE]
  }
  sh<-Reduce(intersect,lapply(parts,rownames)); cols<-Reduce(union,lapply(parts,colnames))
  mk<-function(M){A<-matrix(0,length(sh),length(cols),dimnames=list(sh,cols)); A[sh,colnames(M)]<-M[sh,]; A}
  X<-rbind(mk(parts$needle),mk(parts$branch))
  meta<-data.frame(tree=rep(sh,2),comp=factor(rep(c("needle","branch"),each=length(sh))),stringsAsFactors=FALSE)
  keep<-rowSums(X)>=250; X<-X[keep,,drop=FALSE]; meta<-meta[keep,]
  ok<-names(which(table(meta$tree)==2)); sel<-meta$tree%in%ok
  X<-X[sel,,drop=FALSE]; meta<-meta[sel,]
  RAWfull <- X
  X<-X[,colSums(X>0)>=5 & colSums(X)/sum(X)>=1e-4,drop=FALSE]
  meta$cord<-factor(ifelse(cord_of[meta$tree]=="Costa","Coastal","Andean"))
  ev<-env[match(meta$tree,env$tree_id),]
  list(X=X, RAW=RAWfull, meta=meta, ev=ev)
}

run <- function(raw, tag){
  b<-build(raw); X<-b$X; meta<-b$meta; ev<-b$ev
  ntree<-length(unique(meta$tree))
  cat(sprintf("\n############### %s ###############\n", tag))
  cat(sprintf("  %d trees x 2 compartments = %d observations | %d features (of %d before the filter)\n",
      ntree,nrow(X),ncol(X),ncol(b$RAW)))
  cat(sprintf("  Coastal %d trees, Andean %d trees | sparsity %.1f%% zeros\n",
      length(unique(meta$tree[meta$cord=="Coastal"])),
      length(unique(meta$tree[meta$cord=="Andean"])), 100*mean(X==0)))
  D<-dist(clr(X))
  hB <- how(plots=Plots(strata=factor(meta$tree),type="free"), within=Within(type="none"), nperm=9999)
  hW <- how(plots=Plots(strata=factor(meta$tree),type="none"), within=Within(type="free"), nperm=9999)
  CLIM<-scale(ev[,c("bio7","bio18","bio3","dem_elev_m")]); colnames(CLIM)<-c("bio7","bio18","bio3","elev")
  Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
  KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
  aC<-adonis2(D~cord,data=meta,permutations=hB)
  cat(sprintf("\n  CORDILLERA (trees permuted)     R2=%.4f F=%.3f p=%.4f\n",aC$R2[1],aC$F[1],aC$`Pr(>F)`[1]))
  bd<-betadisper(D,meta$cord); pb<-permutest(bd,permutations=9999)$tab
  cat(sprintf("    betadisper F=%.3f p=%.4f -> %s\n",pb[1,4],pb[1,6],
      if(pb[1,6]<0.05) "spread differs" else "spread homogeneous; location effect"))
  aS<-adonis2(D~site_code,data=data.frame(site_code=factor(ev$site_code)),permutations=hB)
  cat(sprintf("  STAND (Q5, trees permuted)      R2=%.4f F=%.3f p=%.4f\n",aS$R2[1],aS$F[1],aS$`Pr(>F)`[1]))
  aP<-adonis2(D~comp,data=meta,permutations=hW)
  cat(sprintf("  COMPARTMENT (within tree)       R2=%.4f F=%.3f p=%.4f\n",aP$R2[1],aP$F[1],aP$`Pr(>F)`[1]))
  f<-function(expr) tryCatch({cp<-eval(expr); an<-anova(cp,permutations=hB)
        sprintf("R2=%.4f F=%.3f p=%.4f",cp$CCA$tot.chi/cp$tot.chi,an$F[1],an$`Pr(>F)`[1])},
        error=function(e) paste("failed:",conditionMessage(e)))
  cat(sprintf("  cordillera | climate            %s\n", f(quote(dbrda(D~meta$cord+Condition(CLIM))))))
  cat(sprintf("  climate | cordillera            %s\n", f(quote(dbrda(D~CLIM+Condition(meta$cord))))))
  cat(sprintf("  HOST GENETICS kin|clim+comp     %s\n", f(quote(dbrda(D~KIN+Condition(CLIM+meta$comp))))))
  cat(sprintf("  HOST GENETICS kin|clim+comp+cord %s\n",f(quote(dbrda(D~KIN+Condition(CLIM+meta$comp+meta$cord))))))
  invisible(list(meta=meta,CLIM=CLIM,KIN=KIN,hB=hB))
}
rA <- run(stI,"ITS at ASV level (reference)")
rF <- run(FAM,"ITS at FAMILY level")

# ---- core sweep at family level ----
cat("\n\n=================== CORE MICROBIOME, FAMILY LEVEL ===================\n")
b<-build(FAM); RAW<-b$RAW; meta<-b$meta; ev<-b$ev
REL<-RAW/rowSums(RAW)
CLIM<-scale(ev[,c("bio7","bio18","bio3","dem_elev_m")]); colnames(CLIM)<-c("bio7","bio18","bio3","elev")
Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
hB<-how(plots=Plots(strata=factor(meta$tree),type="free"),within=Within(type="none"),nperm=9999)
cat(sprintf("paired design: %d trees x 2 = %d samples, %d families before coring\n\n",
    length(unique(meta$tree)),nrow(RAW),ncol(RAW)))
cat(sprintf("%-9s %-6s %6s %8s %9s %8s %11s %11s %9s %8s\n",
    "detection","prev","fams","%reads","cord R2","cord p","cord|clim p","clim|cord p","kin R2","kin p"))
for (det in c(1e-4,1e-3)){
  for (prv in c(0.25,0.40,0.50,0.60,0.70,0.80,0.90)){
    inc<-colSums(REL>=det)>=prv*nrow(REL)
    if (sum(inc)<8){cat(sprintf("%-9.0e %-6.2f %6d   (<8 families, skipped)\n",det,prv,sum(inc))); next}
    X<-RAW[,inc,drop=FALSE]; D<-dist(clr(X))
    aC<-adonis2(D~cord,data=meta,permutations=hB)
    r1<-tryCatch({cp<-dbrda(D~meta$cord+Condition(CLIM)); anova(cp,permutations=hB)$`Pr(>F)`[1]},error=function(e) NA)
    r2<-tryCatch({cp<-dbrda(D~CLIM+Condition(meta$cord)); anova(cp,permutations=hB)$`Pr(>F)`[1]},error=function(e) NA)
    rk<-tryCatch({cp<-dbrda(D~KIN+Condition(CLIM+meta$comp)); a<-anova(cp,permutations=hB)
                  c(cp$CCA$tot.chi/cp$tot.chi,a$`Pr(>F)`[1])},error=function(e) c(NA,NA))
    cat(sprintf("%-9.0e %-6.2f %6d %8.1f %9.4f %8.4f %11.4f %11.4f %9.4f %8.4f\n",
        det,prv,sum(inc),100*sum(X)/sum(RAW),aC$R2[1],aC$`Pr(>F)`[1],r1,r2,rk[1],rk[2]))
  }
  cat("\n")
}
cat("DONE\n")
