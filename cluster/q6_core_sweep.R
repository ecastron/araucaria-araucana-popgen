#!/usr/bin/env Rscript
# CORE MICROBIOME sweep on the paired needle+branch ITS design.
# core(detection, prevalence) = taxa reaching >= detection relative abundance
# in >= prevalence fraction of the 86 samples. Tree is the permutation unit throughout.
suppressPackageStartupMessages({library(vegan); library(permute)})
set.seed(2026); P<-"/home/ecastron/araucaria"
sink(file.path(P,"results/Q6/varpart/q6_core_sweep.txt"), split=TRUE)
stI<-readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))
rel<-read.delim(file.path(P,"results/Q2/ngsrelate.tsv"),check.names=FALSE)
radman<-read.csv(file.path(P,"metadata/manifest rad runs.csv"),stringsAsFactors=FALSE)
bare<-function(v) sub("\\.[12]$","",v)
tree_of<-setNames(radman$tree_id,radman$run_accession); cord_of<-setNames(radman$cordillera,radman$tree_id)
ta<-tree_of[bare(as.character(rel$ida))]; tb<-tree_of[bare(as.character(rel$idb))]
gd<-!is.na(ta)&!is.na(tb)&ta!=tb&rel$nSites>=2000&is.finite(rel$theta)
trees<-sort(unique(c(ta[gd],tb[gd])))
K<-matrix(NA_real_,length(trees),length(trees),dimnames=list(trees,trees))
for(i in which(gd)){K[ta[i],tb[i]]<-rel$theta[i];K[tb[i],ta[i]]<-rel$theta[i]}
diag(K)<-max(rel$theta[gd],na.rm=TRUE)
ch<-read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"),check.names=FALSE)
te<-read.delim(file.path(P,"metadata/terrain_per_tree.tsv"),check.names=FALSE)
env<-merge(ch[,c("tree_id","site_code","bio7","bio18","bio3")],te[,c("tree_id","dem_elev_m")],by="tree_id",all.x=TRUE)
clr<-function(m){m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0]); t(apply(m,1,function(x) log(x)-mean(log(x))))}

parts<-list()
for(comp in c("needle","branch")){
  o<-readRDS(file.path(P,sprintf("results/Q5/q5_ITS_%s_objects.rds",comp)))
  m<-o$manifest; m<-m[!is.na(m$tree_id),]
  ids<-intersect(rownames(stI),m$run_accession); m<-m[match(ids,m$run_accession),]
  x<-stI[ids,,drop=FALSE]; rownames(x)<-m$tree_id
  parts[[comp]]<-x[rownames(x)%in%trees,,drop=FALSE]
}
sh<-intersect(rownames(parts$needle),rownames(parts$branch))
cols<-union(colnames(parts$needle),colnames(parts$branch))
mk<-function(M){A<-matrix(0,length(sh),length(cols),dimnames=list(sh,cols)); A[sh,colnames(M)]<-M[sh,]; A}
RAW<-rbind(mk(parts$needle),mk(parts$branch))
meta<-data.frame(tree=rep(sh,2),comp=factor(rep(c("needle","branch"),each=length(sh))),stringsAsFactors=FALSE)
keep<-rowSums(RAW)>=250; RAW<-RAW[keep,,drop=FALSE]; meta<-meta[keep,]
ok<-names(which(table(meta$tree)==2)); sel<-meta$tree%in%ok
RAW<-RAW[sel,,drop=FALSE]; meta<-meta[sel,]
meta$cord<-factor(ifelse(cord_of[meta$tree]=="Costa","Coastal","Andean"))
ev<-env[match(meta$tree,env$tree_id),]
CLIM<-scale(ev[,c("bio7","bio18","bio3","dem_elev_m")]); colnames(CLIM)<-c("bio7","bio18","bio3","elev")
Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
REL<-RAW/rowSums(RAW)
hB<-how(plots=Plots(strata=factor(meta$tree),type="free"),within=Within(type="none"),nperm=9999)
cat(sprintf("paired ITS design: %d trees x 2 = %d samples, %d ASVs before coring\n",
            length(unique(meta$tree)), nrow(RAW), ncol(RAW)))
cat(sprintf("Coastal %d trees, Andean %d trees\n\n",
  length(unique(meta$tree[meta$cord=="Coastal"])), length(unique(meta$tree[meta$cord=="Andean"]))))

cat(sprintf("%-9s %-7s %6s %8s %10s %9s %10s %9s %10s %9s\n",
  "detection","prev","ASVs","%reads","cord R2","cord p","cord|clim p","clim|cord p","kin R2","kin p"))
for (det in c(1e-4, 1e-3)) {
  for (prv in c(0.25,0.40,0.50,0.60,0.70,0.80,0.90)) {
    inc <- colSums(REL >= det) >= prv*nrow(REL)
    if (sum(inc) < 8) { cat(sprintf("%-9.0e %-7.2f %6d   (<8 ASVs, skipped)\n",det,prv,sum(inc))); next }
    X<-RAW[,inc,drop=FALSE]
    pr<-100*sum(X)/sum(RAW)
    D<-dist(clr(X))
    aC<-adonis2(D~cord,data=meta,permutations=hB)
    r1<-tryCatch({cp<-dbrda(D~meta$cord+Condition(CLIM)); anova(cp,permutations=hB)$`Pr(>F)`[1]},error=function(e) NA)
    r2<-tryCatch({cp<-dbrda(D~CLIM+Condition(meta$cord)); anova(cp,permutations=hB)$`Pr(>F)`[1]},error=function(e) NA)
    rk<-tryCatch({cp<-dbrda(D~KIN+Condition(CLIM+meta$comp)); a<-anova(cp,permutations=hB)
                  c(cp$CCA$tot.chi/cp$tot.chi, a$`Pr(>F)`[1])},error=function(e) c(NA,NA))
    cat(sprintf("%-9.0e %-7.2f %6d %8.1f %10.4f %9.4f %10.4f %11.4f %10.4f %9.4f\n",
        det,prv,sum(inc),pr,aC$R2[1],aC$`Pr(>F)`[1],r1,r2,rk[1],rk[2]))
  }
  cat("\n")
}
cat("DONE\n")
