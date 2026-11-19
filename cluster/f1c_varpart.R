#!/usr/bin/env Rscript
# What IS the 0.0448 that kinship shows when soil is left out of the model?
# Variance partitioning of kinship against soil, both conditioned on climate + compartment.
# If the kinship signal is soil-shaped, most of it sits in the SHARED fraction and its unique
# fraction is near zero.
suppressPackageStartupMessages({library(vegan); library(permute)})
set.seed(2026); P<-"/home/ecastron/araucaria"
sink(file.path(P,"results/F1C/f1c_varpart.txt"),split=TRUE)
stI<-readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))
tax<-readRDS(file.path(P,"microbiome/ITS/taxonomy.rds"))
lab<-as.character(tax[,"Family"])
for (r in c("Order","Class","Phylum","Kingdom")){
  need<-is.na(lab); if(!any(need)) break
  src<-as.character(tax[,r]); ok<-need&!is.na(src)
  lab[ok]<-paste0(tolower(substr(r,1,1)),"__",src[ok],"_unclassified") }
lab[is.na(lab)]<-"k__Fungi_unclassified"
FAMP<-t(rowsum(t(stI),lab))
rel<-read.delim(file.path(P,"results/Q2/ngsrelate.tsv"),check.names=FALSE)
radman<-read.csv(file.path(P,"metadata/manifest rad runs.csv"),stringsAsFactors=FALSE)
bare<-function(v) sub("\\.[12]$","",v)
tree_of<-setNames(radman$tree_id,radman$run_accession)
ta<-tree_of[bare(as.character(rel$ida))]; tb<-tree_of[bare(as.character(rel$idb))]
gd<-!is.na(ta)&!is.na(tb)&ta!=tb&rel$nSites>=2000&is.finite(rel$theta)
trees<-sort(unique(c(ta[gd],tb[gd])))
K<-matrix(NA_real_,length(trees),length(trees),dimnames=list(trees,trees))
for (i in which(gd)){K[ta[i],tb[i]]<-rel$theta[i]; K[tb[i],ta[i]]<-rel$theta[i]}
diag(K)<-max(rel$theta[gd],na.rm=TRUE)
ch<-read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"),check.names=FALSE)
env<-ch[,c("tree_id","site_code","lon","lat","bio5","bio18","bio6")]
so<-read.delim(file.path(P,"metadata/soil_per_stand.tsv"),check.names=FALSE)
env<-merge(env,so[,c("site_code","clay_0-5cm","cec_0-5cm","phh2o_0-5cm","nitrogen_0-5cm")],by="site_code",all.x=TRUE)
names(env)[names(env)=="clay_0-5cm"]<-"clay"; names(env)[names(env)=="cec_0-5cm"]<-"cec"
names(env)[names(env)=="phh2o_0-5cm"]<-"ph"; names(env)[names(env)=="nitrogen_0-5cm"]<-"nitro"
clr<-function(m){m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0]); t(apply(m,1,function(x) log(x)-mean(log(x))))}
parts<-list()
for (comp in c("needle","branch")){
  o<-readRDS(file.path(P,sprintf("results/Q5/q5_ITS_%s_objects.rds",comp)))
  m<-o$manifest; m<-m[!is.na(m$tree_id),]
  ids<-intersect(rownames(FAMP),m$run_accession); m<-m[match(ids,m$run_accession),]
  x<-FAMP[ids,,drop=FALSE]; rownames(x)<-m$tree_id
  parts[[comp]]<-x[rownames(x)%in%trees,,drop=FALSE] }
sh<-Reduce(intersect,lapply(parts,rownames)); cols<-Reduce(union,lapply(parts,colnames))
mk<-function(M){A<-matrix(0,length(sh),length(cols),dimnames=list(sh,cols)); A[sh,colnames(M)]<-M[sh,]; A}
X<-rbind(mk(parts$needle),mk(parts$branch))
meta<-data.frame(tree=rep(sh,2),comp=factor(rep(c("needle","branch"),each=length(sh))))
keep<-rowSums(X)>=250; X<-X[keep,,drop=FALSE]; meta<-meta[keep,]
ok<-names(which(table(meta$tree)==2)); sel<-meta$tree%in%ok
X<-X[sel,,drop=FALSE]; meta<-meta[sel,]
X<-X[,colSums(X>0)>=5 & colSums(X)/sum(X)>=1e-4,drop=FALSE]
ev<-env[match(meta$tree,env$tree_id),]
D<-dist(clr(X))
CL<-scale(ev[,c("bio5","bio18","bio6")]); SO<-scale(ev[,c("cec","ph","nitro","clay")])
Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
CO<-model.matrix(~meta$comp)[,-1,drop=FALSE]
R2<-function(f) {cp<-tryCatch(eval(f),error=function(e) NULL)
  if(is.null(cp)||is.null(cp$CCA)) NA_real_ else cp$CCA$tot.chi/cp$tot.chi}
a<-R2(quote(dbrda(D~KIN+Condition(cbind(CL,CO)))))          # kinship | climate + comp
b<-R2(quote(dbrda(D~SO +Condition(cbind(CL,CO)))))          # soil    | climate + comp
ab<-R2(quote(dbrda(D~cbind(KIN,SO)+Condition(cbind(CL,CO))))) # both together
uk<-ab-b; us<-ab-a; shared<-a+b-ab
cat("Variance in fungal composition, both conditioned on climate + compartment\n")
cat(sprintf("  kinship, alone in the model      %.4f\n",a))
cat(sprintf("  soil,    alone in the model      %.4f\n",b))
cat(sprintf("  both together                    %.4f\n\n",ab))
cat(sprintf("  kinship UNIQUE (soil also fitted) %.4f\n",uk))
cat(sprintf("  soil    UNIQUE (kin also fitted)  %.4f\n",us))
cat(sprintf("  SHARED by kinship and soil        %.4f  (%.0f%% of the kinship term)\n",
    shared,100*shared/a))
h<-how(plots=Plots(strata=factor(meta$tree),type="free"),within=Within(type="none"),nperm=9999)
pk<-tryCatch(anova(dbrda(D~KIN+Condition(cbind(SO,CL,CO))),permutations=h)$`Pr(>F)`[1],error=function(e) NA)
ps<-tryCatch(anova(dbrda(D~SO +Condition(cbind(KIN,CL,CO))),permutations=h)$`Pr(>F)`[1],error=function(e) NA)
cat(sprintf("\n  p for kinship UNIQUE fraction     %.4f\n",pk))
cat(sprintf("  p for soil    UNIQUE fraction     %.4f\n",ps))
# how collinear are they at the level that matters
cat("\nCollinearity check: soil variables predicted by the 3 kinship axes (tree level)\n")
for (v in colnames(SO)) cat(sprintf("  %-6s R2 = %.3f\n",v,summary(lm(SO[,v]~KIN))$r.squared))
cat("\nDONE\n"); sink()
