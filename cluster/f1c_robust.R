#!/usr/bin/env Rscript
# Is the host-genetic null robust to the conditioning set? The 4 Sep panel value (R2 0.0335,
# p 0.658) used the RETIRED climate set. Under the final 3-variable climate set the same term at
# the same resolution gives 0.0448, p 0.0564 -- above its own null median. That is a different
# statement, so every plausible specification is run here with its OWN shuffled-label null.
suppressPackageStartupMessages({library(vegan); library(permute)})
set.seed(2026); P <- "/home/ecastron/araucaria"
sink(file.path(P,"results/F1C/f1c_robust.txt"), split=TRUE)
stI <- readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))
tax <- readRDS(file.path(P,"microbiome/ITS/taxonomy.rds"))
fam <- tax[,"Family"]; lab <- as.character(fam)
for (r in c("Order","Class","Phylum","Kingdom")){
  need <- is.na(lab); if(!any(need)) break
  src <- as.character(tax[,r]); ok <- need & !is.na(src)
  lab[ok] <- paste0(tolower(substr(r,1,1)),"__",src[ok],"_unclassified") }
lab[is.na(lab)] <- "k__Fungi_unclassified"
FAMP <- t(rowsum(t(stI), lab))
rel<-read.delim(file.path(P,"results/Q2/ngsrelate.tsv"),check.names=FALSE)
radman<-read.csv(file.path(P,"metadata/manifest rad runs.csv"),stringsAsFactors=FALSE)
bare<-function(v) sub("\\.[12]$","",v)
tree_of<-setNames(radman$tree_id,radman$run_accession); cord_of<-setNames(radman$cordillera,radman$tree_id)
ta<-tree_of[bare(as.character(rel$ida))]; tb<-tree_of[bare(as.character(rel$idb))]
gd<-!is.na(ta)&!is.na(tb)&ta!=tb&rel$nSites>=2000&is.finite(rel$theta)
trees<-sort(unique(c(ta[gd],tb[gd])))
K<-matrix(NA_real_,length(trees),length(trees),dimnames=list(trees,trees))
for (i in which(gd)){K[ta[i],tb[i]]<-rel$theta[i]; K[tb[i],ta[i]]<-rel$theta[i]}
diag(K)<-max(rel$theta[gd],na.rm=TRUE)
ch<-read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"),check.names=FALSE)
env<-ch[,c("tree_id","site_code","bio5","bio18","bio6","bio7","bio3")]
te<-read.delim(file.path(P,"metadata/terrain_per_tree.tsv"),check.names=FALSE)
env<-merge(env,te[,c("tree_id","dem_elev_m")],by="tree_id",all.x=TRUE)
so<-read.delim(file.path(P,"metadata/soil_per_stand.tsv"),check.names=FALSE)
env<-merge(env,so[,c("site_code","clay_0-5cm","cec_0-5cm","phh2o_0-5cm","nitrogen_0-5cm")],by="site_code",all.x=TRUE)
names(env)[names(env)=="clay_0-5cm"]<-"clay"; names(env)[names(env)=="cec_0-5cm"]<-"cec"
names(env)[names(env)=="phh2o_0-5cm"]<-"ph"; names(env)[names(env)=="nitrogen_0-5cm"]<-"nitro"
clr<-function(m){m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0]); t(apply(m,1,function(x) log(x)-mean(log(x))))}
build <- function(raw){
  parts<-list()
  for (comp in c("needle","branch")){
    o<-readRDS(file.path(P,sprintf("results/Q5/q5_ITS_%s_objects.rds",comp)))
    m<-o$manifest; m<-m[!is.na(m$tree_id),]
    ids<-intersect(rownames(raw),m$run_accession); if(!length(ids)) next
    m<-m[match(ids,m$run_accession),]
    x<-raw[ids,,drop=FALSE]; rownames(x)<-m$tree_id
    parts[[comp]]<-x[rownames(x)%in%trees,,drop=FALSE] }
  sh<-Reduce(intersect,lapply(parts,rownames)); cols<-Reduce(union,lapply(parts,colnames))
  mk<-function(M){A<-matrix(0,length(sh),length(cols),dimnames=list(sh,cols)); A[sh,colnames(M)]<-M[sh,]; A}
  X<-rbind(mk(parts$needle),mk(parts$branch))
  meta<-data.frame(tree=rep(sh,2),comp=factor(rep(c("needle","branch"),each=length(sh))),stringsAsFactors=FALSE)
  keep<-rowSums(X)>=250; X<-X[keep,,drop=FALSE]; meta<-meta[keep,]
  ok<-names(which(table(meta$tree)==2)); sel<-meta$tree%in%ok
  X<-X[sel,,drop=FALSE]; meta<-meta[sel,]
  X<-X[,colSums(X>0)>=5 & colSums(X)/sum(X)>=1e-4,drop=FALSE]
  list(X=X,meta=meta,ev=env[match(meta$tree,env$tree_id),]) }
SPECS<-list(
 "climate3 (bio5,bio18,bio6)          "=c("bio5","bio18","bio6"),
 "literature7 (climate3 + soil4)      "=c("bio5","bio18","bio6","cec","ph","nitro","clay"),
 "RETIRED set (bio7,bio18,bio3,elev)  "=c("bio7","bio18","bio3","dem_elev_m"),
 "soil4 only (cec,ph,nitro,clay)      "=c("cec","ph","nitro","clay"))
kin<-function(raw,vars,perm=TRUE){
  b<-build(raw); if(ncol(b$X)<3) return(c(NA,NA))
  D<-dist(clr(b$X)); E<-scale(b$ev[,vars,drop=FALSE])
  Ks<-K[b$meta$tree,b$meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
  KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
  cp<-tryCatch(dbrda(D~KIN+Condition(cbind(E,b$meta$comp))),error=function(e) NULL)
  if(is.null(cp)||is.null(cp$CCA)) return(c(NA,NA))
  r2<-cp$CCA$tot.chi/cp$tot.chi
  if(!perm) return(c(r2,NA))
  h<-how(plots=Plots(strata=factor(b$meta$tree),type="free"),within=Within(type="none"),nperm=9999)
  c(r2, anova(cp,permutations=h)$`Pr(>F)`[1]) }
cat("host genetics | environment + compartment, family-propagated, 86 observations\n")
cat(sprintf("%-38s %8s %9s %10s %10s %10s\n","conditioning set","R2","p(perm)","null med","null 97.5","emp p"))
for (nm in names(SPECS)){
  v<-SPECS[[nm]]; o<-kin(FAMP,v)
  nulls<-replicate(40, kin(t(rowsum(t(stI), sample(lab))), v, perm=FALSE)[1])
  nulls<-nulls[is.finite(nulls)]
  ep<-(1+sum(nulls>=o[1]))/(1+length(nulls))
  cat(sprintf("%-38s %8.4f %9.4f %10.4f %10.4f %10.3f%s\n",nm,o[1],o[2],
      median(nulls),quantile(nulls,0.975),ep, if(ep<0.05) "  <-- exceeds null" else "")) }
cat("\nDONE\n"); sink()
