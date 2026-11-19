#!/usr/bin/env Rscript
# Figure 1C, recomputed so all four bars come from ONE analysis.
#
# The published panel mixed sources: stand at ASV level, climate from the RETIRED 4-variable set
# (bio7/bio18/bio3/elevation, dropped 31 Aug because elevation is redundant with bio6 at |r| =
# 0.91), soil from the final literature set, and host genetics at family-propagated resolution.
# Four bars, three analyses.
#
# Fixed here at FAMILY-PROPAGATED resolution throughout, because the shuffled-label null that
# makes the host-genetics bar readable only exists for an agglomerated feature set -- there are
# no labels to permute at ASV level. Conditioning uses the final 7-variable literature set and
# its climate(3)/soil(4) split everywhere.
suppressPackageStartupMessages({library(vegan); library(permute)})
set.seed(2026); P <- "/home/ecastron/araucaria"
dir.create(file.path(P,"results/F1C"), recursive=TRUE, showWarnings=FALSE)
sink(file.path(P,"results/F1C/f1c_family.txt"), split=TRUE)

stI <- readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))
tax <- readRDS(file.path(P,"microbiome/ITS/taxonomy.rds"))
stopifnot(identical(rownames(tax), colnames(stI)))

# rank-propagated family agglomeration: an ASV with no Family inherits its finest assigned rank,
# so no reads are discarded (matches the 26 Aug "propagated" variant exactly)
fam <- tax[,"Family"]; lab <- as.character(fam)
for (r in c("Order","Class","Phylum","Kingdom")){
  need <- is.na(lab); if(!any(need)) break
  src <- as.character(tax[,r]); ok <- need & !is.na(src)
  lab[ok] <- paste0(tolower(substr(r,1,1)), "__", src[ok], "_unclassified")
}
lab[is.na(lab)] <- "k__Fungi_unclassified"
FAMP <- t(rowsum(t(stI), lab))
cat(sprintf("family-propagated: %d bins from %d ASVs, 100%% of reads\n\n", ncol(FAMP), ncol(stI)))

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
env<-ch[,c("tree_id","site_code","lon","lat","bio5","bio18","bio6")]
so<-read.delim(file.path(P,"metadata/soil_per_stand.tsv"),check.names=FALSE)
env<-merge(env, so[,c("site_code","clay_0-5cm","cec_0-5cm","phh2o_0-5cm","nitrogen_0-5cm")],
           by="site_code", all.x=TRUE)
names(env)[names(env)=="clay_0-5cm"]<-"clay"; names(env)[names(env)=="cec_0-5cm"]<-"cec"
names(env)[names(env)=="phh2o_0-5cm"]<-"ph";  names(env)[names(env)=="nitrogen_0-5cm"]<-"nitro"

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
  X<-X[,colSums(X>0)>=5 & colSums(X)/sum(X)>=1e-4,drop=FALSE]
  meta$cord<-factor(ifelse(cord_of[meta$tree]=="Costa","Coastal","Andean"))
  list(X=X, meta=meta, ev=env[match(meta$tree,env$tree_id),])
}

b<-build(FAMP); X<-b$X; meta<-b$meta; ev<-b$ev
cat(sprintf("%d trees x 2 compartments = %d observations | %d features | %.1f%% zeros\n\n",
    length(unique(meta$tree)), nrow(X), ncol(X), 100*mean(X==0)))
D<-dist(clr(X))
hB<-how(plots=Plots(strata=factor(meta$tree),type="free"),within=Within(type="none"),nperm=9999)
hW<-how(plots=Plots(strata=factor(meta$tree),type="none"),within=Within(type="free"),nperm=9999)
LIT<-scale(ev[,c("bio5","bio18","bio6","cec","ph","nitro","clay")])
CL<-LIT[,c("bio5","bio18","bio6")]; SO<-LIT[,c("cec","ph","nitro","clay")]
CD<-model.matrix(~meta$cord)[,-1,drop=FALSE]
Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
f<-function(expr){cp<-eval(expr); an<-anova(cp,permutations=hB)
  c(cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$`Pr(>F)`[1])}

cat("=========== THE FOUR BARS, one analysis, final variable set ===========\n")
aS<-adonis2(D~site_code,data=data.frame(site_code=factor(ev$site_code)),permutations=hB)
cat(sprintf("  stand                                R2 = %.4f  F = %.3f  p = %.4f\n",aS$R2[1],aS$F[1],aS$`Pr(>F)`[1]))
r<-f(quote(dbrda(D~CL+Condition(CD))))
cat(sprintf("  climate | cordillera                 R2 = %.4f  F = %.3f  p = %.4f\n",r[1],r[2],r[3]))
r<-f(quote(dbrda(D~SO+Condition(cbind(CL,CD)))))
cat(sprintf("  soil | climate + cordillera          R2 = %.4f  F = %.3f  p = %.4f\n",r[1],r[2],r[3]))
r<-f(quote(dbrda(D~KIN+Condition(cbind(CL,meta$comp)))))
cat(sprintf("  host genetics | climate + compartment R2 = %.4f  F = %.3f  p = %.4f\n",r[1],r[2],r[3]))
OBS<-r[1]
cat("\n  (context, not plotted)\n")
aP<-adonis2(D~comp,data=meta,permutations=hW)
cat(sprintf("  compartment (within tree)            R2 = %.4f  p = %.4f\n",aP$R2[1],aP$`Pr(>F)`[1]))
aC<-adonis2(D~cord,data=meta,permutations=hB)
cat(sprintf("  cordillera                           R2 = %.4f  p = %.4f\n",aC$R2[1],aC$`Pr(>F)`[1]))

cat("\n=========== SHUFFLED-LABEL NULL for the host-genetics bar ===========\n")
kinfit <- function(raw){
  bb<-build(raw); if(ncol(bb$X)<3) return(NA_real_)
  DD<-dist(clr(bb$X)); ee<-bb$ev
  CLx<-scale(ee[,c("bio5","bio18","bio6")])
  Kx<-K[bb$meta$tree,bb$meta$tree]; Kx[is.na(Kx)]<-mean(Kx[upper.tri(Kx)],na.rm=TRUE); diag(Kx)<-max(Kx,na.rm=TRUE)
  KN<-cmdscale(as.dist(max(Kx)-Kx),k=3)
  tryCatch({cp<-dbrda(DD~KN+Condition(cbind(CLx,bb$meta$comp))); cp$CCA$tot.chi/cp$tot.chi},
           error=function(e) NA_real_)
}
NP<-50; nullR<-rep(NA_real_,NP)
for (i in seq_len(NP)){
  M<-t(rowsum(t(stI), sample(lab)))
  nullR[i]<-kinfit(M)
  if(i%%10==0) cat(sprintf("  null %d/%d  median so far %.4f\n",i,NP,median(nullR[1:i],na.rm=TRUE)))
}
nullR<-nullR[is.finite(nullR)]
cat(sprintf("\n  null R2 over %d shuffles: median %.4f, 95%% range [%.4f, %.4f]\n",
    length(nullR),median(nullR),quantile(nullR,0.025),quantile(nullR,0.975)))
cat(sprintf("  observed %.4f -> empirical p = %.4f\n",OBS,(1+sum(nullR>=OBS))/(1+length(nullR))))
write.table(data.frame(null_r2=nullR),file.path(P,"results/F1C/null_family.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE)
cat("\nDONE\n"); sink()
