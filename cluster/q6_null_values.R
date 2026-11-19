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
sink(file.path(P,"results/Q6/family/null_rerun.txt"), split=TRUE)

stI <- readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))
tax <- readRDS(file.path(P,"microbiome/ITS/taxonomy.rds"))
stopifnot(identical(rownames(tax), colnames(stI)))

fam <- tax[,"Family"]
cat("=================== FAMILY AGGLOMERATION ===================\n")
cat(sprintf("ITS ASV table : %d samples x %d ASVs, %s reads\n", nrow(stI), ncol(stI), format(sum(stI),big.mark=",")))
keepf <- !is.na(fam)
cat(sprintf("Family assigned: %d/%d ASVs (%.1f%%), %.1f%% of reads\n",
    sum(keepf), length(fam), 100*mean(keepf), 100*sum(stI[,keepf])/sum(stI)))

# (a) STRICT: drop ASVs with no Family, as tax_glom(NArm=TRUE) would
cat(sprintf("STRICT  drops %d ASVs carrying %.1f%% of reads\n",
    sum(!keepf), 100*sum(stI[,!keepf])/sum(stI)))
FAM <- t(rowsum(t(stI[,keepf,drop=FALSE]), fam[keepf]))
cat(sprintf("  -> %d families\n", ncol(FAM)))

# (b) PROPAGATED: an ASV with no Family inherits its finest assigned rank as a pseudo-family
# (o__Helotiales_unclassified etc.), so no reads are discarded. These bins are heterogeneous --
# one pseudo-family may span several true families -- so they are coarser than a real Family,
# and that is reported rather than hidden.
lab <- as.character(fam)
for (r in c("Order","Class","Phylum","Kingdom")){
  need <- is.na(lab); if(!any(need)) break
  src <- as.character(tax[,r]); ok <- need & !is.na(src)
  lab[ok] <- paste0(tolower(substr(r,1,1)), "__", src[ok], "_unclassified")
}
lab[is.na(lab)] <- "k__Fungi_unclassified"
FAMP <- t(rowsum(t(stI), lab))
npseudo <- sum(grepl("_unclassified$", colnames(FAMP)))
pseudo_reads <- 100*sum(FAMP[,grepl("_unclassified$",colnames(FAMP)),drop=FALSE])/sum(FAMP)
cat(sprintf("PROPAGATED keeps all %d ASVs and 100%% of reads -> %d bins\n", ncol(stI), ncol(FAMP)))
cat(sprintf("  of which %d are rank-propagated pseudo-families carrying %.1f%% of reads\n\n",
    npseudo, pseudo_reads))

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
#rA <- run(stI,"ITS at ASV level (reference)")
#rF <- run(FAM,"ITS at FAMILY level (strict, unassigned dropped)")
#rP <- run(FAMP,"ITS at FAMILY level (rank-propagated, all reads kept)")

# ---- NULL CONTROL -------------------------------------------------------------
# Agglomeration cuts the feature count ~10x and sharply reduces sparsity, which raises power on
# its own. To tell a biological effect from a resolution artefact, shuffle the family labels
# across ASVs: group SIZES (hence dimensionality and sparsity) are preserved exactly, taxonomic
# coherence is destroyed. If the real signal is not stronger than this null, it is an artefact.
cat("\n\n=================== NULL: SHUFFLED AGGLOMERATION ===================\n")
cat("family labels permuted across ASVs; group sizes preserved, taxonomy destroyed\n")
realfit <- function(raw){
  b<-build(raw); X<-b$X; meta<-b$meta; ev<-b$ev
  D<-dist(clr(X))
  CLIM<-scale(ev[,c("bio7","bio18","bio3","dem_elev_m")])
  Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
  KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
  hB<-how(plots=Plots(strata=factor(meta$tree),type="free"),within=Within(type="none"),nperm=999)
  tryCatch({cp<-dbrda(D~KIN+Condition(CLIM+meta$comp)); a<-anova(cp,permutations=hB)
            c(cp$CCA$tot.chi/cp$tot.chi, a$`Pr(>F)`[1])}, error=function(e) c(NA,NA))
}
obsF <- realfit(FAMP)
cat(sprintf("\nOBSERVED family-level host genetics: R2=%.4f p=%.4f\n", obsF[1], obsF[2]))
NPERM <- 50
nullR <- numeric(NPERM)
for (i in seq_len(NPERM)){
  shuf <- sample(lab)
  Fi <- t(rowsum(t(stI), shuf))
  v <- realfit(Fi); nullR[i] <- v[1]
  if (i %% 10 == 0) cat(sprintf("  null %d/%d  median R2 so far %.4f\n", i, NPERM, median(nullR[1:i],na.rm=TRUE)))
}
nullR <- nullR[is.finite(nullR)]
write.csv(data.frame(shuffle=seq_along(nullR), R2=nullR),
          file.path(P,"results/Q6/family/null_shuffles.csv"), row.names=FALSE)
cat(sprintf("\nwrote %d individual shuffle values -> results/Q6/family/null_shuffles.csv\n", length(nullR)))
cat(sprintf("\nnull R2 over %d shuffles: median %.4f, 95%% range [%.4f, %.4f]\n",
    length(nullR), median(nullR), quantile(nullR,0.025), quantile(nullR,0.975)))
cat(sprintf("observed R2 = %.4f -> empirical p = %.4f  (%s)\n", obsF[1],
    (1+sum(nullR>=obsF[1]))/(1+length(nullR)),
    if ((1+sum(nullR>=obsF[1]))/(1+length(nullR)) < 0.05)
      "exceeds the resolution artefact -- taxonomic coherence matters"
    else "NOT distinguishable from a same-dimensionality artefact"))

