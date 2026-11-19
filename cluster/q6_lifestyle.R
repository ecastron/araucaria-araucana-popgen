#!/usr/bin/env Rscript
# ITS recoded to FUNCTIONAL LIFESTYLE (FungalTraits v1.2, Polme et al. 2020) vs ASV level,
# same paired needle+branch design throughout.
#
# HYPOTHESIS, declared before running (ECN, 4 Sep 2026): host genotype may select on what a
# fungus DOES rather than on what it IS. This is NOT the family-agglomeration test of 26 Aug --
# that merged WITHIN clades and left kinship R2 unchanged (0.0341/0.0340/0.0335). Lifestyle cuts
# ACROSS clades: two congeners can be an asymptomatic endophyte and a necrotroph, and no
# taxonomic rank will ever separate them. It is an orthogonal projection, not a coarser one.
#
# REPORTED WHICHEVER WAY IT LANDS. The taxonomic null is printed alongside, not replaced.
#
# Coverage measured beforehand: primary_lifestyle reaches 2,137/4,061 ASVs (52.6%) carrying
# 68.1% of reads. The 31.9% of reads with no lifestyle are NOT a random sample -- they are
# led by Teratosphaeriaceae, Coryneliaceae, Phacidiaceae and Tympanidaceae, southern-conifer
# associated lineages that are exactly where a host-specific effect would be expected to sit.
# The propagated variant below keeps them as their own bins for that reason.
suppressPackageStartupMessages({library(vegan); library(permute)})
set.seed(2026); P <- "/home/ecastron/araucaria"
dir.create(file.path(P,"results/Qlife"), recursive=TRUE, showWarnings=FALSE)
sink(file.path(P,"results/Qlife/q6_lifestyle.txt"), split=TRUE)

stI <- readRDS(file.path(P,"microbiome/ITS/seqtab_nochim.rds"))
tax <- readRDS(file.path(P,"microbiome/ITS/taxonomy.rds"))
stopifnot(identical(rownames(tax), colnames(stI)))
ft <- read.delim(file.path(P,"refs/fungaltraits_genus.tsv"), stringsAsFactors=FALSE, quote="")
ft$genus <- trimws(ft$genus); ft$primary_lifestyle[ft$primary_lifestyle==""] <- NA

gc_ <- trimws(sub("^g__","",as.character(tax[,"Genus"])))
gc_[gc_ %in% c("","NA","unidentified")] <- NA
life <- ft$primary_lifestyle[match(gc_, ft$genus)]
rd <- colSums(stI); tot <- sum(rd); has <- !is.na(life)

cat("=================== LIFESTYLE RECODING ===================\n")
cat(sprintf("ITS ASV table   : %d samples x %d ASVs, %s reads\n", nrow(stI), ncol(stI), format(tot,big.mark=",")))
cat(sprintf("genus assigned  : %d ASVs (%.1f%%), %.1f%% of reads\n", sum(!is.na(gc_)),100*mean(!is.na(gc_)),100*sum(rd[!is.na(gc_)])/tot))
cat(sprintf("lifestyle known : %d ASVs (%.1f%%), %.1f%% of reads -> %d lifestyles\n\n",
    sum(has), 100*mean(has), 100*sum(rd[has])/tot, length(unique(life[has]))))

# (a) STRICT: unassigned dropped
LIFE <- t(rowsum(t(stI[,has,drop=FALSE]), life[has]))
cat(sprintf("STRICT     : %d bins, drops %d ASVs carrying %.1f%% of reads\n",
    ncol(LIFE), sum(!has), 100*sum(rd[!has])/tot))

# (b) PROPAGATED: an ASV with no lifestyle keeps its finest taxonomic rank as its own bin,
# so no reads are discarded and the conifer-associated unassigned fraction stays in the model.
lab <- life
fin <- rep("k__Fungi", ncol(stI))
for (r in c("Kingdom","Phylum","Class","Order","Family")) if (r %in% colnames(tax)) {
  v <- as.character(tax[,r]); ok <- !is.na(v) & v!=""
  fin[ok] <- paste0(substr(tolower(r),1,1),"__",sub("^[kpcof]__","",v[ok]),"_nolifestyle") }
lab[!has] <- fin[!has]
LIFEP <- t(rowsum(t(stI), lab))
cat(sprintf("PROPAGATED : %d bins, keeps all ASVs and 100%% of reads (%d are _nolifestyle bins, %.1f%% of reads)\n",
    ncol(LIFEP), sum(grepl("_nolifestyle$",colnames(LIFEP))),
    100*sum(LIFEP[,grepl("_nolifestyle$",colnames(LIFEP)),drop=FALSE])/tot))

# (c) COARSE trophic mode -- the classic FUNGuild-style 4-way axis. Declared here rather than
# added later, and reported whatever it shows.
tro <- rep(NA_character_, length(life))
tro[grepl("saprotroph|sooty_mold", life)] <- "saprotroph"
tro[grepl("pathogen|parasite", life)] <- "pathotroph"
tro[grepl("mycorrhizal|endophyte|lichenized|symbiont|epiphyte", life)] <- "symbiotroph"
tro[has & is.na(tro)] <- "other"
TRO <- t(rowsum(t(stI[,!is.na(tro),drop=FALSE]), tro[!is.na(tro)]))
cat(sprintf("COARSE     : %d trophic modes (%s)\n\n", ncol(TRO), paste(colnames(TRO),collapse=", ")))

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
  nbefore<-ncol(X)
  X<-X[,colSums(X>0)>=5 & colSums(X)/sum(X)>=1e-4,drop=FALSE]
  meta$cord<-factor(ifelse(cord_of[meta$tree]=="Costa","Coastal","Andean"))
  ev<-env[match(meta$tree,env$tree_id),]
  list(X=X, meta=meta, ev=ev, nbefore=nbefore)
}

run <- function(raw, tag){
  b<-build(raw); X<-b$X; meta<-b$meta; ev<-b$ev
  cat(sprintf("\n############### %s ###############\n", tag))
  cat(sprintf("  %d trees x 2 compartments = %d observations | %d features (of %d before the filter) | %.1f%% zeros\n",
      length(unique(meta$tree)), nrow(X), ncol(X), b$nbefore, 100*mean(X==0)))
  if (ncol(X) < 3){cat("  <3 features survive; not fitted\n"); return(invisible())}
  D<-dist(clr(X))
  hB<-how(plots=Plots(strata=factor(meta$tree),type="free"),within=Within(type="none"),nperm=9999)
  hW<-how(plots=Plots(strata=factor(meta$tree),type="none"),within=Within(type="free"),nperm=9999)
  LIT<-scale(ev[,c("bio5","bio18","bio6","cec","ph","nitro","clay")])
  colnames(LIT)<-c("bio5","bio18","bio6","cec","ph","nitro","clay")
  CL<-LIT[,c("bio5","bio18","bio6")]; SO<-LIT[,c("cec","ph","nitro","clay")]
  Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
  KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
  CD<-model.matrix(~meta$cord)[,-1,drop=FALSE]
  f<-function(expr) tryCatch({cp<-eval(expr)
        if(is.null(cp$CCA)||cp$CCA$rank<1) return("not identifiable (constrained rank 0)")
        an<-anova(cp,permutations=hB)
        sprintf("R2 = %.4f  F = %.3f  p = %.4f", cp$CCA$tot.chi/cp$tot.chi, an$F[1], an$`Pr(>F)`[1])},
        error=function(e) paste("failed:",conditionMessage(e)))
  cat("\n  --- POSITIVE CONTROLS (must survive, or the recoding destroyed information) ---\n")
  aS<-adonis2(D~site_code,data=data.frame(site_code=factor(ev$site_code)),permutations=hB)
  cat(sprintf("   STAND                          R2 = %.4f  F = %.3f  p = %.4f\n",aS$R2[1],aS$F[1],aS$`Pr(>F)`[1]))
  aP<-adonis2(D~comp,data=meta,permutations=hW)
  cat(sprintf("   COMPARTMENT (within tree)      R2 = %.4f  F = %.3f  p = %.4f\n",aP$R2[1],aP$F[1],aP$`Pr(>F)`[1]))
  cat(sprintf("   climate | cordillera           %s\n", f(quote(dbrda(D~CL+Condition(CD))))))
  cat(sprintf("   soil(4) | climate(3) + cord    %s\n", f(quote(dbrda(D~SO+Condition(cbind(CL,CD)))))))
  aC<-adonis2(D~cord,data=meta,permutations=hB)
  cat(sprintf("   CORDILLERA                     R2 = %.4f  F = %.3f  p = %.4f\n",aC$R2[1],aC$F[1],aC$`Pr(>F)`[1]))
  cat("\n  --- THE TEST ---\n")
  cat(sprintf("   HOST GENETICS kin | clim + comp        %s\n", f(quote(dbrda(D~KIN+Condition(cbind(CL,meta$comp)))))))
  cat(sprintf("   HOST GENETICS kin | lit7 + comp        %s\n", f(quote(dbrda(D~KIN+Condition(cbind(LIT,meta$comp)))))))
  cat(sprintf("   HOST GENETICS kin | clim + comp + cord %s\n", f(quote(dbrda(D~KIN+Condition(cbind(CL,meta$comp,CD)))))))
  invisible(NULL)
}

run(stI,  "ITS at ASV level (reference -- reproduces the published null)")
run(LIFE, "LIFESTYLE, strict (unassigned dropped)")
run(LIFEP,"LIFESTYLE, propagated (unassigned kept as own bins)")
run(TRO,  "COARSE trophic mode (saprotroph / pathotroph / symbiotroph / other)")

# ---------------- NULL CONTROL ----------------
# Recoding collapses ~4,000 features to ~20 against the same 86 observations and the same 3
# kinship axes. That changes dimensionality and sparsity, which moves the statistic on its own.
# The 26 Aug family run showed the ~3.4% kinship R2 is the FLOOR from fitting 3 axes to 86 obs,
# not a weak signal. Shuffling lifestyle labels across ASVs preserves bin sizes exactly --
# hence dimensionality and sparsity -- while destroying functional coherence.
cat("\n\n=================== NULL: SHUFFLED LIFESTYLE LABELS ===================\n")
kinfit <- function(raw){
  b<-build(raw); X<-b$X; meta<-b$meta; ev<-b$ev
  if (ncol(X)<3) return(c(NA,NA))
  D<-dist(clr(X))
  CL<-scale(ev[,c("bio5","bio18","bio6")])
  Ks<-K[meta$tree,meta$tree]; Ks[is.na(Ks)]<-mean(Ks[upper.tri(Ks)],na.rm=TRUE); diag(Ks)<-max(Ks,na.rm=TRUE)
  KIN<-cmdscale(as.dist(max(Ks)-Ks),k=3)
  hB<-how(plots=Plots(strata=factor(meta$tree),type="free"),within=Within(type="none"),nperm=999)
  tryCatch({cp<-dbrda(D~KIN+Condition(cbind(CL,meta$comp)))
            if(is.null(cp$CCA)||cp$CCA$rank<1) return(c(NA,NA))
            a<-anova(cp,permutations=hB); c(cp$CCA$tot.chi/cp$tot.chi, a$`Pr(>F)`[1])},
           error=function(e) c(NA,NA))
}
for (nm in c("strict","propagated")) {
  obs <- kinfit(if (nm=="strict") LIFE else LIFEP)
  cat(sprintf("\nOBSERVED %s : kinship R2 = %.4f  p = %.4f\n", nm, obs[1], obs[2]))
  NP<-50; nullR<-rep(NA_real_,NP)
  for (i in seq_len(NP)){
    if (nm=="strict") { sh<-sample(life[has]); M<-t(rowsum(t(stI[,has,drop=FALSE]), sh)) }
    else              { sh<-sample(lab);       M<-t(rowsum(t(stI), sh)) }
    nullR[i]<-kinfit(M)[1]
    if (i%%10==0) cat(sprintf("   null %d/%d  median so far %.4f\n", i, NP, median(nullR[1:i],na.rm=TRUE)))
  }
  nullR<-nullR[is.finite(nullR)]
  ep<-(1+sum(nullR>=obs[1]))/(1+length(nullR))
  cat(sprintf("  null R2 over %d shuffles: median %.4f, 95%% range [%.4f, %.4f]\n",
      length(nullR), median(nullR), quantile(nullR,0.025), quantile(nullR,0.975)))
  cat(sprintf("  observed %.4f -> empirical p = %.4f  (%s)\n", obs[1], ep,
      if (ep<0.05) "EXCEEDS the same-dimensionality artefact -- functional coherence matters"
      else "NOT distinguishable from a same-dimensionality artefact"))
  write.table(data.frame(variant=nm, null_r2=nullR), file.path(P,sprintf("results/Qlife/null_%s.tsv",nm)),
              sep="\t", row.names=FALSE, quote=FALSE)
}
cat("\nDONE\n")
sink()
