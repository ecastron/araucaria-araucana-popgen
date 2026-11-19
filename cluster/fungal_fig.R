#!/usr/bin/env Rscript
# Fungal figure: A dbRDA constrained on range (Bray-Curtis), B phylum stacked bars by site,
# C Shannon by range x compartment with Wilcoxon, D host genetics from the SNP matrix directly.
#
# NOTE ON A: "RDA on Bray-Curtis" is necessarily DISTANCE-BASED RDA -- Bray-Curtis is not
# Euclidean, so classical RDA cannot take it. vegan::dbrda is used throughout.
# NOTE ON B: ITS covers NINE stands (Coastal PN, TR; Andean BP, LM, LR, MC, PC, RA, TH).
# VA, HL and MM have no ITS libraries, so a 12-site panel is not possible.
# NOTE ON D: kinship estimates are replaced by the SNP matrix itself -- admix.bed is read
# directly, dosages are mean-imputed and centred, and principal components of that matrix are
# used as the constraint. No relatedness estimator is involved at any point.
.libPaths(c("~/R/library",.libPaths()))
suppressPackageStartupMessages({library(vegan); library(permute); library(genio)})
set.seed(2026); P<-"/home/ecastron/araucaria"
# ITS TABLE SELECTION (added 19 Nov 2026). "merged" = dada2_step1.R merged-pair table, the
# default, so every run before this date reproduces unchanged. "r1" = dada2_its_r1.R forward-read
# table, the pipeline the project chose on 4 Aug 2026 (ANALYSIS_SPEC.md: R1 keeps 84.7% of reads
# against 52.1%). Choose with the first command-line argument or the ITS_TABLE environment
# variable. OUT_DIR redirects output so an R1 run never overwrites the merged results.
ITS_TABLE<-{a<-commandArgs(trailingOnly=TRUE); if(length(a)) a[1] else Sys.getenv("ITS_TABLE","merged")}
stopifnot(ITS_TABLE %in% c("merged","r1"))
SFX<-if(ITS_TABLE=="r1") "_r1" else ""
cat(sprintf("ITS table: %s (seqtab_nochim%s.rds, taxonomy%s.rds)\n",ITS_TABLE,SFX,SFX))
OUT<-Sys.getenv("OUT_DIR",file.path(P,"results/fungfig")); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
sink(file.path(OUT,"fungal_fig.txt"),split=TRUE)
cat(sprintf("ITS table: %s\n",ITS_TABLE))

st<-readRDS(file.path(P,sprintf("microbiome/ITS/seqtab_nochim%s.rds",SFX)))
tax<-as.data.frame(readRDS(file.path(P,sprintf("microbiome/ITS/taxonomy%s.rds",SFX))),stringsAsFactors=FALSE)
man<-read.csv(file.path(P,"metadata/manifest microbiome runs.csv"),stringsAsFactors=FALSE,check.names=FALSE)
m<-man[match(rownames(st),man$run_accession),]
keep<-!is.na(m$tree_id) & m$compartment %in% c("needle","branch") & rowSums(st)>=1000
st<-st[keep,,drop=FALSE]; m<-m[keep,]
COASTAL<-c("PN","TR","VA")
m$range<-ifelse(m$site_code %in% COASTAL,"Coastal","Andean")
cat(sprintf("above-ground ITS libraries with tree_id and >=1000 reads: %d\n",nrow(st)))
print(table(m$site_code,m$compartment))
cat(sprintf("\nrange: %s\n\n",paste(names(table(m$range)),table(m$range),sep="=",collapse="  ")))

REL<-st/rowSums(st)

# TREE IS THE UNIT. Range is a property of a tree, and trees were sampled to span the
# distribution; stands are where the trees happen to be, not assigned treatments. Each tree
# therefore contributes ONE row, the mean of its relative abundances over the compartments
# sampled, which also removes the one non-independence nobody disputes -- needle and branch
# from the same individual. Stand is reported as a term, not used as the unit of analysis.
ntree<-rowsum(REL,m$tree_id)/as.vector(table(m$tree_id))
ntree<-ntree/rowSums(ntree)
mt0<-m[!duplicated(m$tree_id),]; mt0<-mt0[match(rownames(ntree),mt0$tree_id),]
cat(sprintf("\nTREE-LEVEL dataset: %d trees (%d Coastal, %d Andean) across %d stands\n\n",
    nrow(ntree),sum(mt0$range=="Coastal"),sum(mt0$range=="Andean"),length(unique(mt0$site_code))))
BC<-vegdist(REL,method="bray")          # library level, kept for the compartment contrast only
BCt<-vegdist(ntree,method="bray")       # tree level, used for every range test

# ---------------- A: dbRDA constrained on range ----------------
cat("=============== A. dbRDA, Bray-Curtis ~ range ===============\n")
# range alone is a 2-level factor, so it yields ONE constrained axis (rank = levels - 1) and
# there is no dbRDA2 to plot. Constraining on range + compartment gives 2 constrained axes,
# which is what the colour/shape scheme of this panel already encodes.
cap<-dbrda(BC~range+compartment,data=m)   # library level: compartment only exists here
an<-anova(cap,permutations=9999)
ant<-anova(cap,permutations=9999,by="terms")
R2<-cap$CCA$tot.chi/cap$tot.chi
cat(sprintf("constrained R2 (range + compartment) = %.4f  rank %d  F = %.3f  p = %.4f\n",
    R2,cap$CCA$rank,an$F[1],an$`Pr(>F)`[1]))
cat("  by term:\n"); print(ant)
capR<-dbrda(BC~range,data=m)
cat(sprintf("  range alone: R2 = %.4f (rank %d -- one axis only)\n",
    capR$CCA$tot.chi/capR$tot.chi,capR$CCA$rank))
# betadisper at BOTH levels. The range claim is made at TREE level, so pb -- written to
# panelA_stats.csv as disp_p and quoted in the legend -- must be the TREE-level value. The
# library-level figure is printed for context only. These differ (0.0021 vs 0.0071) and quoting
# the wrong one was a real inconsistency between the stats file and the documents.
bdL<-betadisper(BC,factor(m$range)); pbL<-permutest(bdL,permutations=9999)$tab
cat(sprintf("betadisper LIBRARY level F = %.3f p = %.4f (context only)\n",pbL[1,4],pbL[1,6]))
bd<-betadisper(BCt,factor(mt0$range)); pb<-permutest(bd,permutations=9999)$tab
mdc<-tapply(bd$distances,bd$group,mean)
cat(sprintf("betadisper TREE level: Coastal %.4f (n=%d), Andean %.4f (n=%d) -> smaller group is %s dispersed\n",
    mdc[["Coastal"]],sum(mt0$range=="Coastal"),mdc[["Andean"]],sum(mt0$range=="Andean"),
    if(mdc[["Coastal"]]>mdc[["Andean"]]) "MORE" else "LESS"))
cat(sprintf("betadisper TREE level F = %.3f p = %.4f -> %s\n",pb[1,4],pb[1,6],
    if(pb[1,6]<0.05) "dispersions differ; location effect is not clean" else "dispersions homogeneous"))
# also with stand as a stratum, since samples cluster within site
cat(sprintf("stand explains R2 = %.4f (p = %.4f)\n",
    {as<-adonis2(BC~site_code,data=m,permutations=9999); as$R2[1]},
    {as<-adonis2(BC~site_code,data=m,permutations=9999); as$`Pr(>F)`[1]}))
# The range contrast at the level the design intends: one row per tree.
at<-adonis2(BCt~range,data=mt0,permutations=9999)
ast<-adonis2(BCt~site_code,data=mt0,permutations=9999)
cat(sprintf("TREE-LEVEL range:  R2 = %.4f  F = %.3f  p = %.4f  (n = %d trees)\n",
    at$R2[1],at$F[1],at$`Pr(>F)`[1],nrow(ntree)))
cat(sprintf("TREE-LEVEL stand:  R2 = %.4f  F = %.3f  p = %.4f  (reported, not the unit)\n",
    ast$R2[1],ast$F[1],ast$`Pr(>F)`[1]))
ahp<-at$`Pr(>F)`[1]

sc<-scores(cap,display="sites",choices=1:2)
A<-data.frame(CAP1=sc[,1],CAP2=sc[,2],range=m$range,compartment=m$compartment,
              site=m$site_code,tree=m$tree_id)
ei<-cap$CCA$eig/cap$tot.chi
cat(sprintf("  dbRDA1 = %.2f%% of total inertia, dbRDA2 = %.2f%%\n",100*ei[1],100*ei[2]))
write.csv(A,file.path(OUT,"panelA_scores.csv"),row.names=FALSE)
write.csv(data.frame(R2=R2,F=an$F[1],p=an$`Pr(>F)`[1],
                     disp_p=pb[1,6],strata_p=ahp,
                     cap1_pct=100*ei[1],cap2_pct=100*ei[2],
                     tree_range_R2=at$R2[1],tree_range_p=at$`Pr(>F)`[1],
                     tree_stand_R2=ast$R2[1],tree_stand_p=ast$`Pr(>F)`[1],n_tree=nrow(ntree),
                     n_coastal=sum(mt0$range=="Coastal"),n_andean=sum(mt0$range=="Andean"),
                     range_R2=ant$SumOfSqs[1]/sum(ant$SumOfSqs),range_p=ant$`Pr(>F)`[1],
                     comp_R2=ant$SumOfSqs[2]/sum(ant$SumOfSqs),comp_p=ant$`Pr(>F)`[2]),
          file.path(OUT,"panelA_stats.csv"),row.names=FALSE)

# ---------------- B: phylum composition per site ----------------
cat("\n=============== B. ORDER composition by site ===============\n")
ph<-tax[colnames(st),"Order"]; ph<-sub("^o__","",ph); ph[is.na(ph)]<-"unassigned"
cat(sprintf("Order assigned for %d of %d ASVs (%.1f%%), %.1f%% of reads | %d orders\n",
    sum(ph!="unassigned"),length(ph),100*mean(ph!="unassigned"),
    100*sum(REL[,ph!="unassigned"])/sum(REL),length(unique(ph[ph!="unassigned"]))))
PHt<-t(rowsum(t(ntree),ph))
sitem<-rowsum(PHt,mt0$site_code)/as.vector(table(mt0$site_code))
sitem<-sitem/rowSums(sitem)
ord<-c("PN","TR",setdiff(rownames(sitem),c("PN","TR")))
sitem<-sitem[ord,,drop=FALSE]
keepph<-names(sort(colMeans(sitem),decreasing=TRUE))
top<-setdiff(keepph,"unassigned")[1:10]
top<-top[!is.na(top)]
oth<-setdiff(colnames(sitem),c(top,"unassigned"))
Bmat<-cbind(sitem[,top,drop=FALSE],
            `other assigned`=if(length(oth)) rowSums(sitem[,oth,drop=FALSE]) else 0,
            unassigned=if("unassigned" %in% colnames(sitem)) sitem[,"unassigned"] else 0)
print(round(100*Bmat,2))
write.csv(data.frame(site=rownames(Bmat),
                     range=ifelse(rownames(Bmat)%in%COASTAL,"Coastal","Andean"),
                     Bmat,check.names=FALSE),
          file.path(OUT,"panelB_phylum.csv"),row.names=FALSE)

# ---------------- C: Shannon by range x compartment ----------------
cat("\n=============== C. Shannon diversity ===============\n")
dep<-min(rowSums(st))
# RAREFY_DEPTH overrides the shallowest-library rule; used only to rerun the R1 table at the
# merged-table depth (1,075) as a sensitivity check. Unset, behavior is unchanged.
if(nzchar(Sys.getenv("RAREFY_DEPTH"))){dep<-as.integer(Sys.getenv("RAREFY_DEPTH")); stopifnot(dep<=min(rowSums(st)))}
cat(sprintf("rarefying to the shallowest retained library: %d reads\n",dep))
rar<-rrarefy(st,dep)
sh<-diversity(rar,index="shannon")
C<-data.frame(shannon=sh,range=m$range,compartment=m$compartment,
              site=m$site_code,tree=m$tree_id)
# range is tested on TREE means; compartment stays at library level because it varies within tree
Ct<-aggregate(shannon~tree+range+site,data=C,FUN=mean)
cat(sprintf("Shannon, tree level: %d trees\n",nrow(Ct)))
wt<-wilcox.test(Ct$shannon[Ct$range=="Coastal"],Ct$shannon[Ct$range=="Andean"])
cat(sprintf("  TREE-LEVEL range: Coastal (n=%d, med %.2f) vs Andean (n=%d, med %.2f): p = %.4f\n",
    sum(Ct$range=="Coastal"),median(Ct$shannon[Ct$range=="Coastal"]),
    sum(Ct$range=="Andean"),median(Ct$shannon[Ct$range=="Andean"]),wt$p.value))
write.csv(Ct,file.path(OUT,"panelC_shannon_tree.csv"),row.names=FALSE)
write.csv(C,file.path(OUT,"panelC_shannon.csv"),row.names=FALSE)
res<-list()
for (cp in c("needle","branch")){
  x<-C$shannon[C$compartment==cp & C$range=="Coastal"]
  y<-C$shannon[C$compartment==cp & C$range=="Andean"]
  w<-wilcox.test(x,y)
  cat(sprintf("  %-7s Coastal (n=%2d, med %.2f) vs Andean (n=%2d, med %.2f): W = %.0f, p = %.4f\n",
      cp,length(x),median(x),length(y),median(y),w$statistic,w$p.value))
  res[[length(res)+1]]<-data.frame(test=paste("range within",cp),n1=length(x),n2=length(y),
      med1=median(x),med2=median(y),W=w$statistic,p=w$p.value)
}
# compartment is a within-tree contrast, so it is tested PAIRED
pw<-merge(C[C$compartment=="needle",c("tree","shannon")],
          C[C$compartment=="branch",c("tree","shannon")],by="tree",suffixes=c(".n",".b"))
pw<-pw[!duplicated(pw$tree),]
w<-wilcox.test(pw$shannon.n,pw$shannon.b,paired=TRUE)
cat(sprintf("  needle vs branch, PAIRED within tree (n=%d pairs): V = %.0f, p = %.4f\n",
    nrow(pw),w$statistic,w$p.value))
res[[length(res)+1]]<-data.frame(test="compartment (paired)",n1=nrow(pw),n2=nrow(pw),
    med1=median(pw$shannon.n),med2=median(pw$shannon.b),W=w$statistic,p=w$p.value)
cat("\n  tree means within each compartment\n")
for (cp in c("needle","branch")){
  sm<-tapply(C$shannon[C$compartment==cp],C$tree[C$compartment==cp],mean)
  trg<-mt0$range[match(names(sm),mt0$tree_id)]
  cо<-sm[trg=="Coastal"]; an_<-sm[trg=="Andean"]
  w<-suppressWarnings(wilcox.test(cо,an_))
  cat(sprintf("  %-7s Coastal trees (n=%d, med %.2f) vs Andean (n=%d, med %.2f): p = %.4f\n",
      cp,length(cо),median(cо),length(an_),median(an_),w$p.value))
  res[[length(res)+1]]<-data.frame(test=paste("range within",cp,"(tree means)"),
      n1=length(cо),n2=length(an_),med1=median(cо),med2=median(an_),W=w$statistic,p=w$p.value)
}
res[[length(res)+1]]<-data.frame(test="range (tree level, both compartments)",
    n1=sum(Ct$range=="Coastal"),n2=sum(Ct$range=="Andean"),
    med1=median(Ct$shannon[Ct$range=="Coastal"]),med2=median(Ct$shannon[Ct$range=="Andean"]),
    W=wt$statistic,p=wt$p.value)
CT_<-do.call(rbind,res)
CT_$p_BH<-p.adjust(CT_$p,method="BH")
CT_$p_bonf<-p.adjust(CT_$p,method="bonferroni")
cat("\n  multiple-testing adjustment over the reported tests\n")
print(CT_[,c("test","p","p_BH","p_bonf")],row.names=FALSE)
write.csv(CT_,file.path(OUT,"panelC_tests.csv"),row.names=FALSE)

# ---------------- D: host genetics straight from the SNP matrix ----------------
cat("\n=============== D. host genetics, SNP matrix ===============\n")
# genio::read_plink is the standard PLINK reader. The earlier hand-written bit-unpacker is gone.
pl<-read_plink(file.path(P,"results/Q1/admix"),verbose=FALSE)
G<-t(pl$X)                                    # genio returns SNPs x individuals
cat(sprintf("SNP matrix: %d individuals x %d SNPs | missing %.1f%%\n",
    nrow(G),ncol(G),100*mean(is.na(G))))
radman<-read.csv(file.path(P,"metadata/manifest rad runs.csv"),stringsAsFactors=FALSE)
tof<-setNames(radman$tree_id,radman$run_accession)
tid<-tof[sub("\\.[12]$","",rownames(G))]
ok<-!is.na(tid); G<-G[ok,,drop=FALSE]; tid<-tid[ok]
G<-rowsum(G,tid,na.rm=TRUE)/rowsum(+(!is.na(G)),tid)   # average replicate libraries per tree
cat(sprintf("collapsed to %d trees\n",nrow(G)))
# The RAD matrix is 60.6% missing and per-TREE call rate spans 0.013 to 0.764. A PCA on a
# mean-imputed matrix that sparse can track WHICH sites a tree has data for -- i.e. library
# quality -- rather than genotype, which is the confound documented for this dataset elsewhere.
# Trees below 20% call rate are dropped, SNPs are then filtered at 50% (no threshold above 0.7
# retains a usable number), and PC-vs-call-rate correlations are reported so the confound is
# visible rather than assumed away.
trcr<-rowMeans(!is.na(G))
cat(sprintf("per-tree call rate: min %.3f median %.3f max %.3f | %d trees below 0.20\n",
    min(trcr),median(trcr),max(trcr),sum(trcr<0.20)))
G<-G[trcr>=0.20,,drop=FALSE]; trcr<-trcr[rownames(G)]
cr<-colMeans(!is.na(G))
G<-G[,cr>=0.5,drop=FALSE]
cat(sprintf("retained %d trees and %d SNPs at call rate >= 0.5 (missing now %.1f%%)\n",
    nrow(G),ncol(G),100*mean(is.na(G))))
mu<-colMeans(G,na.rm=TRUE)
G<-sweep(G,2,mu); G[is.na(G)]<-0                    # centre, then impute residual gaps at 0
sdv<-apply(G,2,sd); G<-G[,sdv>0,drop=FALSE]
pcG<-prcomp(G,center=FALSE,scale.=FALSE)
varex<-pcG$sdev^2/sum(pcG$sdev^2)
cat(sprintf("SNP PCA: PC1 %.1f%%  PC2 %.1f%%  PC3 %.1f%%  PC4 %.1f%%\n",
    100*varex[1],100*varex[2],100*varex[3],100*varex[4]))
cat("  CONFOUND CHECK -- correlation of each genetic PC with that tree's call rate:\n")
for (i in 1:6) cat(sprintf("    PC%-2d r = %+.3f%s\n",i,cor(pcG$x[,i],trcr[rownames(pcG$x)]),
    if(abs(cor(pcG$x[,i],trcr[rownames(pcG$x)]))>0.5) "   <-- tracks library quality" else ""))
# does the SNP PCA agree with the project PLINK PCA, which used the same genotypes?
ev<-tryCatch(read.table(file.path(P,"results/Q1/q1_pca_plink.eigenvec"),stringsAsFactors=FALSE),
             error=function(e) NULL)
if(!is.null(ev)){
  evid<-tof[sub("\\.[12]$","",ev$V1)]; ok2<-!is.na(evid)
  pl<-tapply(ev$V3[ok2],evid[ok2],mean); sh0<-intersect(names(pl),rownames(pcG$x))
  cat(sprintf("  agreement with the PLINK PCA on PC1 (%d shared trees): |r| = %.3f\n",
      length(sh0),abs(cor(pl[sh0],pcG$x[sh0,1]))))
}
shared<-intersect(rownames(pcG$x),m$tree_id)
cat(sprintf("trees with BOTH fungal and SNP data: %d\n",length(shared)))
sel<-m$tree_id %in% shared
Dm<-vegdist(REL[sel,,drop=FALSE],method="bray"); md<-m[sel,]
# Genetic PCs in a spatially structured population are partly a PROXY FOR STAND: trees in the
# same stand are genetically similar, and stand alone explains 11% of fungal composition here.
# A test of genetic PCs that does not hold stand constant therefore restates "stands differ".
# Samples are also paired within tree. Both are handled below: the model ladder adds stand, and
# whole trees are the permutation unit. Trees are restricted to those with BOTH compartments so
# the strata are balanced.
# Requiring both compartments left only 22 trees, and at that n the positive controls fail, so
# the test is uninformative. Instead each TREE is represented once, by the mean of its relative
# abundances across the compartments sampled. That removes the within-tree pseudoreplication
# entirely, makes every row independent, and keeps all trees that have both data types.
relT<-rowsum(REL[sel,,drop=FALSE],md$tree_id)/as.vector(table(md$tree_id))
relT<-relT/rowSums(relT)
mt<-md[!duplicated(md$tree_id),]; mt<-mt[match(rownames(relT),mt$tree_id),]
D2<-vegdist(relT,method="bray")
cat(sprintf("\n  tree-level design: %d trees, one row each (no pseudoreplication)\n",nrow(relT)))
cat(sprintf("  stands represented: %s\n",paste(sort(unique(mt$site_code)),collapse=" ")))
rows<-list()
LAD<-list("genetic PCs, no conditioning"="dbrda(D2~X)",
          "genetic PCs | stand"         ="dbrda(D2~X+Condition(mt$site_code))",
          "genetic PCs | range"         ="dbrda(D2~X+Condition(mt$range))")
for (k in c(2,4,8)){
  X<-pcG$x[mt$tree_id,1:k,drop=FALSE]
  for (nm in names(LAD)){
    cp<-tryCatch(eval(parse(text=LAD[[nm]])),error=function(e) NULL)
    if(is.null(cp)||is.null(cp$CCA)||cp$CCA$rank<1){cat(sprintf("  k=%-2d %-30s not identifiable\n",k,nm)); next}
    a<-anova(cp,permutations=9999); r<-cp$CCA$tot.chi/cp$tot.chi
    cat(sprintf("  k=%-2d %-30s R2 = %.4f  F = %.3f  p = %.4f\n",k,nm,r,a$F[1],a$`Pr(>F)`[1]))
    rows[[length(rows)+1]]<-data.frame(k=k,model=nm,R2=r,F=a$F[1],p=a$`Pr(>F)`[1]) } }
cat("\n  how much of each genetic axis is simply stand identity\n")
for (i in 1:4) cat(sprintf("    PC%-2d ~ stand  R2 = %.3f\n",i,
    summary(lm(pcG$x[mt$tree_id,i]~factor(mt$site_code)))$r.squared))
cat("\n  positive controls on the same tree-level data\n")
for (trm in c("range","site_code")){
  aa<-adonis2(as.formula(paste("D2 ~",trm)),data=mt,permutations=9999)
  cat(sprintf("  %-11s R2 = %.4f  p = %.4f\n",trm,aa$R2[1],aa$`Pr(>F)`[1]))
  rows[[length(rows)+1]]<-data.frame(k=NA,model=trm,R2=aa$R2[1],F=aa$F[1],p=aa$`Pr(>F)`[1]) }
D<-do.call(rbind,rows)
write.csv(D,file.path(OUT,"panelD_results.csv"),row.names=FALSE)
write.csv(data.frame(tree=rownames(pcG$x),PC1=pcG$x[,1],PC2=pcG$x[,2]),
          file.path(OUT,"panelD_snp_pca.csv"),row.names=FALSE)
cat("\nDONE\n"); sink()
