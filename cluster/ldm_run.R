#!/usr/bin/env Rscript
# LDM (Hu & Satten 2020) as a formal alternative to PERMANOVA for the range effect.
# Run at BOTH levels: tree level (the agreed unit) and library level with tree as a cluster,
# which is what LDM's cluster.id argument is for.
.libPaths(c("~/R/library",.libPaths()))
suppressPackageStartupMessages({library(vegan); library(LDM)})
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
OUT<-Sys.getenv("OUT_DIR",file.path(P,"results/fungfig")); dir.create(OUT,recursive=TRUE,showWarnings=FALSE); sink(file.path(OUT,"ldm.txt"),split=TRUE)
cat(sprintf("ITS table: %s\n",ITS_TABLE))
cat(sprintf("LDM version %s\n\n",as.character(packageVersion("LDM"))))

st<-readRDS(file.path(P,sprintf("microbiome/ITS/seqtab_nochim%s.rds",SFX)))
man<-read.csv(file.path(P,"metadata/manifest microbiome runs.csv"),stringsAsFactors=FALSE,check.names=FALSE)
m<-man[match(rownames(st),man$run_accession),]
keep<-!is.na(m$tree_id)&m$compartment%in%c("needle","branch")&rowSums(st)>=1000
st<-st[keep,,drop=FALSE]; m<-m[keep,]
m$range<-factor(ifelse(m$site_code%in%c("PN","TR","VA"),"Coastal","Andean"))
m$compartment<-factor(m$compartment); m$site_code<-factor(m$site_code)
st<-st[,colSums(st>0)>=3,drop=FALSE]
cat(sprintf("%d libraries x %d ASVs (prevalence >= 3)\n",nrow(st),ncol(st)))

# ---------- tree level ----------
REL<-st/rowSums(st)
cnt<-rowsum(st,m$tree_id)                       # summed counts per tree
mt<-m[!duplicated(m$tree_id),]; mt<-mt[match(rownames(cnt),mt$tree_id),]
mt<-droplevels(mt)
cat(sprintf("\nTREE LEVEL: %d trees (%d Coastal, %d Andean)\n",
    nrow(cnt),sum(mt$range=="Coastal"),sum(mt$range=="Andean")))
r1<-ldm(formula=cnt~range,data=mt,dist.method="bray",n.perm.max=10000,seed=2026)
cat("  LDM global test, range:\n")
cat(sprintf("    p.global.freq = %s | p.global.tran = %s | p.global.omni = %s\n",
    format(r1$p.global.freq,digits=4),format(r1$p.global.tran,digits=4),
    format(r1$p.global.omni,digits=4)))
cat(sprintf("    VE (freq) = %s\n",paste(format(r1$VE.global.freq.submodels,digits=4),collapse=", ")))
nd<-tryCatch(sum(r1$q.otu.omni<=0.1,na.rm=TRUE),error=function(e) NA)
cat(sprintf("    ASVs with omnibus q <= 0.10: %s\n",nd))
p2<-permanovaFL(formula=cnt~range,data=mt,dist.method="bray",n.perm.max=10000,seed=2026)
cat(sprintf("  permanovaFL (LDM's PERMANOVA), range: p = %s\n",format(p2$p.permanova,digits=4)))

# ---------- library level, tree as cluster ----------
cat(sprintf("\nLIBRARY LEVEL with tree as cluster: %d libraries in %d trees\n",
    nrow(st),length(unique(m$tree_id))))
r2<-tryCatch(ldm(formula=st~range+compartment,data=m,dist.method="bray",
                 cluster.id=m$tree_id,n.perm.max=10000,seed=2026),
             error=function(e){cat("  ldm clustered failed:",conditionMessage(e),"\n"); NULL})
if(!is.null(r2)){
  cat("  LDM global test, range then compartment (tree as cluster):\n")
  cat(sprintf("    p.global.omni = %s\n",paste(format(r2$p.global.omni,digits=4),collapse=", ")))
}
p3<-tryCatch(permanovaFL(formula=st~range+compartment,data=m,dist.method="bray",
                         cluster.id=m$tree_id,n.perm.max=10000,seed=2026),
             error=function(e){cat("  permanovaFL clustered failed:",conditionMessage(e),"\n"); NULL})
if(!is.null(p3)) cat(sprintf("  permanovaFL clustered: p = %s\n",
                             paste(format(p3$p.permanova,digits=4),collapse=", ")))
cat("\nDONE\n"); sink()
