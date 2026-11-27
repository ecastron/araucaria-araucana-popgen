.libPaths(c("~/R/library",.libPaths()))
suppressPackageStartupMessages({library(vegan); library(LDM)})
set.seed(2026); P<-"/home/ecastron/araucaria"
# ITS TABLE SELECTION. "merged" = dada2_step1.R merged-pair table (the default). "r1" =
# dada2_its_r1.R forward-read table, the pipeline chosen for the ITS analyses (R1 keeps 84.7%
# of reads against 52.1% for merged pairs). Choose with the first command-line argument or the ITS_TABLE environment
# variable. OUT_DIR redirects output so an R1 run never overwrites the merged results.
ITS_TABLE<-{a<-commandArgs(trailingOnly=TRUE); if(length(a)) a[1] else Sys.getenv("ITS_TABLE","merged")}
stopifnot(ITS_TABLE %in% c("merged","r1"))
SFX<-if(ITS_TABLE=="r1") "_r1" else ""
cat(sprintf("ITS table: %s (seqtab_nochim%s.rds, taxonomy%s.rds)\n",ITS_TABLE,SFX,SFX))
OUT<-Sys.getenv("OUT_DIR",file.path(P,"results/fungfig")); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
sink(file.path(OUT,"ldm_taxa.txt"),split=TRUE)
cat(sprintf("ITS table: %s\n",ITS_TABLE))
st<-readRDS(file.path(P,sprintf("microbiome/ITS/seqtab_nochim%s.rds",SFX)))
tax<-as.data.frame(readRDS(file.path(P,sprintf("microbiome/ITS/taxonomy%s.rds",SFX))),stringsAsFactors=FALSE)
man<-read.csv(file.path(P,"metadata/manifest microbiome runs.csv"),stringsAsFactors=FALSE,check.names=FALSE)
m<-man[match(rownames(st),man$run_accession),]
keep<-!is.na(m$tree_id)&m$compartment%in%c("needle","branch")&rowSums(st)>=1000
st<-st[keep,,drop=FALSE]; m<-m[keep,]
m$range<-factor(ifelse(m$site_code%in%c("PN","TR","VA"),"Coastal","Andean"))
st<-st[,colSums(st>0)>=3,drop=FALSE]
cnt<-rowsum(st,m$tree_id)
mt<-m[!duplicated(m$tree_id),]; mt<-droplevels(mt[match(rownames(cnt),mt$tree_id),])
r<-ldm(formula=cnt~range,data=mt,dist.method="bray",n.perm.max=10000,seed=2026)
cat("=== what LDM's VE field actually is ===\n")
cat("names(r) containing VE:\n"); print(grep("VE",names(r),value=TRUE))
for (nm in grep("VE",names(r),value=TRUE)) cat(sprintf("  %-32s %s\n",nm,paste(format(r[[nm]],digits=4),collapse=", ")))
rel<-cnt/rowSums(cnt)
a<-adonis2(vegdist(rel,method="bray")~range,data=mt,permutations=999)
cat(sprintf("\n  PERMANOVA R2 on the same table = %.4f\n",a$R2[1]))
cat("\n=== ASVs with omnibus q <= 0.10 ===\n")
q<-r$q.otu.omni
if(is.matrix(q)) q<-q[1,]
sig<-which(q<=0.10)
cat(sprintf("n = %d\n",length(sig)))
if(length(sig)){
  ids<-colnames(cnt)[sig]
  co<-colMeans(rel[mt$range=="Coastal",ids,drop=FALSE])
  an<-colMeans(rel[mt$range=="Andean", ids,drop=FALSE])
  tt<-tax[ids,c("Class","Order","Family","Genus")]
  out<-data.frame(q=signif(q[sig],3),
                  Coastal_pct=round(100*co,3),Andean_pct=round(100*an,3),
                  higher=ifelse(co>an,"Coastal","Andean"),
                  Order=sub("^o__","",tt$Order),Family=sub("^f__","",tt$Family),
                  Genus=sub("^g__","",tt$Genus),row.names=NULL)
  out<-out[order(out$q),]
  print(out,row.names=FALSE)
  write.csv(out,file.path(OUT,"ldm_taxa.csv"),row.names=FALSE)
}
cat("\nDONE\n"); sink()
