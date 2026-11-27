.libPaths(c("~/R/library",.libPaths()))
suppressPackageStartupMessages({library(vegan)})
set.seed(2026); P<-"/home/ecastron/araucaria"
# ITS TABLE SELECTION. "merged" = dada2_step1.R merged-pair table (the default). "r1" =
# dada2_its_r1.R forward-read table, the pipeline chosen for the ITS analyses (R1 keeps 84.7%
# of reads against 52.1% for merged pairs). Choose with the first command-line argument or the ITS_TABLE environment
# variable. OUT_DIR redirects output so an R1 run never overwrites the merged results.
ITS_TABLE<-{a<-commandArgs(trailingOnly=TRUE); if(length(a)) a[1] else Sys.getenv("ITS_TABLE","merged")}
stopifnot(ITS_TABLE %in% c("merged","r1"))
SFX<-if(ITS_TABLE=="r1") "_r1" else ""
cat(sprintf("ITS table: %s (seqtab_nochim%s.rds, taxonomy%s.rds)\n",ITS_TABLE,SFX,SFX))
st<-readRDS(file.path(P,sprintf("microbiome/ITS/seqtab_nochim%s.rds",SFX)))
man<-read.csv(file.path(P,"metadata/manifest microbiome runs.csv"),stringsAsFactors=FALSE,check.names=FALSE)
m<-man[match(rownames(st),man$run_accession),]
keep<-!is.na(m$tree_id)&m$compartment%in%c("needle","branch")&rowSums(st)>=1000
st<-st[keep,,drop=FALSE]; m<-m[keep,]
m$range<-ifelse(m$site_code%in%c("PN","TR","VA"),"Coastal","Andean")
REL<-st/rowSums(st)
# tree level, the agreed unit
nt<-rowsum(REL,m$tree_id)/as.vector(table(m$tree_id)); nt<-nt/rowSums(nt)
mt<-m[!duplicated(m$tree_id),]; mt<-mt[match(rownames(nt),mt$tree_id),]
BCt<-vegdist(nt,method="bray")
bd<-betadisper(BCt,factor(mt$range))
cat("TREE LEVEL, Bray-Curtis\n")
cat(sprintf("  n: Coastal %d, Andean %d  (unbalanced %.1f:1)\n",
    sum(mt$range=="Coastal"),sum(mt$range=="Andean"),
    sum(mt$range=="Andean")/sum(mt$range=="Coastal")))
mn<-tapply(bd$distances,bd$group,mean)
cat(sprintf("  mean distance-to-centroid: Coastal %.4f | Andean %.4f\n",mn["Coastal"],mn["Andean"]))
pt<-permutest(bd,permutations=9999)$tab
cat(sprintf("  PERMDISP F = %.3f  p = %.4f\n",pt[1,4],pt[1,6]))
smaller<-"Coastal"
cat(sprintf("\n  the SMALLER group (%s, n=%d) is %s dispersed than the larger\n",
    smaller,sum(mt$range==smaller),
    if(mn["Coastal"]>mn["Andean"]) "MORE" else "LESS"))
cat(sprintf("  -> under Anderson & Walsh (2013) an unbalanced design with the smaller group\n"))
cat(sprintf("     %s dispersed makes PERMANOVA %s, so the range p is %s\n",
    if(mn["Coastal"]>mn["Andean"]) "MORE" else "LESS",
    if(mn["Coastal"]>mn["Andean"]) "CONSERVATIVE" else "LIBERAL",
    if(mn["Coastal"]>mn["Andean"]) "if anything understated" else "INFLATED - treat with caution"))
a<-adonis2(BCt~range,data=mt,permutations=9999)
cat(sprintf("\n  PERMANOVA range: R2 = %.4f  p = %.4f\n",a$R2[1],a$`Pr(>F)`[1]))
# a dispersion-insensitive cross-check: ANOSIM works on ranks, and a balanced subsample
cat("\n  cross-checks\n")
an<-anosim(BCt,mt$range,permutations=9999)
cat(sprintf("  ANOSIM R = %.4f p = %.4f\n",an$statistic,an$signif))
ps<-replicate(200,{i<-c(which(mt$range=="Coastal"),sample(which(mt$range=="Andean"),sum(mt$range=="Coastal")))
  adonis2(as.dist(as.matrix(BCt)[i,i])~range,data=mt[i,],permutations=999)$`Pr(>F)`[1]})
cat(sprintf("  balanced subsampling 20 vs 20, 200 draws: median p = %.4f, %.0f%% below 0.05\n",
    median(ps),100*mean(ps<0.05)))
