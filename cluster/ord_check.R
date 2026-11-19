.libPaths(c("~/R/library",.libPaths()))
suppressPackageStartupMessages({library(vegan)})
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
st<-readRDS(file.path(P,sprintf("microbiome/ITS/seqtab_nochim%s.rds",SFX)))
man<-read.csv(file.path(P,"metadata/manifest microbiome runs.csv"),stringsAsFactors=FALSE,check.names=FALSE)
m<-man[match(rownames(st),man$run_accession),]
keep<-!is.na(m$tree_id) & m$compartment %in% c("needle","branch") & rowSums(st)>=1000
st<-st[keep,,drop=FALSE]; m<-m[keep,]
m$range<-ifelse(m$site_code %in% c("PN","TR","VA"),"Coastal","Andean")
BC<-vegdist(st/rowSums(st),method="bray")
cat("design balance: range x compartment\n"); print(table(m$range,m$compartment))
cat("\n--- SEQUENTIAL, range first ---\n")
print(anova(dbrda(BC~range+compartment,data=m),permutations=9999,by="terms"))
cat("\n--- SEQUENTIAL, compartment first ---\n")
print(anova(dbrda(BC~compartment+range,data=m),permutations=9999,by="terms"))
cat("\n--- MARGINAL (each term adjusted for the other) ---\n")
print(anova(dbrda(BC~range+compartment,data=m),permutations=9999,by="margin"))
cat("\n--- how much the two terms overlap ---\n")
cat(sprintf("  Cramer-style association, range vs compartment: chisq p = %.3f\n",
    suppressWarnings(chisq.test(table(m$range,m$compartment))$p.value)))
cap<-dbrda(BC~range+compartment,data=m)
cat(sprintf("  joint constrained R2 = %.4f (rank %d)\n",cap$CCA$tot.chi/cap$tot.chi,cap$CCA$rank))
r1<-dbrda(BC~range,data=m); r2<-dbrda(BC~compartment,data=m)
cat(sprintf("  range alone %.4f | compartment alone %.4f | sum %.4f vs joint %.4f -> shared %.5f\n",
    r1$CCA$tot.chi/r1$tot.chi, r2$CCA$tot.chi/r2$tot.chi,
    r1$CCA$tot.chi/r1$tot.chi + r2$CCA$tot.chi/r2$tot.chi,
    cap$CCA$tot.chi/cap$tot.chi,
    r1$CCA$tot.chi/r1$tot.chi + r2$CCA$tot.chi/r2$tot.chi - cap$CCA$tot.chi/cap$tot.chi))
