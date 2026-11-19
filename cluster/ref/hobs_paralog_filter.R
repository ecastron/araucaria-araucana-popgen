#!/usr/bin/env Rscript
# Species-wide Hobs filter for fixed ancestral duplications, on the polyRAD object.
#
# FIX 1: GetWeightedMeanGenotypes returns TAXA x LOCI (177 x 34967). rowMeans()
#        gives 177 per-tree values; the per-locus statistic needs colMeans().
# FIX 2: the scale is 0-1, not 0-2 (discrete 0/1/2 -> 0.052/0.471/0.927). The
#        het window is 0.2-0.8; "<1.8" would count hom-alt genotypes as hets.
# FIX 3: columns are ALLELES (34967) not loci (34674); aggregate to locus before
#        subsetting, or SubsetByLocus gets allele indices.
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
suppressPackageStartupMessages(library(polyRAD))
OUT <- "/home/ecastron/araucaria/results/polyrad"
r <- readRDS(file.path(OUT,"rad_genotyped.rds"))

W <- GetWeightedMeanGenotypes(r)
a2l <- r$alleles2loc[match(colnames(W), colnames(r$alleleDepth))]
het <- W > 0.2 & W < 0.8
allele_hobs <- colMeans(het, na.rm = TRUE)
locus_hobs  <- tapply(allele_hobs, a2l, max)      # worst allele defines the locus
locus_ids   <- as.integer(names(locus_hobs))
cat(sprintf("loci scored: %d\n", length(locus_hobs)))
print(round(quantile(locus_hobs, c(.05,.25,.5,.75,.9,.95,.99)),3))

G <- GetProbableGenotypes(r, omit1allelePerLocus=TRUE)$genotypes
gl <- r$alleles2loc[match(colnames(G), colnames(r$alleleDepth))]
D <- matrix(0,nrow(G),ncol(G),dimnames=dimnames(G))
for (L in unique(gl)) D[,which(gl==L)] <- rowSums(r$alleleDepth[,which(r$alleles2loc==L),drop=FALSE])
pop <- read.table("/home/ecastron/araucaria/metadata/popmap_177.tsv", stringsAsFactors=FALSE)
stv <- setNames(pop[[2]],pop[[1]])[sub("\\.[12]$","",rownames(G))]
EDG <- c(6,10,15,20,30,50,100,Inf); NM <- c("6-10","10-15","15-20","20-30","30-50","50-100","100+")

Fdepth <- function(keepLoci){
  cols <- which(gl %in% keepLoci)
  ho<-he<-ng<-numeric(length(NM))
  for (s in unique(stv[!is.na(stv)])) {
    rows <- which(stv==s); if(length(rows)<6) next
    Gs<-G[rows,cols,drop=FALSE]; Ds<-D[rows,cols,drop=FALSE]
    v <- !is.na(Gs)&Ds>=6; n_ok<-colSums(v)
    q <- ifelse(n_ok>0, colSums(ifelse(v,Gs,0),na.rm=TRUE)/(2*n_ok), NA)
    ok <- which(is.finite(q)&q>0&q<1&n_ok>=6)
    for(b in seq_along(NM)){
      inb <- Ds[,ok,drop=FALSE]>=EDG[b]&Ds[,ok,drop=FALSE]<EDG[b+1]&!is.na(Gs[,ok,drop=FALSE])
      ho[b]<-ho[b]+sum((Gs[,ok,drop=FALSE]==1)&inb,na.rm=TRUE)
      he[b]<-he[b]+sum(sweep(inb,2,2*q[ok]*(1-q[ok]),"*")); ng[b]<-ng[b]+sum(inb)
    }
  }
  list(F=1-ho/he, n=ng)
}
cat("\n%-10s %7s %s\n")
cat(sprintf("%-26s %7s  %s\n","threshold","loci",paste(sprintf("%8s",NM),collapse="")))
for (t in c(1.01, 0.45, 0.35, 0.30, 0.20)) {
  keep <- locus_ids[locus_hobs <= t]
  z <- Fdepth(keep); good <- z$n>2000
  lab <- if (t>1) "no Hobs filter" else sprintf("Hobs <= %.2f", t)
  cat(sprintf("%-26s %7d  %s   spread %.3f\n", lab, length(keep),
      paste(sprintf("%8.3f", z$F), collapse=""), max(z$F[good])-min(z$F[good])))
}
keep35 <- locus_ids[locus_hobs <= 0.35]
final <- SubsetByLocus(r, keep35)
saveRDS(final, file.path(OUT,"rad_hobs_filtered.rds"))
cat(sprintf("\nSubsetByLocus at Hobs<=0.35 -> saved rad_hobs_filtered.rds (%d loci)\n", length(keep35)))
