#!/usr/bin/env Rscript
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
suppressPackageStartupMessages(library(polyRAD))
P <- "/home/ecastron/araucaria/results/polyrad_fb"
r <- readRDS(file.path(P,"rad_fb_genotyped.rds"))
G <- GetProbableGenotypes(r, omit1allelePerLocus=TRUE)$genotypes
a2l <- r$alleles2loc; idx <- match(colnames(G), colnames(r$alleleDepth)); locOf <- a2l[idx]
D <- matrix(0,nrow(G),ncol(G),dimnames=dimnames(G))
for (L in unique(locOf)) D[,which(locOf==L)] <- rowSums(r$alleleDepth[,which(a2l==L),drop=FALSE])
pop <- read.table("/home/ecastron/araucaria/metadata/popmap_177.tsv", stringsAsFactors=FALSE)
stv <- setNames(pop[[2]],pop[[1]])[sub("\\.[12]$","",rownames(G))]
cat("taxa matched to a stand:", sum(!is.na(stv)), "of", nrow(G), "\n")
EDG <- c(6,10,15,20,30,50,100,Inf); NM <- c("6-10","10-15","15-20","20-30","30-50","50-100","100+")
ho<-he<-ng<-numeric(length(NM))
for (s in unique(stv[!is.na(stv)])) {
  rows <- which(stv==s); if(length(rows)<6) next
  Gs<-G[rows,,drop=FALSE]; Ds<-D[rows,,drop=FALSE]
  valid <- !is.na(Gs) & Ds>=6; n_ok <- colSums(valid)
  q <- ifelse(n_ok>0, colSums(ifelse(valid,Gs,0),na.rm=TRUE)/(2*n_ok), NA)
  ok <- which(is.finite(q) & q>0 & q<1 & n_ok>=6)
  for (b in seq_along(NM)) {
    inb <- Ds[,ok,drop=FALSE]>=EDG[b] & Ds[,ok,drop=FALSE]<EDG[b+1] & !is.na(Gs[,ok,drop=FALSE])
    ho[b]<-ho[b]+sum((Gs[,ok,drop=FALSE]==1)&inb,na.rm=TRUE)
    he[b]<-he[b]+sum(sweep(inb,2,2*q[ok]*(1-q[ok]),"*")); ng[b]<-ng[b]+sum(inb)
  }
}
res <- data.frame(depth=NM, n_geno=ng, Ho=ho/ng, He=he/ng, F_IS=1-ho/he)
res <- res[res$n_geno>1000,]
cat("\n##### FreeBayes + polyRAD (d=7): F by depth #####\n")
print(res, row.names=FALSE, digits=4)
write.csv(res, file.path(P,"fb_F_by_depth.csv"), row.names=FALSE)
cat("\ntotal genotypes:", sum(ng), "  spread:", round(max(res$F_IS)-min(res$F_IS),3), "\n")
