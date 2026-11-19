#!/usr/bin/env Rscript
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
suppressPackageStartupMessages(library(polyRAD))
OUT <- "/home/ecastron/araucaria/results/polyrad"
rad2 <- readRDS(file.path(OUT,"rad_genotyped.rds"))
gp <- GetProbableGenotypes(rad2, omit1allelePerLocus = TRUE)
G  <- gp$genotypes
a2l <- rad2$alleles2loc
idx <- match(colnames(G), colnames(rad2$alleleDepth))
locOf <- a2l[idx]
D <- matrix(0, nrow(G), ncol(G), dimnames=dimnames(G))
for (L in unique(locOf)) {
  cols <- which(a2l == L)
  D[, which(locOf == L)] <- rowSums(rad2$alleleDepth[, cols, drop=FALSE])
}
pop <- read.table("/home/ecastron/araucaria/metadata/popmap_177.tsv", stringsAsFactors=FALSE)
st  <- setNames(pop[[2]], pop[[1]]); stv <- st[sub("\\.[12]$","",rownames(G))]

EDGES <- c(6,10,15,20,30,50,100,Inf)
NAMES <- c("6-10","10-15","15-20","20-30","30-50","50-100","100+")
ho <- he <- ng <- numeric(length(NAMES))
for (s in unique(stv[!is.na(stv)])) {
  rows <- which(stv == s); if (length(rows) < 6) next
  Gs <- G[rows,,drop=FALSE]; Ds <- D[rows,,drop=FALSE]
  # allele frequency per locus per stand from ALL depths -> He is a locus property
  cnt <- colSums(Gs, na.rm=TRUE); n_ok <- colSums(!is.na(Gs) & Ds>=6)
  q <- ifelse(n_ok>0, cnt/(2*n_ok), NA)
  ok <- which(is.finite(q) & q>0 & q<1 & n_ok>=6)
  for (b in seq_along(NAMES)) {
    inb <- Ds[,ok,drop=FALSE] >= EDGES[b] & Ds[,ok,drop=FALSE] < EDGES[b+1]
    gsub_ <- Gs[,ok,drop=FALSE]
    hets  <- (gsub_ == 1) & inb & !is.na(gsub_)
    nb    <- inb & !is.na(gsub_)
    ho[b] <- ho[b] + sum(hets)
    he[b] <- he[b] + sum(sweep(nb, 2, 2*q[ok]*(1-q[ok]), "*"))
    ng[b] <- ng[b] + sum(nb)
  }
}
res <- data.frame(depth=NAMES, n_geno=ng, Ho=ho/ng, He=he/ng, F_IS=1-ho/he)
res <- res[res$n_geno>500,]
cat("\n##### polyRAD (d=5, PopStructLD) -- F by depth, locus-level allele freq #####\n")
print(res, row.names=FALSE, digits=4)
write.csv(res, file.path(OUT,"polyrad_F_by_depth_fixed.csv"), row.names=FALSE)
cat("\ntotal genotypes used:", sum(ng), "\n")
cat("spread across depth  :", round(max(res$F_IS)-min(res$F_IS),3), "\n")
