#!/usr/bin/env Rscript
# polyRAD step 2: genotype with population-structure + LD priors, then F by DEPTH.
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
suppressPackageStartupMessages({library(polyRAD)})
OUT <- "/home/ecastron/araucaria/results/polyrad"
rad <- readRDS(file.path(OUT,"rad_raw.rds"))
od  <- readRDS(file.path(OUT,"overdispersion_test.rds"))

# TestOverdispersion returns $optimal directly
bestd <- NA
try({ if (!is.null(od$optimal)) bestd <- as.numeric(od$optimal) }, silent=TRUE)
if (!is.finite(bestd)) bestd <- 9
cat("### overdispersion d chosen:", bestd, "\n")

cat("### IteratePopStructLD (LDdist = 150 bp = one RAD tag)\n")
set.seed(42)
rad2 <- IteratePopStructLD(rad, LDdist = 150, overdispersion = bestd)
saveRDS(rad2, file.path(OUT,"rad_genotyped.rds"))

gp <- GetProbableGenotypes(rad2, omit1allelePerLocus = TRUE)
G  <- gp$genotypes                      # taxa x alleles, copy number 0/1/2
cat("genotype matrix:", nrow(G), "taxa x", ncol(G), "alleles\n")

# per-genotype read depth = total depth at that locus for that taxon
a2l <- rad2$alleles2loc
kept <- colnames(G)
idx  <- match(kept, colnames(rad2$alleleDepth))
locOf <- a2l[idx]
D <- matrix(0, nrow(G), ncol(G), dimnames=dimnames(G))
for (L in unique(locOf)) {
  cols <- which(a2l == L)
  tot  <- rowSums(rad2$alleleDepth[, cols, drop=FALSE])
  D[, which(locOf == L)] <- tot
}

pop <- read.table("/home/ecastron/araucaria/metadata/popmap_177.tsv",
                  stringsAsFactors=FALSE)
st  <- setNames(pop[[2]], pop[[1]])
taxa <- sub("\\.[12]$","",rownames(G))
stv  <- st[taxa]

EDGES <- c(6,10,15,20,30,50,100,Inf)
NAMES <- c("6-10","10-15","15-20","20-30","30-50","50-100","100+")
res <- data.frame()
for (b in seq_len(length(EDGES)-1)) {
  ho <- 0; he <- 0; n <- 0
  for (s in unique(stv[!is.na(stv)])) {
    rows <- which(stv == s)
    if (length(rows) < 6) next
    for (j in seq_len(ncol(G))) {
      d <- D[rows, j]; g <- G[rows, j]
      keep <- which(d >= EDGES[b] & d < EDGES[b+1] & !is.na(g))
      if (length(keep) < 6) next
      gg <- g[keep]; m <- length(gg)
      q <- sum(gg)/(2*m)
      if (q <= 0 || q >= 1) next
      ho <- ho + sum(gg == 1); he <- he + 2*q*(1-q)*m; n <- n + m
    }
  }
  if (he > 0) res <- rbind(res, data.frame(depth=NAMES[b], n_geno=n,
                    Ho=ho/n, He=he/n, F_IS=1-ho/he))
}
cat("\n########## polyRAD posterior genotypes: F by DEPTH ##########\n")
print(res, row.names=FALSE)
write.csv(res, file.path(OUT,"polyrad_F_by_depth.csv"), row.names=FALSE)
cat("\nCRITERION: a correct model gives F FLAT across depth.\n")
cat("spread =", if (nrow(res)) round(max(res$F_IS)-min(res$F_IS),3) else NA, "\n")
cat("(binomial/ipyrad baseline spread was 1.06-1.13)\n")
