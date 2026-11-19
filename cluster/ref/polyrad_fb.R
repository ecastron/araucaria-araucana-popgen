#!/usr/bin/env Rscript
# polyRAD on the FreeBayes/reference calls -- independent replicate of the de novo run.
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
suppressPackageStartupMessages({library(polyRAD); library(VariantAnnotation)})
P <- "/home/ecastron/araucaria/results/polyrad_fb"
bg <- file.path(P,"fb_AD.vcf.gz")

# same site filters as the de novo run so the two are comparable; no per-genotype
# depth threshold, so shallow genotypes stay in.
rad <- VCF2RADdata(bg, possiblePloidies=list(2), tagsize=150,
                   min.ind.with.reads=20, min.ind.with.minor.allele=5,
                   phaseSNPs=FALSE)
cat("\n--- RADdata (FreeBayes) ---\n"); print(rad)
saveRDS(rad, file.path(P,"rad_fb_raw.rds"))
set.seed(42)
od <- TestOverdispersion(rad, to_test=c(2,3,5,7,9,12,15,20))
cat("\noverdispersion optimal:", od$optimal, "  (de novo run gave 5)\n")
saveRDS(od, file.path(P,"od_fb.rds"))
rad2 <- IteratePopStructLD(rad, LDdist=150, overdispersion=as.numeric(od$optimal))
saveRDS(rad2, file.path(P,"rad_fb_genotyped.rds"))
cat("\ngenotyping complete\n")
