#!/usr/bin/env Rscript
# polyRAD genotype calling for ARAUCARIA (diploid, de novo RAD, 12 stands)
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
suppressPackageStartupMessages({library(polyRAD); library(VariantAnnotation)})
OUT <- "/home/ecastron/araucaria/results/polyrad"
VCF <- file.path(OUT, "ct091_cr30_AD.vcf")

bg <- paste0(VCF, ".bgz")
if (!file.exists(bg)) { bg <- bgzip(VCF, bg); indexTabix(bg, format="vcf") }
cat("bgzipped:", bg, "\n")

# min.ind.with.reads / min.ind.with.minor.allele are SITE filters, not depth filters.
# There is deliberately NO per-genotype depth threshold: low-depth genotypes are kept
# so the beta-binomial + priors can use them (min.depth.threshold = 0 in spirit).
rad <- VCF2RADdata(bg,
                   possiblePloidies = list(2),      # A. araucana is diploid
                   tagsize          = 150,          # RAD tags run to 131 bp
                   min.ind.with.reads        = 20,
                   min.ind.with.minor.allele = 5,
                   phaseSNPs = FALSE)
cat("\nRADdata object:\n"); print(rad)
saveRDS(rad, file.path(OUT, "rad_raw.rds"))

# fit the overdispersion parameter d from OUR data instead of trusting the default 9
cat("\n===== TestOverdispersion =====\n")
set.seed(42)
od <- TestOverdispersion(rad, to_test = c(2,3,5,7,9,12,15,20))
print(od)
saveRDS(od, file.path(OUT, "overdispersion_test.rds"))
