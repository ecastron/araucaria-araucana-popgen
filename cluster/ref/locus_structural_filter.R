# Structural locus filtering for polyRAD RADdata  (corrected)
library(polyRAD)
my_polyrad <- readRDS("/home/ecastron/araucaria/results/polyrad/rad_genotyped.rds")

# FIX 1 -- locDepth is TAXA x LOCI (177 x 34674), not loci x taxa. Transpose,
#          otherwise rowMeans() returns 177 per-TREE values and the subsequent
#          which() indexes trees, silently keeping <=177 arbitrary loci.
locus_depths <- t(my_polyrad$locDepth)                 # loci x taxa

# FIX 2 -- a 0 means MISSING, not "sequenced at depth zero". Leaving zeros in
#          makes mean/sd/cv measure missingness: with ~63% zeros every CV came
#          out >1.3 and mean depth was deflated from 21.8 to 8.0.
present <- locus_depths; present[present == 0] <- NA
locus_mean_depth <- rowMeans(present, na.rm = TRUE)
locus_sd_depth   <- apply(present, 1, sd, na.rm = TRUE)
locus_cv_depth   <- locus_sd_depth / locus_mean_depth
locus_missingness <- rowSums(locus_depths == 0, na.rm = TRUE) / ncol(locus_depths)

# FIX 3 -- thresholds must come from THIS dataset. It is cr30, so a locus needs
#          data in only 30% of trees; missingness runs 0.475-0.689 (median 0.633)
#          and miss<=0.05 retains ZERO loci. Observed: cv p50 0.628 p95 0.787;
#          mean depth p50 21.8 p95 43.2 p99 521.6.
max_cv_threshold <- 0.80      # ~p85
max_locus_mean   <- 64.0      # ~3x the median present-genotype depth
max_missing_rate <- 0.65      # ~p60 for this cr30 matrix

good_loci_indices <- which(locus_cv_depth   <= max_cv_threshold &
                           locus_mean_depth <= max_locus_mean   &
                           locus_missingness<= max_missing_rate)
cat(sprintf("Retained %d of %d loci (%.1f%%)\n", length(good_loci_indices),
            nrow(locus_depths), 100*length(good_loci_indices)/nrow(locus_depths)))

clean_polyrad <- SubsetByLocus(my_polyrad, good_loci_indices)
weighted_genotypes <- GetWeightedMeanGenotypes(clean_polyrad)
saveRDS(clean_polyrad, "/home/ecastron/araucaria/results/polyrad/rad_locusfiltered.rds")

# NOTE (measured 27 Aug 2026): this filter does NOT flatten the F-vs-depth curve.
# spread 1.207 unfiltered -> 1.549 at 58% retention -> 0.878 at 19% retention,
# and the 6-20x end is invariant at 0.76-0.78 under every threshold tried.
