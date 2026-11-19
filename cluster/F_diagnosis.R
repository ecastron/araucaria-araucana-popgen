#!/usr/bin/env Rscript
# WHERE DOES F COME FROM? Two questions, neither requiring new data.
#
# (1) RECONCILE. Several F-like quantities exist and differ wildly: ngsRelate per-individual F
#     from genotype likelihoods (~0.70), pooled F_IS (+0.575), within-stand F_IS (+0.100), and
#     ipyrad's step-4 heterozygosity estimated from RAW CLUSTERS before any consensus call
#     (H = 0.0121). They are not the same estimator and the gaps between them localise the
#     problem: a normal raw H with a high called F places the artefact DOWNSTREAM of the reads.
#
# (2) DOES F TRACK LIBRARY QUALITY? A real inbreeding coefficient is a property of the tree. A
#     technical one is a property of the library. If F correlates with per-library sequencing
#     metrics after stand is held constant, that residual correlation IS residual artefact --
#     there is no biological reason a better-sequenced tree should be less inbred.
#     Stand is conditioned out because stands genuinely differ in F (Wahlund), and ignoring it
#     would let real biology masquerade as a quality effect.

P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q1/F_diagnosis.txt"), split = TRUE)

# ---- per-individual F from ngsRelate (each library appears in many pairs) -------------------
rel <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
ok <- rel$nSites >= 2000
bare0 <- function(v) sub("\\.[12]$", "", as.character(v))
fa <- data.frame(lib = bare0(c(rel$ida[ok], rel$idb[ok])),
                 F  = c(rel$Fa[ok], rel$Fb[ok]))
fa <- fa[is.finite(fa$F), ]
Fi <- aggregate(F ~ lib, fa, median)
cat(sprintf("ngsRelate per-individual F: %d libraries\n", nrow(Fi)))
cat(sprintf("  median %.4f | mean %.4f | IQR %.4f-%.4f | range %.4f-%.4f\n\n",
    median(Fi$F), mean(Fi$F), quantile(Fi$F,.25), quantile(Fi$F,.75), min(Fi$F), max(Fi$F)))

# ---- ipyrad step 4: heterozygosity and error from raw clusters ------------------------------
s4 <- read.table(file.path(P, "rad/ipyrad_full/full_ct091_clust_0.91/s4_joint_estimate.txt"),
                 header = TRUE)
bare <- function(v) sub("\\.[12]$", "", as.character(v))
s4$lib <- bare(rownames(s4))
cat(sprintf("ipyrad step-4 (raw clusters, pre-consensus): mean H = %.5f, mean error = %.5f\n",
            mean(s4$hetero_est), mean(s4$error_est)))
cat(sprintf("  H and error correlate: r = %+.3f  (joint estimation can trade them off)\n\n",
            cor(s4$hetero_est, s4$error_est)))

# ---- per-library call rate from the working VCF ---------------------------------------------
vcf <- file.path(P, "rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")
con <- file(vcf, "r"); samples <- NULL; called <- NULL; n <- 0
repeat {
  ln <- readLines(con, n = 20000); if (!length(ln)) break
  ln <- ln[!startsWith(ln, "##")]; if (!length(ln)) next
  h <- startsWith(ln, "#CHROM")
  if (any(h)) { samples <- strsplit(ln[which(h)[1]], "\t")[[1]][-(1:9)]
                called <- numeric(length(samples)); ln <- ln[!h] }
  if (!length(ln)) next
  m <- do.call(rbind, strsplit(ln, "\t", fixed = TRUE))
  gt <- sub(":.*", "", m[, -(1:9), drop = FALSE])
  called <- called + colSums(gt != "./." & gt != ".")
  n <- n + nrow(m)
}
close(con)
cr <- data.frame(lib = bare(samples), callrate = called / n, stringsAsFactors = FALSE)
cat(sprintf("call rate from %d sites: median %.3f (range %.3f-%.3f)\n\n",
            n, median(cr$callrate), min(cr$callrate), max(cr$callrate)))

# ---- read counts ----------------------------------------------------------------------------
sp <- tryCatch(read.delim(file.path(P, "metadata/library_spots.tsv")), error = function(e) NULL)
pop <- read.delim(file.path(P, "metadata/popmap_177.tsv"), header = FALSE,
                  col.names = c("lib", "stand"))

d <- merge(Fi, s4[, c("lib","hetero_est","error_est")], by = "lib")
d <- merge(d, cr, by = "lib")
d <- merge(d, pop, by = "lib")
if (!is.null(sp)) { names(sp)[1:2] <- c("lib","spots")
  sp$lib <- bare(sp$lib); d <- merge(d, sp[,c("lib","spots")], by = "lib", all.x = TRUE) }
if (nrow(d) < 20) stop(sprintf("merge collapsed to %d rows -- key mismatch, not a result", nrow(d)))
cat(sprintf("libraries with F + quality metrics + stand: %d across %d stands\n\n", nrow(d), length(unique(d$stand))))

cat("=========== (2) DOES F TRACK LIBRARY QUALITY? ===========\n")
cat("Raw correlation, then the same slope with STAND held as a stratum.\n\n")
cat(sprintf("  %-22s %10s %12s %12s\n", "metric", "raw r", "slope|stand", "p|stand"))
metrics <- list("ipyrad H (raw clusters)" = d$hetero_est,
                "ipyrad error estimate"   = d$error_est,
                "VCF call rate"           = d$callrate)
if (!is.null(sp) && "spots" %in% names(d)) metrics[["log10 read count"]] <- log10(d$spots)
for (nm in names(metrics)) {
  x <- metrics[[nm]]; keep <- is.finite(x) & is.finite(d$F)
  if (sum(keep) < 20) next
  r <- cor(x[keep], d$F[keep])
  fit <- lm(d$F[keep] ~ x[keep] + factor(d$stand[keep]))
  co <- summary(fit)$coefficients[2, ]
  cat(sprintf("  %-22s %+10.3f %+12.4f %12.2g %s\n", nm, r, co[1], co[4],
              if (co[4] < 0.05) "<< residual artefact" else ""))
}
cat("\nA significant slope AFTER conditioning on stand means F is partly a property of the\n")
cat("library rather than of the tree. There is no biological mechanism by which a\n")
cat("better-sequenced tree is less inbred.\n")

saveRDS(d, file.path(P, "results/Q1/F_diagnosis.rds"))
sink()
