#!/usr/bin/env Rscript
# Does sequencing depth track WITHIN-STAND distance between trees?
#
# This is the premise the Q2 robustness argument rests on. Kinship estimates are biased by
# heterozygote under-calling at low depth (median inferred F = 0.70, impossible in a dioecious
# species). That bias is harmless to a SLOPE only if it is unrelated to distance:
#
#   depth unrelated to distance -> bias adds NOISE, which cannot manufacture a slope.
#                                  Q2's Sp stands, and Q6 is attenuated but conservative.
#   depth tracks distance       -> bias tracks the predictor, so it CAN manufacture a slope,
#                                  and both Q2 and Q6 need rethinking.
#
# Depth was already shown not to be confounded with STAND (F = 1.202, p = 0.288). This asks the
# finer question: within a stand, are trees that sit close together also sequenced at similar
# depth? Plausible if field collection order tracked spatial order and libraries were prepared
# in that order.
#
# Three tests, because the bias could enter through either member of a pair or through their
# similarity:
#   1. mean pair depth        vs distance
#   2. minimum pair depth     vs distance   (the weaker library limits a pairwise estimate)
#   3. |depth difference|     vs distance   (similarity in depth, not level)

set.seed(1)
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/Q2")
sink(file.path(OUT, "depth_distance_check.txt"), split = TRUE)

# per-library call rate from the same VCF the kinship used
vcf <- file.path(P, "rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")
con <- file(vcf, "r"); samples <- NULL; called <- NULL; n <- 0
repeat {
  ln <- readLines(con, n = 20000)
  if (!length(ln)) break
  ln <- ln[!startsWith(ln, "##")]
  if (!length(ln)) next
  h <- startsWith(ln, "#CHROM")
  if (any(h)) {
    samples <- strsplit(ln[which(h)[1]], "\t")[[1]][-(1:9)]
    samples <- sub("\\.[12]$", "", samples)
    called <- numeric(length(samples))
    ln <- ln[!h]
  }
  if (!length(ln)) next
  m <- do.call(rbind, strsplit(ln, "\t", fixed = TRUE))
  g <- m[, -(1:9), drop = FALSE]
  gt <- sub(":.*", "", g)
  called <- called + colSums(gt != "./." & gt != ".")
  n <- n + nrow(m)
}
close(con)
cr <- setNames(called / n, samples)
cat(sprintf("call rate from %d sites, %d libraries: median %.3f (range %.3f-%.3f)\n\n",
            n, length(cr), median(cr), min(cr), max(cr)))

x <- readRDS(file.path(OUT, "q2_sgs.rds"))
rel <- read.delim(file.path(OUT, "ngsrelate.tsv"), check.names = FALSE)
radman <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
pop <- read.delim(file.path(P, "metadata/popmap_177.tsv"), header = FALSE,
                  col.names = c("lib", "stand"))
ch <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)

bare <- function(v) sub("\\.[12]$", "", v)
tree_of <- setNames(radman$tree_id, radman$run_accession)
stand_of <- setNames(pop$stand, pop$lib)
lat <- setNames(ch$lat, ch$tree_id); lon <- setNames(ch$lon, ch$tree_id)

A <- bare(as.character(rel$ida)); B <- bare(as.character(rel$idb))
ta <- tree_of[A]; tb <- tree_of[B]; sa <- stand_of[A]; sb <- stand_of[B]
ok <- !is.na(ta) & !is.na(tb) & ta != tb & !is.na(sa) & sa == sb &
      !is.na(lat[ta]) & !is.na(lat[tb]) & rel$nSites >= 2000
hav <- function(la1, lo1, la2, lo2) {
  R <- 6371000; p1 <- la1*pi/180; p2 <- la2*pi/180
  dp <- p2-p1; dl <- (lo2-lo1)*pi/180
  2*R*asin(sqrt(sin(dp/2)^2 + cos(p1)*cos(p2)*sin(dl/2)^2))
}
dd <- hav(lat[ta[ok]], lon[ta[ok]], lat[tb[ok]], lon[tb[ok]])
ca <- cr[A[ok]]; cb <- cr[B[ok]]; stv <- sa[ok]
keep <- dd > 0 & is.finite(dd) & !is.na(ca) & !is.na(cb)
dd <- dd[keep]; ca <- ca[keep]; cb <- cb[keep]; stv <- stv[keep]
cat(sprintf("within-stand pairs analysed: %d across %d stands\n\n", length(dd),
            length(unique(stv))))

ld <- log(dd)
tests <- list("mean pair call rate" = (ca + cb) / 2,
              "minimum pair call rate" = pmin(ca, cb),
              "|call rate difference|" = abs(ca - cb))
cat("--- correlation with ln(distance), stand held as a stratum ---\n")
for (nm in names(tests)) {
  y <- tests[[nm]]
  fit <- lm(y ~ ld + factor(stv))
  co <- summary(fit)$coefficients["ld", ]
  r <- cor(ld, y)
  cat(sprintf("  %-24s slope %+.5f  p = %.4f   (raw r = %+.3f)\n",
              nm, co[1], co[4], r))
}

cat("\n--- verdict ---\n")
ps <- sapply(tests, function(y) summary(lm(y ~ ld + factor(stv)))$coefficients["ld", 4])
if (all(ps > 0.05)) {
  cat("Depth does NOT track within-stand distance on any measure. The kinship bias therefore\n")
  cat("adds noise rather than a trend: it cannot manufacture the Q2 slope, and it attenuates\n")
  cat("Q6 toward zero rather than inflating it. Q2 stands; Q6 is conservative.\n")
} else {
  cat("Depth DOES track within-stand distance on at least one measure. The kinship bias can\n")
  cat("therefore align with the predictor, so the Q2 slope and any Q6 result must be\n")
  cat("re-examined with depth conditioned out before either is reported.\n")
}
sink()
