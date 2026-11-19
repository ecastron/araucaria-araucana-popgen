#!/usr/bin/env Rscript
# HOW MUCH F SURVIVES AT HIGH LIBRARY QUALITY?
#
# F correlates with every quality metric (|r| = 0.23-0.89 after conditioning on stand), so a
# substantial part of it is technical. That does NOT establish that all of it is. Two ways to
# ask what remains:
#   (a) STRATIFY  -- F within quality quartiles. If it flattens at the top, the residual is the
#                    biological part; if it keeps falling, the artefact is not saturated.
#   (b) EXTRAPOLATE -- regress F on call rate with stand held, and read off the intercept at a
#                    perfect call rate of 1.0. This is an extrapolation beyond the observed range
#                    (max 0.692) and is therefore indicative, NOT an estimate; it is reported with
#                    that stated rather than quoted as a corrected F.
# A dioecious species cannot have F much above 0. Whatever the residual is, it also has to be
# reconciled with the within-stand F_IS of +0.100 already measured from called genotypes.

P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q1/F_extrapolate.txt"), split = TRUE)
d <- readRDS(file.path(P, "results/Q1/F_diagnosis.rds"))
cat(sprintf("libraries: %d | call rate %.3f-%.3f | F %.3f-%.3f\n\n",
            nrow(d), min(d$callrate), max(d$callrate), min(d$F), max(d$F)))

cat("--- (a) F by call-rate quartile ---\n")
q <- cut(d$callrate, quantile(d$callrate, 0:4/4), include.lowest = TRUE,
         labels = c("Q1 worst","Q2","Q3","Q4 best"))
cat(sprintf("  %-10s %5s %14s %14s %12s\n","quartile","n","call rate","median F","mean error"))
for (l in levels(q)) {
  s <- d[q == l, ]
  cat(sprintf("  %-10s %5d %14.3f %14.4f %12.5f\n", l, nrow(s),
              median(s$callrate), median(s$F), median(s$error_est)))
}

cat("\n--- (b) extrapolation to a perfect library ---\n")
fit <- lm(F ~ callrate + factor(stand), data = d)
b <- coef(fit)
pred1 <- mean(predict(fit, transform(d, callrate = 1)))
pred_obsmax <- mean(predict(fit, transform(d, callrate = max(d$callrate))))
cat(sprintf("  slope on call rate: %+.3f per unit (p = %.2g)\n",
            b["callrate"], summary(fit)$coefficients["callrate",4]))
cat(sprintf("  predicted mean F at observed best call rate (%.3f): %.3f\n",
            max(d$callrate), pred_obsmax))
cat(sprintf("  predicted mean F extrapolated to call rate 1.000 : %.3f   <- EXTRAPOLATION\n", pred1))
cat("  (observed call rate never exceeds 0.692, so the last line is indicative only)\n")

cat("\n--- how much of the variance in F is quality vs stand? ---\n")
m0 <- lm(F ~ factor(stand), data = d)
m1 <- lm(F ~ callrate + error_est + hetero_est + factor(stand), data = d)
cat(sprintf("  stand alone            : R2 = %.3f\n", summary(m0)$r.squared))
cat(sprintf("  stand + quality metrics: R2 = %.3f\n", summary(m1)$r.squared))
cat(sprintf("  variance attributable to quality, beyond stand: %.1f%%\n",
            100*(summary(m1)$r.squared - summary(m0)$r.squared)))

cat("\n--- reconciliation with the called-genotype estimates ---\n")
cat("  ngsRelate F (genotype likelihoods, this file) : median 0.665\n")
cat("  pooled F_IS   (called genotypes, vcftools)    : +0.575\n")
cat("  within-stand F_IS (called genotypes)          : +0.100\n")
cat("  ipyrad step-4 H (raw clusters, pre-consensus) :  0.0121\n")
cat("  Expected for a dioecious outcrossing conifer  :  ~0\n")
sink()
