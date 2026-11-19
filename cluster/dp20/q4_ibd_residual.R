#!/usr/bin/env Rscript
# Q4: is the Coastal-Andean contrast a DIVERGENCE, or just the far end of a distance cline?
#
# The question matters because ADMIXTURE returned K = 1 -- no discrete genetic clusters -- while
# FST reaches 0.159 and the coastal stands (PN, TR, VA) are the most differentiated of the twelve.
# Both can be true if differentiation accumulates smoothly with geographic distance, in which case
# there is no divergence event to reconstruct and Q4's premise changes.
#
# THE TEST. Fit the standard isolation-by-distance regression across ALL stand pairs:
#
#     FST / (1 - FST)  ~  ln(geographic distance)          (Rousset 1997)
#
# then ask whether Coastal-Andean pairs sit ABOVE the fitted line. Under a pure cline every pair,
# whatever its cordillera, falls on the same relationship and the residuals are centred on zero.
# A genuine barrier or historical split makes between-cordillera pairs MORE differentiated than
# their distance alone predicts, i.e. systematically positive residuals.
#
#   between-cordillera residuals > 0, significantly -> real divergence beyond distance;
#                                                      Q4's colonisation question stands
#   residuals indistinguishable from zero           -> a distance cline; the "Coastal-Andean
#                                                      divergence" is geography, and Q4 should be
#                                                      reframed rather than answered
#
# Permutation null: shuffle the cordillera labels across stands, preserving the geography.

set.seed(1)
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/dp20/Q4")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "q4_summary.txt"), split = TRUE)

fst <- read.delim(file.path(P, "results/dp20/Q1/q1_fst_vcftools.tsv"), header = FALSE,
                  col.names = c("a", "b", "fst"))
rad <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
ch  <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)

cord <- tapply(rad$cordillera, rad$site_code, function(x) names(sort(table(x), TRUE))[1])
cent <- aggregate(cbind(lat, lon) ~ site_code, data = ch, FUN = mean)
rownames(cent) <- cent$site_code
cat("stands by cordillera:\n")
print(table(cord))

hav <- function(la1, lo1, la2, lo2) {
  R <- 6371; p1 <- la1*pi/180; p2 <- la2*pi/180
  dp <- p2-p1; dl <- (lo2-lo1)*pi/180
  2*R*asin(sqrt(sin(dp/2)^2 + cos(p1)*cos(p2)*sin(dl/2)^2))
}
fst$km <- mapply(function(a, b) hav(cent[a,"lat"], cent[a,"lon"], cent[b,"lat"], cent[b,"lon"]),
                 fst$a, fst$b)
fst$ca <- cord[fst$a]; fst$cb <- cord[fst$b]
fst$type <- ifelse(fst$ca == fst$cb,
                   ifelse(fst$ca == "Costa", "within-Costa", "within-Andes"), "BETWEEN")
fst <- fst[is.finite(fst$fst) & fst$km > 0, ]
cat(sprintf("\n%d stand pairs, distance %.1f-%.1f km\n", nrow(fst), min(fst$km), max(fst$km)))
print(table(fst$type))

fst$y <- fst$fst / (1 - fst$fst)
m <- lm(y ~ log(km), data = fst)
cat("\n--- isolation-by-distance regression (all pairs) ---\n")
print(summary(m)$coefficients)
cat(sprintf("R2 = %.3f\n", summary(m)$r.squared))

fst$resid <- resid(m)
cat("\n--- residuals by pair type (positive = more differentiated than distance predicts) ---\n")
agg <- aggregate(resid ~ type, data = fst, FUN = function(x) c(n = length(x), mean = mean(x)))
print(data.frame(type = agg$type, n = agg$resid[, "n"], mean_resid = round(agg$resid[, "mean"], 4)))

obs <- mean(fst$resid[fst$type == "BETWEEN"]) - mean(fst$resid[fst$type != "BETWEEN"])
cat(sprintf("\nbetween-cordillera minus within-cordillera mean residual = %+.4f\n", obs))

stands <- names(cord); nperm <- 9999; cnt <- 0
for (i in seq_len(nperm)) {
  pc <- setNames(sample(cord), stands)
  tp <- ifelse(pc[fst$a] == pc[fst$b], "within", "BETWEEN")
  d <- mean(fst$resid[tp == "BETWEEN"]) - mean(fst$resid[tp != "BETWEEN"])
  if (!is.na(d) && d >= obs) cnt <- cnt + 1
}
pval <- (cnt + 1) / (nperm + 1)
cat(sprintf("permutation p (cordillera labels shuffled across stands) = %.4f\n", pval))

cat("\nVERDICT: ")
if (pval < 0.05 && obs > 0) {
  cat("Coastal-Andean pairs ARE more differentiated than distance alone predicts.\n")
  cat("A real divergence sits on top of the cline, so Q4's colonisation question stands.\n")
} else {
  cat("Coastal-Andean pairs are NOT more differentiated than distance predicts.\n")
  cat("The contrast is a distance cline, consistent with ADMIXTURE K = 1. There is no discrete\n")
  cat("divergence event to reconstruct, and Q4 should be REFRAMED as a question about\n")
  cat("continuous structure rather than answered as a colonisation history.\n")
}

write.csv(fst[, c("a","b","fst","km","type","resid")],
          file.path(OUT, "q4_fst_distance.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %s\n", file.path(OUT, "q4_fst_distance.csv")))
sink()
