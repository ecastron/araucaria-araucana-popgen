#!/usr/bin/env Rscript
# AMOVA for Q1: partition genetic variance among stands versus among trees within stands.
#
# adonis2 on the genetic distance matrix is a distance-based AMOVA and needs only vegan.
#
# NA HANDLING. The first attempt required libraries with ZERO undefined pairs and only 1 of 177
# qualified, so it failed outright. Undefined pairs arise where two libraries share fewer than
# 200 called sites -- 4.2% of pairs overall, but concentrated in a few poor libraries: 170 of 177
# have under 10% missing. Rather than impute (which invents data at exactly the pairs we know
# least about), greedily drop the library with the most undefined pairs until the matrix is
# complete, and report how many were removed.

suppressPackageStartupMessages(library(vegan))
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/Q1")
d <- as.matrix(read.delim(file.path(OUT, "q1_genetic_distance.tsv"), row.names = 1,
                          check.names = FALSE))
g <- read.delim(file.path(OUT, "q1_stands.tsv"))
stopifnot(nrow(d) == nrow(g))
cat(sprintf("input: %d libraries, %.2f%% undefined pairs\n", nrow(d), 100 * mean(is.na(d))))

keep <- rep(TRUE, nrow(d))
while (any(is.na(d[keep, keep]))) {
  sub <- d[keep, keep]
  worst <- which.max(rowSums(is.na(sub)))
  keep[which(keep)[worst]] <- FALSE
}
cat(sprintf("dropped %d libraries to obtain a complete matrix; %d retained\n",
            sum(!keep), sum(keep)))
d <- d[keep, keep]; g <- g[keep, ]
tb <- table(g$stand)
cat("stands retained: ", paste(names(tb), tb, sep = "=", collapse = "  "), "\n\n")
if (length(tb) < 2) stop("fewer than 2 stands survive; AMOVA not possible")

a <- adonis2(as.dist(d) ~ stand, data = g, permutations = 9999)
print(a)
r2 <- a$R2[1]
cat(sprintf("\namong stands  : %.2f%% of genetic variance (p = %.4f)\n", 100 * r2,
            a$`Pr(>F)`[1]))
cat(sprintf("within stands : %.2f%%\n", 100 * (1 - r2)))
cat("\nThe spec's expectation -- to be tested, not assumed -- is that MOST variance sits WITHIN\n")
cat("stands. If it does not, the within-stand design behind Q6 is not identifiable, and that\n")
cat("must be stated before proceeding rather than discovered later.\n")

# Per-stand beta dispersion: adonis2 assumes comparable within-group spread, so an unequal
# dispersion would make a significant stand effect ambiguous between location and spread.
bd <- betadisper(as.dist(d), factor(g$stand))
pt <- permutest(bd, permutations = 999)
cat(sprintf("\nbeta dispersion homogeneity: p = %.4f%s\n", pt$tab$`Pr(>F)`[1],
            if (pt$tab$`Pr(>F)`[1] < 0.05)
              "  <-- UNEQUAL; the stand effect is partly dispersion, not only location" else
              "  (equal spread; the stand effect is a location difference)"))

capture.output(a, pt, file = file.path(OUT, "q1_amova.txt"))
cat("\nwritten:", file.path(OUT, "q1_amova.txt"), "\n")
