#!/usr/bin/env Rscript
# Compares results/reconstructed/its_* (from its_figdata.R) with the 26 Aug 2026 originals in
# results/figdata/. Reports byte identity (md5) and, per numeric column, max |difference|.
# Reconstructed 18-19 Nov 2026 to replace lost inline code (verification helper).
P <- "/home/ecastron/araucaria"; A <- file.path(P, "results/figdata"); B <- file.path(P, "results/reconstructed")
for (f in c("its_prevalence.csv", "its_richness.csv", "its_composition.csv", "its_ordination.csv", "its_stand_pairs.tsv")) {
  fa <- file.path(A, f); fb <- file.path(B, f)
  same <- unname(tools::md5sum(fa) == tools::md5sum(fb))
  rd <- if (grepl("tsv$", f)) function(p) read.delim(p) else function(p) read.csv(p, check.names = FALSE)
  a <- rd(fa); b <- rd(fb)
  cat(sprintf("\n%s  byte-identical: %s  dims %s vs %s\n", f, same, paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x")))
  if (!identical(dim(a), dim(b))) next
  for (cn in names(a)) {
    if (is.numeric(a[[cn]])) {
      d <- abs(a[[cn]] - b[[cn]])
      # PCA axes are sign-indeterminate: report the better of the two orientations too
      dflip <- abs(a[[cn]] + b[[cn]])
      cat(sprintf("  %-12s max|diff| %.3g  (sign-flipped %.3g)  r = %.6f\n", cn, max(d), max(dflip), cor(a[[cn]], b[[cn]])))
    } else cat(sprintf("  %-12s identical: %s\n", cn, identical(as.character(a[[cn]]), as.character(b[[cn]]))))
  }
}
