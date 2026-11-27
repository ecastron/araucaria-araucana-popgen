#!/usr/bin/env Rscript
# Revision (26 Nov 2026): multiple regression on distance matrices with the published
# ecodist::MRM (Lichstein 2007, Plant Ecology 188:117-131), used in the paper.
#   FST/(1-FST) ~ ln(km) + between-range indicator; 9,999 matrix permutations (seed 1),
#   full data and without VA.
# Input : figdata/host_fst_pairs.tsv (repository)
# Output: figdata/revision/ibd_mrm.json
# Run from the repository root: Rscript analysis/ibd_mrm.R
suppressPackageStartupMessages({library(ecodist); library(jsonlite)})
d0 <- read.delim("figdata/host_fst_pairs.tsv")
stopifnot(nrow(d0) == 66); d0$between <- as.logical(d0$between)
run <- function(d) {
  set.seed(1)
  st <- sort(unique(c(d$stand_a, d$stand_b))); n <- length(st)
  M <- function(v) { m <- matrix(0, n, n, dimnames = list(st, st))
    m[cbind(d$stand_a, d$stand_b)] <- v; m[cbind(d$stand_b, d$stand_a)] <- v; as.dist(m) }
  Y <- M(d$fst / (1 - d$fst)); X <- M(log(d$km)); B <- M(as.integer(d$between))
  f <- MRM(Y ~ X + B, nperm = 9999)
  list(coef_lnkm = f$coef["X", 1], p_lnkm = f$coef["X", 2],
       coef_between = f$coef["B", 1], p_between = f$coef["B", 2],
       r2 = unname(f$r.squared[1]), p_r2 = unname(f$r.squared[2]))
}
out <- list(all = run(d0), without_VA = run(d0[d0$stand_a != "VA" & d0$stand_b != "VA", ]))
str(out)
dir.create("figdata/revision", showWarnings = FALSE, recursive = TRUE)
write_json(out, "figdata/revision/ibd_mrm.json", auto_unbox = TRUE, digits = 6, pretty = TRUE)
