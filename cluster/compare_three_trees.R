#!/usr/bin/env Rscript
# Three 16S trees now exist over the same 5,214 ASVs:
#
#   asv_16S_sepp_gg2_pruned.nwk  SEPP into Greengenes2 2022.10 (whole-genome WoL backbone)
#   asv_16S_sepp_pruned.nwk      SEPP into Greengenes 13_8     (16S-only backbone, 2013)
#   asv_16S_denovo.nwk           MAFFT + FastTree from our reads alone
#
# GG2 is the declared primary. This quantifies what that choice is worth: if UniFrac is
# effectively identical across all three, the decision is inconsequential and should be reported
# as such rather than defended. If GG2 and GG 13_8 agree with each other but both differ from de
# novo, the reference is doing real work and the de novo tree would have misled us.
#
# Unweighted UniFrac is included deliberately alongside weighted: it is the metric most exposed
# to topology, so it is where any difference will show. Weighted UniFrac is dominated by
# abundant taxa and tends to agree regardless of tree.

suppressPackageStartupMessages({library(ape); library(vegan); library(phyloseq)})
set.seed(1)

PROJ <- "/home/ecastron/araucaria"
P    <- file.path(PROJ, "microbiome/16S/phylo")
seqtab <- readRDS(file.path(PROJ, "microbiome/16S/seqtab_nochim_md5.rds"))

trees <- list(
  gg2    = file.path(P, "asv_16S_sepp_gg2_pruned.nwk"),
  gg13_8 = file.path(P, "asv_16S_sepp_pruned.nwk"),
  denovo = file.path(P, "asv_16S_denovo.nwk")
)
trees <- lapply(trees[file.exists(unlist(trees))], read.tree)
cat("trees loaded:", paste(names(trees), sapply(trees, Ntip), sep = "=", collapse = "  "), "\n")

common <- Reduce(intersect, lapply(trees, function(t) t$tip.label))
common <- intersect(common, colnames(seqtab))
cat(sprintf("ASVs common to all trees and the table: %d\n\n", length(common)))

# ---- topology agreement -------------------------------------------------------------------
sub <- sample(common, min(500, length(common)))
cat("pairwise cophenetic-distance correlation (Mantel r) on", length(sub), "ASVs\n")
cop <- lapply(trees, function(t) cophenetic(keep.tip(t, sub)))
nm <- names(trees)
for (i in seq_along(nm)) for (j in seq_along(nm)) if (i < j) {
  a <- cop[[i]]; b <- cop[[j]][rownames(cop[[i]]), colnames(cop[[i]])]
  cat(sprintf("  %-8s vs %-8s r = %.3f\n", nm[i], nm[j],
              mantel(as.dist(a), as.dist(b), permutations = 199)$statistic))
}

# ---- does the tree change the ecology? ----------------------------------------------------
ss <- sample(rownames(seqtab), min(60, nrow(seqtab)))
tb <- seqtab[ss, , drop = FALSE]
tb <- tb[, colnames(tb) %in% common, drop = FALSE]
tb <- tb[, colSums(tb) > 0, drop = FALSE]
tb <- tb[rowSums(tb) > 0, , drop = FALSE]
tips <- colnames(tb)
cat(sprintf("\nUniFrac on %d samples x %d ASVs\n", nrow(tb), ncol(tb)))

ps <- lapply(trees, function(t)
  phyloseq(otu_table(tb, taxa_are_rows = FALSE), phy_tree(keep.tip(t, tips))))

for (w in c(FALSE, TRUE)) {
  cat(sprintf("\n%s UniFrac\n", ifelse(w, "weighted", "unweighted")))
  u <- lapply(ps, function(x) UniFrac(x, weighted = w))
  for (i in seq_along(nm)) for (j in seq_along(nm)) if (i < j) {
    m <- mantel(u[[i]], u[[j]], permutations = 999)
    p <- protest(cmdscale(u[[i]], k = 2), cmdscale(u[[j]], k = 2), permutations = 999)
    cat(sprintf("  %-8s vs %-8s Mantel r = %.3f (p = %.3f) | Procrustes = %.3f\n",
                nm[i], nm[j], m$statistic, m$signif, p$scale))
  }
}
cat("\nHigh correlations across all three => the tree choice does not drive conclusions, and\n")
cat("that should be REPORTED rather than the primary being defended on authority.\n")
