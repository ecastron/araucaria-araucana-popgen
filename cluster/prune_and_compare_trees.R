#!/usr/bin/env Rscript
# Prune the SEPP insertion tree to our own ASVs, then ask whether the reference-guided tree and
# the de novo FastTree tree actually lead to different answers.
#
# Two separate issues, deliberately kept apart:
#
# 1. TREE SIZE. The raw SEPP output carries the whole reference backbone (208,666 tips) because
#    insertion adds our fragments to the reference tree rather than extracting a subtree. That
#    is an artefact of the output format, not a property of the phylogeny. Pruning to our 5,214
#    ASVs preserves the topology and branch lengths the reference supplied while giving a tree
#    whose tips are exactly our features. UniFrac restricts itself to the subtree spanning
#    observed features anyway, so pruning should not change any distance -- which is testable,
#    and is tested below rather than assumed.
#
# 2. WHETHER THE REFERENCE IS EARNING ITS KEEP. The argument for insertion is that 253 bp of V4
#    cannot resolve deep branching, so a de novo tree invents it. The argument against is that
#    the reference is external, dated (GG 13_8, 2013), and could misplace lineages it does not
#    contain. Rather than argue, compare: if UniFrac under the two trees gives the same
#    ordination, the choice is inconsequential and the simpler de novo tree is defensible. If
#    it does not, the difference is exactly what the sensitivity check exists to expose.
#
# Note on "5,214 / 5,214 placed": SEPP places every fragment somewhere. Full placement is NOT
# evidence of good placement, and should not be reported as if it were.

suppressPackageStartupMessages({library(ape); library(vegan)})

PROJ <- "/home/ecastron/araucaria"
P    <- file.path(PROJ, "microbiome/16S/phylo")
seqtab <- readRDS(file.path(PROJ, "microbiome/16S/seqtab_nochim_md5.rds"))
ids  <- colnames(seqtab)

sepp <- read.tree(file.path(P, "asv_16S_sepp.nwk"))
den  <- read.tree(file.path(P, "asv_16S_denovo.nwk"))
cat(sprintf("SEPP raw tips: %d | de novo tips: %d | our ASVs: %d\n",
            Ntip(sepp), Ntip(den), length(ids)))

# ---- 1. prune ---------------------------------------------------------------------------
keep <- intersect(sepp$tip.label, ids)
stopifnot(length(keep) == length(ids))
sepp_p <- keep.tip(sepp, keep)
write.tree(sepp_p, file.path(P, "asv_16S_sepp_pruned.nwk"))
cat(sprintf("pruned SEPP tree: %d tips -> %s\n", Ntip(sepp_p),
            file.path(P, "asv_16S_sepp_pruned.nwk")))

# Did pruning change the phylogenetic distances among our ASVs? It must not.
sub <- sample(keep, 400)
d_raw <- cophenetic(keep.tip(sepp, sub)); d_prn <- cophenetic(keep.tip(sepp_p, sub))
d_prn <- d_prn[rownames(d_raw), colnames(d_raw)]
cat(sprintf("pruning check: max abs difference in pairwise distance = %.3g (should be ~0)\n",
            max(abs(d_raw - d_prn))))

# ---- 2. do the two trees describe the same phylogenetic structure? -----------------------
common <- intersect(sepp_p$tip.label, den$tip.label)
s2 <- keep.tip(sepp_p, common); d2 <- keep.tip(den, common)
rf <- dist.topo(s2, d2, method = "PH85")
cat(sprintf("\nRobinson-Foulds distance: %.0f (max possible ~%d)\n", rf, 2 * (length(common) - 3)))

sub2 <- sample(common, 500)
cs <- cophenetic(keep.tip(s2, sub2)); cd <- cophenetic(keep.tip(d2, sub2))
cd <- cd[rownames(cs), colnames(cs)]
cat(sprintf("cophenetic distance correlation (Mantel r): %.3f\n",
            mantel(as.dist(cs), as.dist(cd), permutations = 199)$statistic))

# ---- 3. the decision-relevant test: does UniFrac change? ---------------------------------
# Presence/absence UniFrac is used because it is the metric most exposed to tree topology;
# abundance weighting would mask exactly the differences we are trying to detect.
if (requireNamespace("phyloseq", quietly = TRUE)) {
  suppressPackageStartupMessages(library(phyloseq))
  ss <- sample(rownames(seqtab), min(60, nrow(seqtab)))
  tb <- seqtab[ss, , drop = FALSE]; tb <- tb[, colSums(tb) > 0, drop = FALSE]
  tips <- intersect(colnames(tb), common)
  tb <- tb[, tips, drop = FALSE]

  ps_s <- phyloseq(otu_table(tb, taxa_are_rows = FALSE), phy_tree(keep.tip(s2, tips)))
  ps_d <- phyloseq(otu_table(tb, taxa_are_rows = FALSE), phy_tree(keep.tip(d2, tips)))

  for (w in c(FALSE, TRUE)) {
    u_s <- UniFrac(ps_s, weighted = w); u_d <- UniFrac(ps_d, weighted = w)
    m <- mantel(u_s, u_d, permutations = 999)
    p <- protest(cmdscale(u_s, k = 2), cmdscale(u_d, k = 2), permutations = 999)
    cat(sprintf("\n%s UniFrac, SEPP vs de novo (n = %d samples, %d taxa)\n  Mantel r = %.3f (p = %.3f) | Procrustes = %.3f (p = %.3f)\n",
                ifelse(w, "weighted", "unweighted"), nrow(tb), ncol(tb),
                m$statistic, m$signif, p$scale, p$signif))
  }
  cat("\nIf these correlations are high the tree choice does not drive the ecology, and the\nsimpler de novo tree is defensible. If they are low, the trees disagree where it matters\nand the reference-guided one is preferred on the published benchmark.\n")
} else {
  cat("\nphyloseq not installed; UniFrac comparison skipped\n")
}
