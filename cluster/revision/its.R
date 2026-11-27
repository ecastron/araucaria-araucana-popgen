#!/usr/bin/env Rscript
# Revision (25 Nov 2026): endophyte range test with the stand as the unit, and
# sensitivity of Shannon diversity to the rarefaction depth.
#
# Inputs (as fungal_fig.R with ITS_TABLE = r1): microbiome/ITS/seqtab_nochim_r1.rds,
#   metadata/manifest microbiome runs.csv; needle and branch libraries with a tree_id and >= 1,000 reads
#   (239 libraries, 79 trees, 9 stands).
# 1. Range effect at tree level (adonis2, Bray-Curtis on tree-mean relative abundances, as Figure 5A),
#    tested against every assignment of 2 of the 9 stands to the Coastal range (36 labelings, exact).
#    Stands are nested in range, so range cannot be conditioned on stand; permuting whole stands is
#    the test that treats the stand as the replicate. Minimum attainable p = 1/36.
# 2. Shannon diversity at rarefaction depths 1,062 (paper), 2,000, 5,000 and 10,000 reads (libraries
#    below the depth dropped; 20 rarefactions, mean Shannon per library), and without rarefaction
#    (Shannon of relative abundances). Same three tests as Figure 5C: range within needles and within
#    branches (Wilcoxon, libraries), compartment paired within tree; Benjamini-Hochberg over the three.
# Outputs: results/revision/its_stand_perm.tsv, its_rarefaction_sensitivity.tsv, printed summary
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages(library(vegan))
set.seed(2026); P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/revision")
st <- readRDS(file.path(P, "microbiome/ITS/seqtab_nochim_r1.rds"))
man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
m <- man[match(rownames(st), man$run_accession), ]
keep <- !is.na(m$tree_id) & m$compartment %in% c("needle", "branch") & rowSums(st) >= 1000
st <- st[keep, , drop = FALSE]; m <- m[keep, ]
COASTAL <- c("PN", "TR", "VA"); m$range <- ifelse(m$site_code %in% COASTAL, "Coastal", "Andean")
stopifnot(nrow(st) == 239)

# ---- 1. stand-level permutation of range ----
REL <- st / rowSums(st)
nt <- rowsum(REL, m$tree_id) / as.vector(table(m$tree_id)); nt <- nt / rowSums(nt)
mt <- m[!duplicated(m$tree_id), ]; mt <- mt[match(rownames(nt), mt$tree_id), ]
BCt <- vegdist(nt, "bray")
r2_of <- function(lab) adonis2(BCt ~ lab, permutations = 0)$R2[1]
obs <- r2_of(factor(mt$range))
stands <- sort(unique(mt$site_code)); nco <- sum(stands %in% COASTAL)
cmb <- combn(stands, nco)
null <- apply(cmb, 2, function(cc) r2_of(factor(mt$site_code %in% cc)))
p_stand <- mean(null >= obs - 1e-12)
labs <- apply(cmb, 2, paste, collapse = "+")
write.table(data.frame(coastal_stands = labs, R2 = round(null, 5), observed = labs == paste(sort(stands[stands %in% COASTAL]), collapse = "+")),
            file.path(OUT, "its_stand_perm.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("range R2 at tree level = %.4f; rank among %d stand labelings = %d; exact stand-level p = %.4f\n",
            obs, ncol(cmb), sum(null >= obs - 1e-12), p_stand))

# ---- 2. rarefaction sensitivity ----
tests <- function(sh) {
  C <- data.frame(sh = sh, range = m$range, comp = m$compartment, tree = m$tree_id)
  p <- c(sapply(c("needle", "branch"), function(cp)
           wilcox.test(C$sh[C$comp == cp & C$range == "Coastal"], C$sh[C$comp == cp & C$range == "Andean"], exact = FALSE)$p.value))
  pw <- merge(C[C$comp == "needle", c("tree", "sh")], C[C$comp == "branch", c("tree", "sh")], by = "tree")
  pw <- pw[!duplicated(pw$tree), ]
  p <- c(p, compartment = wilcox.test(pw$sh.x, pw$sh.y, paired = TRUE, exact = FALSE)$p.value)
  c(n_lib = nrow(C), setNames(p.adjust(p, "BH"), paste0("pBH_", c("range_needle", "range_branch", "compartment"))),
    med_coastal = median(C$sh[C$range == "Coastal"]), med_andean = median(C$sh[C$range == "Andean"]))
}
rows <- list()
for (dep in c(1062, 2000, 5000, 10000)) {
  ok <- rowSums(st) >= dep
  sh <- rowMeans(sapply(1:20, function(i) diversity(rrarefy(st[ok, ], dep), "shannon")))
  full <- rep(NA_real_, nrow(st)); full[ok] <- sh
  mm <- m; m <- m[ok, ]; r <- tests(sh); m <- mm
  rows[[length(rows) + 1]] <- data.frame(depth = dep, libraries_dropped = sum(!ok), t(r))
}
r <- tests(diversity(REL, "shannon"))
rows[[length(rows) + 1]] <- data.frame(depth = "none (relative abundance)", libraries_dropped = 0, t(r))
tab <- do.call(rbind, rows)
write.table(tab, file.path(OUT, "its_rarefaction_sensitivity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(tab, digits = 3)
