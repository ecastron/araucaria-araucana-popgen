#!/usr/bin/env Rscript
# FEASIBILITY: can the host organellar reads in the 16S metabarcoding be used as haplotype markers?
#
# 91.6% of our 16S reads are Araucaria chloroplast (75.2%) and mitochondria (16.4%). Those are
# host sequence, uniparentally inherited, and therefore a potential INDEPENDENT marker system --
# different inheritance from the nuclear RAD data, so cpDNA and mtDNA can in principle separate
# pollen flow from seed flow.
#
# THE QUESTION THAT DECIDES VIABILITY, asked before any structure is looked at: is the region
# VARIABLE within the species at all? 16S rRNA is the most conserved gene available. If cp16S V4
# is invariant in A. araucana then the 161 chloroplast ASVs are sequencing error, not haplotypes,
# and the idea dies here. The diagnostic is decisive because organelles are HAPLOID and
# uniparental: a tree should carry ONE dominant haplotype. If instead every tree shows a smear of
# ASVs at similar abundance, that is error; if trees carry one clear dominant each and those
# dominants differ between trees, that is real haplotype variation.
#
# Also separated here: SILVA files plant mitochondria under Rickettsiales, so some ASVs called
# "Mitochondria" may be genuine Rickettsia-like BACTERIA rather than host organelle.

P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q4/organelle_feasibility.txt"), split = TRUE)
dir.create(file.path(P, "results/Q4"), showWarnings = FALSE, recursive = TRUE)

st <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
man <- read.csv(file.path(P, "metadata/manifest microbiome runs.csv"), stringsAsFactors = FALSE)
ch  <- read.delim(file.path(P, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)

hit <- function(pat) { m <- apply(tx, 1, function(r) any(grepl(pat, r, ignore.case = TRUE)))
                       m[is.na(m)] <- FALSE; m }
cp <- hit("Chloroplast"); mt <- hit("Mitochondri")
cat(sprintf("chloroplast ASVs: %d (%.1f%% of reads)\n", sum(cp), 100*sum(st[,cp])/sum(st)))
cat(sprintf("mitochondrial ASVs: %d (%.1f%% of reads)\n\n", sum(mt), 100*sum(st[,mt])/sum(st)))

# are the "Mitochondria" ASVs really host, or Rickettsiales bacteria? host mtDNA should be
# hugely abundant and present in nearly every library; a real bacterium should not be.
prevm <- colSums(st[, mt, drop=FALSE] > 0); relm <- colSums(st[, mt, drop=FALSE])/sum(st[, mt, drop=FALSE])
cat(sprintf("mito ASVs in >90%% of libraries: %d | carrying %.1f%% of mito reads\n",
            sum(prevm > 0.9*nrow(st)), 100*sum(relm[prevm > 0.9*nrow(st)])))

i <- match(rownames(st), man$run_accession)
tid <- man$tree_id[i]; comp <- man$compartment[i]
j <- match(tid, ch$tree_id)
cord <- ifelse(ch$lon[j] < -72.5, "Costa", "Andes")

report <- function(sel, lab) {
  cat(sprintf("\n########## %s ##########\n", lab))
  sub <- st[, sel, drop = FALSE]
  keep <- rowSums(sub) >= 500 & !is.na(tid)
  sub <- sub[keep, , drop = FALSE]
  tt <- tid[keep]; cc <- cord[keep]; kk <- comp[keep]
  cat(sprintf("libraries with >=500 %s reads: %d across %d trees\n", lab, nrow(sub), length(unique(tt))))
  if (nrow(sub) < 20) { cat("  too few; skipped\n"); return(invisible()) }

  # pool per tree (a tree's organellar haplotype must be identical in every tissue)
  trees <- sort(unique(tt))
  pool <- t(sapply(trees, function(x) colSums(sub[tt == x, , drop = FALSE])))
  frac <- pool / rowSums(pool)
  dom  <- colnames(pool)[apply(frac, 1, which.max)]
  domf <- apply(frac, 1, max)

  cat(sprintf("\n  dominance of the top ASV per tree: median %.3f (range %.3f-%.3f)\n",
              median(domf), min(domf), max(domf)))
  cat("  -> a real haploid haplotype should dominate its tree; a smear means error\n")
  cat(sprintf("  trees where the top ASV holds >50%% of reads: %d of %d (%.0f%%)\n",
              sum(domf > .5), length(domf), 100*mean(domf > .5)))
  cat(sprintf("  distinct dominant ASVs across trees: %d\n", length(unique(dom))))

  if (length(unique(dom)) < 2) { cat("  INVARIANT -- one haplotype only; not usable as a marker\n"); return(invisible()) }

  crd <- cord[match(trees, tid)]
  tab <- table(dom, crd)
  keeph <- rowSums(tab) >= 3
  cat("\n  dominant haplotype x cordillera (haplotypes in >=3 trees):\n")
  print(tab[keeph, , drop = FALSE])
  if (sum(keeph) >= 2 && ncol(tab) == 2) {
    ft <- fisher.test(tab[keeph, , drop = FALSE], simulate.p.value = TRUE, B = 9999)
    cat(sprintf("\n  Fisher exact, haplotype vs cordillera: p = %.4f  %s\n", ft$p.value,
        if (ft$p.value < 0.05) "<< haplotypes partition by cordillera" else "no cordillera partition"))
  }
  # consistency check: does a tree give the same dominant haplotype in different tissues?
  multi <- names(which(table(tt) > 1))
  if (length(multi) >= 5) {
    agree <- sapply(multi, function(x) {
      r <- sub[tt == x, , drop = FALSE]
      length(unique(colnames(r)[apply(r, 1, which.max)])) == 1 })
    cat(sprintf("\n  same dominant ASV across tissues of the same tree: %d of %d (%.0f%%)\n",
                sum(agree), length(agree), 100*mean(agree)))
    cat("  -> low agreement means the signal is noise, not a heritable haplotype\n")
  }
}
report(cp, "CHLOROPLAST")
report(mt, "MITOCHONDRIA")
sink()
