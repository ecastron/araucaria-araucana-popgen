#!/usr/bin/env Rscript
# Q5: how does endophyte composition partition among environment, space and host, at what scale?
# Analysis plan is Step 5 of the spec.
#
# CLR / robust Aitchison throughout, so 16S and ITS sit on identical, phylogeny-free footing.
# ITS has no tree (not alignable across fungi at the observed 76-487 bp span); UniFrac is a
# 16S-side supplement only and never carries a 16S-versus-ITS contrast.
#
# STRATIFIED BY COMPARTMENT, not pooled. The deposit holds branch (137), needle (136) and root
# (23) samples. Compartment is a first-order driver of endophyte composition, so pooling would
# put the largest source of variation into the residual and flatter everything else. Root is
# reported DESCRIPTIVELY only: all 23 root samples are seedlings, so compartment is perfectly
# confounded with age class there and no inference can separate them (spec, Section 7).
#
# THE POSITIVE CONTROL IS NOT OPTIONAL AND RUNS FIRST. Environment must dominate BETWEEN stands
# -- 98% of mean annual temperature variance is among stands, over a 577 m elevation and 2.5-fold
# precipitation span. If an environmental signal cannot be recovered where one demonstrably
# exists, a weak within-stand result means nothing. The control's verdict is printed with every
# downstream result so the two are never read apart.

suppressPackageStartupMessages({library(vegan); library(phyloseq)})
set.seed(1)

PROJ <- "/home/ecastron/araucaria"
OUT  <- file.path(PROJ, "results/Q5")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "q5_summary.txt"), split = TRUE)
cat("Q5 community analysis\n\n")

man <- read.csv(file.path(PROJ, "metadata/manifest microbiome runs.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
chelsa <- read.delim(file.path(PROJ, "metadata/chelsa_per_tree.tsv"), check.names = FALSE)
terrain <- tryCatch(read.delim(file.path(PROJ, "metadata/terrain_per_tree.tsv"),
                               check.names = FALSE), error = function(e) NULL)
cat(sprintf("manifest %d runs | chelsa %d trees | terrain %s\n",
            nrow(man), nrow(chelsa), if (is.null(terrain)) "absent" else nrow(terrain)))
cat("compartments: ", paste(names(table(man$compartment)),
                            table(man$compartment), sep = "=", collapse = "  "), "\n\n")

# CLR with multiplicative replacement: zeros are structural in compositional data, so replace
# only the zeros rather than adding a pseudocount to every cell (which distorts rare taxa).
clr <- function(m) {
  m <- as.matrix(m)
  if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0])
  t(apply(m, 1, function(x) log(x) - mean(log(x))))
}

run_one <- function(amp, comp, st_all, man_all) {
  cat(sprintf("\n---------- %s / %s ----------\n", amp, comp))
  sel <- man_all$compartment == comp
  if (sum(sel) < 20) { cat("  <20 libraries; descriptive only\n"); return(invisible(NULL)) }
  st <- st_all[sel, , drop = FALSE]; m <- man_all[sel, ]
  st <- st[, colSums(st > 0) >= 2, drop = FALSE]
  cat(sprintf("  %d libraries x %d ASVs\n", nrow(st), ncol(st)))

  env <- chelsa[match(m$tree_id, chelsa$tree_id), ]
  if (!is.null(terrain)) {
    tr <- terrain[match(m$tree_id, terrain$tree_id), ]
    add <- setdiff(names(tr), names(env))
    if (length(add)) env <- cbind(env, tr[, add, drop = FALSE])
  }
  stand <- factor(m$site_code)
  bio <- env[, grep("^bio[0-9]+$", names(env)), drop = FALSE]
  bio <- bio[, apply(bio, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  ok <- complete.cases(bio) & !is.na(m$latitude) & !is.na(m$longitude)
  cat(sprintf("  complete environment + coordinates: %d | stands: %d\n",
              sum(ok), nlevels(droplevels(stand[ok]))))
  if (sum(ok) < 20) { cat("  too few complete cases\n"); return(invisible(NULL)) }

  X <- clr(st[ok, , drop = FALSE])
  D <- dist(X)
  Dm <- as.matrix(D)
  sd_ <- droplevels(stand[ok])
  logdepth <- log(rowSums(st[ok, , drop = FALSE]))
  pc <- prcomp(scale(bio[ok, , drop = FALSE]))$x[, 1:3]

  cat("\n  [positive control] environment BETWEEN stands\n")
  a1 <- adonis2(D ~ sd_, permutations = 999)
  a2 <- adonis2(D ~ pc, permutations = 999)
  cat(sprintf("    stand identity  : R2 = %.4f  p = %.4f\n", a1$R2[1], a1$`Pr(>F)`[1]))
  cat(sprintf("    climate PC1-3   : R2 = %.4f  p = %.4f\n", a2$R2[1], a2$`Pr(>F)`[1]))
  passed <- isTRUE(a2$`Pr(>F)`[1] < 0.05)
  cat(sprintf("    CONTROL %s\n", if (passed)
      "PASSED - downstream results interpretable" else
      "FAILED - do NOT interpret anything below"))

  mem <- tryCatch({
    p <- pcnm(dist(cbind(m$longitude[ok], m$latitude[ok])))
    p$vectors[, seq_len(min(8, ncol(p$vectors))), drop = FALSE]
  }, error = function(e) { cat("    spatial eigenvectors failed:", conditionMessage(e), "\n"); NULL })

  cat("\n  [variance partitioning] climate | space | library size\n")
  if (!is.null(mem)) {
    vp <- try(varpart(D, ~ pc, ~ mem, ~ logdepth), silent = TRUE)
    if (!inherits(vp, "try-error")) {
      print(vp$part$indfract[, c("Adj.R.square")], digits = 3)
      capture.output(vp, file = file.path(OUT, sprintf("q5_%s_%s_varpart.txt", amp, comp)))
    }
  }

  cat("\n  [dbRDA] climate + space, conditioned on log library size\n")
  dat <- data.frame(pc, logdepth = logdepth)
  if (!is.null(mem)) dat <- cbind(dat, as.data.frame(mem))
  pred <- setdiff(names(dat), "logdepth")
  cap <- try(dbrda(as.formula(paste("D ~", paste(pred, collapse = " + "),
                                    "+ Condition(logdepth)")), data = dat), silent = TRUE)
  if (!inherits(cap, "try-error")) {
    print(anova(cap, permutations = 999))
    cat(sprintf("    constrained: %.1f%% of inertia\n", 100 * cap$CCA$tot.chi / cap$tot.chi))
  }

  cat("\n  [within stands] mean Aitchison distance among trees inside a stand\n")
  tb <- table(sd_)
  for (s in names(tb)[tb >= 6]) {
    i <- which(sd_ == s)
    cat(sprintf("    %-6s n=%2d  %.2f\n", s, length(i), mean(as.dist(Dm[i, i]))))
  }

  saveRDS(list(clr = X, dist = D, env = env[ok, ], stand = sd_, manifest = m[ok, ],
               control_passed = passed),
          file.path(OUT, sprintf("q5_%s_%s_objects.rds", amp, comp)))
  invisible(NULL)
}

for (a in list(c("16S", "seqtab_nochim.rds"), c("ITS", "seqtab_nochim_r1.rds"))) {
  amp <- a[1]; f <- file.path(PROJ, "microbiome", amp, a[2])
  cat(sprintf("\n================ %s ================\n", amp))
  if (!file.exists(f)) { cat("  seqtab missing:", f, "\n"); next }
  st <- readRDS(f)
  m <- man[match(rownames(st), man$run_accession), ]
  keep <- !is.na(m$tree_id) & m$tree_id %in% chelsa$tree_id & rowSums(st) >= 1000
  cat(sprintf("  %d libraries; with tree_id, coordinates and >=1000 reads: %d\n",
              nrow(st), sum(keep)))
  st <- st[keep, , drop = FALSE]; m <- m[keep, ]

  cat("\n  [descriptive] compartment effect, all compartments pooled\n")
  Xa <- clr(st[, colSums(st > 0) >= 2, drop = FALSE])
  Da <- dist(Xa)
  ad <- try(adonis2(Da ~ compartment + site_code, data = m, permutations = 999), silent = TRUE)
  if (!inherits(ad, "try-error")) print(ad)
  cat("  (root = 23 seedlings only; compartment is confounded with age class there,\n")
  cat("   so root is DESCRIPTIVE only and carries no inference -- spec Section 7)\n")

  for (comp in c("branch", "needle", "root")) {
    tryCatch(run_one(amp, comp, st, m),
             error = function(e) cat("\n  ERROR", amp, comp, ":", conditionMessage(e), "\n"))
  }
}
cat("\nDONE\n")
sink()
