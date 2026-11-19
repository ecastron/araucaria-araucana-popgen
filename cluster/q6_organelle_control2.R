#!/usr/bin/env Rscript
# Corrected organellar control. The first attempt rebuilt the ASV set from raw counts and so
# did NOT reproduce the published Q6 number (got p = 0.089 against a published p = 0.007).
# The published CLR was built on the FULL compartment library set, so its prevalence filter and
# its geometric means differ from anything rebuilt on the 49 genotyped trees. This version
# starts from the published matrix itself and only PARTITIONS ITS COLUMNS, so the baseline is
# guaranteed to match before anything is concluded from the split.

suppressPackageStartupMessages({library(vegan)})
set.seed(1)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q6/organelle_control2.txt"), split = TRUE)

st <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org_seq <- rep(FALSE, ncol(st))
for (rk in colnames(tx)) org_seq <- org_seq | grepl("Chloroplast|Mitochondri", tx[[rk]], ignore.case = TRUE)
org_seq[is.na(org_seq)] <- FALSE
names(org_seq) <- colnames(st)

rel <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
radman <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
tree_of <- setNames(radman$tree_id, radman$run_accession)
ta <- tree_of[bare(as.character(rel$ida))]; tb <- tree_of[bare(as.character(rel$idb))]
good <- !is.na(ta) & !is.na(tb) & ta != tb & rel$nSites >= 2000 & is.finite(rel$theta)
trees <- sort(unique(c(ta[good], tb[good])))
K <- matrix(NA_real_, length(trees), length(trees), dimnames = list(trees, trees))
for (i in which(good)) { K[ta[i], tb[i]] <- rel$theta[i]; K[tb[i], ta[i]] <- rel$theta[i] }
diag(K) <- max(rel$theta[good], na.rm = TRUE)

fit_one <- function(X, m, tag) {
  tr <- m$tree_id
  Ks <- K[tr, tr]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm = TRUE)
  diag(Ks) <- max(Ks, na.rm = TRUE)
  pc <- cmdscale(as.dist(max(Ks) - Ks), k = min(4, length(tr) - 2))
  colnames(pc) <- paste0("G", seq_len(ncol(pc)))
  mem <- tryCatch({ p <- pcnm(dist(cbind(m$longitude, m$latitude)))
                    p$vectors[, seq_len(min(4, ncol(p$vectors))), drop = FALSE] },
                  error = function(e) NULL)
  dat <- data.frame(pc, stand = factor(m$site_code), age = factor(m$age_class),
                    logdepth = log(rowSums(exp(X))))
  s_ <- character(0)
  if (!is.null(mem)) { colnames(mem) <- paste0("S", seq_len(ncol(mem))); dat <- cbind(dat, mem); s_ <- colnames(mem) }
  D <- dist(X)
  fo <- as.formula(paste("D ~", paste(colnames(pc), collapse = " + "),
        "+ Condition(stand + logdepth", if (length(s_)) paste("+", paste(s_, collapse = " + ")) else "", "+ age)"))
  cap <- try(dbrda(fo, data = dat), silent = TRUE)
  if (inherits(cap, "try-error") || is.null(cap$CCA)) { cat(sprintf("  %-30s FAILED\n", tag)); return(invisible()) }
  a <- anova(cap, permutations = 999)
  cat(sprintf("  %-30s ASVs=%4d  partial R2 = %.4f  F = %.3f  p = %.4f %s\n",
              tag, ncol(X), cap$CCA$tot.chi / cap$tot.chi, a$F[1], a$`Pr(>F)`[1],
              if (isTRUE(a$`Pr(>F)`[1] < 0.05)) "  << DETECTED" else ""))
}

for (comp in c("needle", "branch")) {
  o <- readRDS(file.path(P, "results/Q5", sprintf("q5_16S_%s_objects.rds", comp)))
  m <- o$manifest; X <- o$clr
  keep <- !is.na(m$tree_id) & m$tree_id %in% rownames(K)
  X <- X[keep, , drop = FALSE]; m <- m[keep, ]
  dup <- !duplicated(m$tree_id); X <- X[dup, , drop = FALSE]; m <- m[dup, ]
  cn <- colnames(X)
  o_ <- if (all(cn %in% names(org_seq))) org_seq[cn] else NULL
  cat(sprintf("\n===== 16S / %s : %d trees, %d ASVs in published matrix =====\n",
              comp, nrow(m), ncol(X)))
  if (is.null(o_)) { cat("  colnames of published CLR are not ASV sequences; cannot partition\n")
                     cat("  first colname: ", substr(cn[1], 1, 60), "\n"); next }
  cat(sprintf("  organellar columns: %d of %d\n", sum(o_), length(o_)))
  fit_one(X, m, "ALL ASVs (published baseline)")
  if (sum(o_)  >= 5) fit_one(X[,  o_, drop = FALSE], m, "ORGANELLAR columns only")
  if (sum(!o_) >= 5) fit_one(X[, !o_, drop = FALSE], m, "NON-organellar columns only")
}
sink()
