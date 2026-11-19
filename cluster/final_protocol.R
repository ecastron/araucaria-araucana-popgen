#!/usr/bin/env Rscript
# FINAL ENDOPHYTE PROTOCOL. Inherits the published MRA 2020 pipeline for this exact dataset
# (PRJNA517193) and adds the two things it lacks.
#
# FROM MRA 2020 (Alarcon et al., 10.1128/mra.00544-20), unchanged:
#   DADA2 maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE   [already matched by our dada2 run]
#   SILVA (16S) / UNITE (ITS) taxonomy
#   samples with <1000 reads not considered
# ADDED, because the announcement did not address them:
#   (1) chloroplast + mitochondria removed. In SILVA these classify INSIDE the bacterial
#       taxonomy -- mitochondria under Alphaproteobacteria/Rickettsiales, chloroplast under
#       Cyanobacteria -- so they pass every read-count filter and enter community models as
#       if they were bacteria. They are 91.63% of 16S reads here.
#   (2) a sample floor of >=500 BACTERIAL reads for 16S, and an ASV filter of prevalence >=5
#       libraries AND >=1e-4 relative abundance. The MRA floor is on TOTAL reads and is
#       non-binding on this subset: all 273 needle+branch libraries already exceed 1000.
#
# CONSEQUENCE ACCEPTED IN ADVANCE: the >=500 bacterial floor removes 16S needle entirely
# (16 libraries survive, below the 20 needed). That is the correct outcome for a compartment
# whose median bacterial content is 27 reads, and it is reported as an exclusion, not hidden.

suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q5/final_protocol.txt"), split = TRUE)

st16 <- readRDS(file.path(P, "microbiome/16S/seqtab_nochim.rds"))
tx16 <- as.data.frame(readRDS(file.path(P, "microbiome/16S/taxonomy.rds")), stringsAsFactors = FALSE)
org <- rep(FALSE, ncol(st16))
for (rk in colnames(tx16)) org <- org | grepl("Chloroplast|Mitochondri", tx16[[rk]], ignore.case = TRUE)
org[is.na(org)] <- FALSE
stI <- readRDS(file.path(P, "microbiome/ITS/seqtab_nochim.rds"))

rel <- read.delim(file.path(P, "results/Q2/ngsrelate.tsv"), check.names = FALSE)
radman <- read.csv(file.path(P, "metadata/manifest rad runs.csv"), stringsAsFactors = FALSE)
bare <- function(v) sub("\\.[12]$", "", v)
tree_of <- setNames(radman$tree_id, radman$run_accession)
ta <- tree_of[bare(as.character(rel$ida))]; tb <- tree_of[bare(as.character(rel$idb))]
gd <- !is.na(ta) & !is.na(tb) & ta != tb & rel$nSites >= 2000 & is.finite(rel$theta)
trees <- sort(unique(c(ta[gd], tb[gd])))
K <- matrix(NA_real_, length(trees), length(trees), dimnames = list(trees, trees))
for (i in which(gd)) { K[ta[i], tb[i]] <- rel$theta[i]; K[tb[i], ta[i]] <- rel$theta[i] }
diag(K) <- max(rel$theta[gd], na.rm = TRUE)

clr <- function(m) { m <- as.matrix(m)
  if (any(m == 0)) m[m == 0] <- 0.65 * min(m[m > 0]); t(apply(m, 1, function(x) log(x) - mean(log(x)))) }

run <- function(amp, comp, raw, decon, floor_) {
  f <- file.path(P, "results/Q5", sprintf("q5_%s_%s_objects.rds", amp, comp))
  if (!file.exists(f)) return(invisible())
  o <- readRDS(f); m <- o$manifest; env <- o$env
  ids <- intersect(rownames(raw), m$run_accession)
  k <- match(ids, m$run_accession); m <- m[k, ]; env <- env[k, ]
  cnt <- raw[ids, , drop = FALSE]
  if (decon) cnt <- cnt[, !org, drop = FALSE]
  sel <- rowSums(cnt) >= floor_
  cat(sprintf("\n---- %s / %s : %d -> %d libraries at >=%d %s reads ----\n", amp, comp,
              length(sel), sum(sel), floor_, if (decon) "BACTERIAL" else "total"))
  if (sum(sel) < 20) { cat("   EXCLUDED: <20 libraries survive\n"); return(invisible()) }
  cnt <- cnt[sel, , drop = FALSE]; m <- m[sel, ]; env <- env[sel, ]
  cnt <- cnt[, colSums(cnt > 0) >= 5 & colSums(cnt) / sum(cnt) >= 1e-4, drop = FALSE]
  cat(sprintf("   ASVs after prevalence>=5 & rel>=1e-4: %d\n", ncol(cnt)))
  if (ncol(cnt) < 10) { cat("   EXCLUDED: <10 ASVs\n"); return(invisible()) }
  X <- clr(cnt); D <- dist(X)
  bio <- env[, grep("^bio[0-9]+$", names(env)), drop = FALSE]
  bio <- bio[, apply(bio, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  ok <- complete.cases(bio); stand <- droplevels(factor(m$site_code))
  if (nlevels(stand) >= 2 && sum(ok) >= 20) {
    a1 <- adonis2(D ~ stand, permutations = 9999)
    a2 <- adonis2(as.dist(as.matrix(D)[ok, ok]) ~ prcomp(scale(bio[ok, ]))$x[, 1:3], permutations = 9999)
    cat(sprintf("   [Q5] stand R2 = %.4f (p=%.4f) | climate R2 = %.4f (p=%.4f) | ratio %.2f\n",
                a1$R2[1], a1$`Pr(>F)`[1], a2$R2[1], a2$`Pr(>F)`[1], a1$R2[1] / a2$R2[1]))
  }
  # Q6
  g <- !is.na(m$tree_id) & m$tree_id %in% rownames(K)
  Xg <- X[g, , drop = FALSE]; mg <- m[g, ]
  d2 <- !duplicated(mg$tree_id); Xg <- Xg[d2, , drop = FALSE]; mg <- mg[d2, ]
  if (nrow(mg) < 25 || nlevels(droplevels(factor(mg$site_code))) < 2) {
    cat(sprintf("   [Q6] %d genotyped trees -- below 25, NOT FITTED\n", nrow(mg))); return(invisible()) }
  Ks <- K[mg$tree_id, mg$tree_id]; Ks[is.na(Ks)] <- mean(Ks[upper.tri(Ks)], na.rm = TRUE)
  diag(Ks) <- max(Ks, na.rm = TRUE)
  pc <- cmdscale(as.dist(max(Ks) - Ks), k = 4); colnames(pc) <- paste0("G", 1:4)
  mem <- pcnm(dist(cbind(mg$longitude, mg$latitude)))$vectors[, 1:4, drop = FALSE]
  colnames(mem) <- paste0("S", 1:4)
  dat <- data.frame(pc, mem, stand = droplevels(factor(mg$site_code)),
                    age = factor(mg$age_class), logdepth = log(rowSums(exp(Xg))))
  cap <- try(dbrda(dist(Xg) ~ G1+G2+G3+G4 + Condition(stand + logdepth + S1+S2+S3+S4 + age),
                   data = dat), silent = TRUE)
  if (inherits(cap, "try-error") || is.null(cap$CCA)) { cat("   [Q6] model failed\n"); return(invisible()) }
  a <- anova(cap, permutations = 9999)
  cat(sprintf("   [Q6] n=%d  host genetics partial R2 = %.4f  F = %.3f  p = %.4f  %s\n",
              nrow(mg), cap$CCA$tot.chi / cap$tot.chi, a$F[1], a$`Pr(>F)`[1],
              if (a$`Pr(>F)`[1] < 0.05) "<< DETECTED" else "not detected"))
}

cat("################ 16S : organelles removed, >=500 bacterial reads ################\n")
for (cc in c("needle", "branch", "root")) run("16S", cc, st16, TRUE, 500)
cat("\n\n################ ITS : MRA floor >=1000 reads ################\n")
for (cc in c("needle", "branch", "root")) run("ITS", cc, stI, FALSE, 1000)
sink()
