#!/usr/bin/env Rscript
# WHAT SIZE OF HOST GENETIC EFFECT COULD Q6 HAVE DETECTED?
#
# Replaces the untestable caveat "a null is ambiguous" with a number. A signal of KNOWN size,
# structured by the real kinship matrix, is injected into data whose own signal has been
# destroyed by permutation; the actual Q6 test is then run and detection recorded. The smallest
# effect detected 80% of the time is the sensitivity of the design.
#
# CONSTRUCTION, so the answer is not rigged:
#   Yr        observed CLR residualised on the Q6 covariates (stand, age, space, depth)
#   Yr[perm]  rows permuted -> any real host-genetic signal destroyed, but the ASV correlation
#             structure, sparsity and variance heterogeneity of the REAL data are preserved.
#             Simulating from a clean multivariate normal would flatter the design.
#   S = G L   signal built from the real kinship PCoA axes with random loadings, so the injected
#             effect has the same geometry as the effect being looked for.
#   Y = sqrt(1-w) Yr[perm] + sqrt(w) S,  scaled so S accounts for a fraction w of total variance.
#
# w = 0 IS INCLUDED AS A CALIBRATION CHECK. If the test is correctly sized, detection at w = 0
# must come out near 0.05. If it does not, the power curve means nothing and is not reported.

suppressPackageStartupMessages({library(vegan); library(parallel)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
CORES <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","16")))
NREP <- 200
sink(file.path(P, "results/Q6/power.txt"), split = TRUE)
cat(sprintf("replicates per effect size: %d | cores %d\n", NREP, CORES))

man <- read.csv(file.path(P,"metadata/manifest microbiome runs.csv"), stringsAsFactors=FALSE)
rad <- read.csv(file.path(P,"metadata/manifest rad runs.csv"), stringsAsFactors=FALSE)
ch  <- read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"), check.names=FALSE)
rf  <- read.delim(file.path(P,"results/Q2/ngsrelate.tsv"), check.names=FALSE)
bare <- function(v) sub("\\.[12]$","",v)
to <- setNames(rad$tree_id, rad$run_accession)
ta <- to[bare(as.character(rf$ida))]; tb <- to[bare(as.character(rf$idb))]
gd <- !is.na(ta)&!is.na(tb)&ta!=tb&rf$nSites>=2000&is.finite(rf$theta)
kt <- sort(unique(c(ta[gd],tb[gd])))
K0 <- matrix(NA_real_,length(kt),length(kt),dimnames=list(kt,kt))
for (z in which(gd)) { K0[ta[z],tb[z]] <- rf$theta[z]; K0[tb[z],ta[z]] <- rf$theta[z] }
diag(K0) <- max(rf$theta[gd],na.rm=TRUE)
clr <- function(m){m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0])
  t(apply(m,1,function(x) log(x)-mean(log(x))))}
pool <- function(amp){
  st <- readRDS(file.path(P,sprintf("microbiome/%s/seqtab_nochim.rds",amp)))
  if (amp=="16S"){tx<-as.data.frame(readRDS(file.path(P,"microbiome/16S/taxonomy.rds")),stringsAsFactors=FALSE)
    o<-rep(FALSE,ncol(st)); for(rk in colnames(tx)) o<-o|grepl("Chloroplast|Mitochondri",tx[[rk]],ignore.case=TRUE)
    o[is.na(o)]<-FALSE; st<-st[,!o,drop=FALSE]}
  i<-match(rownames(st),man$run_accession)
  k<-man$compartment[i]%in%c("branch","needle")&!is.na(man$tree_id[i])
  st<-st[k,,drop=FALSE]; tid<-man$tree_id[i][k]; tr<-sort(unique(tid))
  list(cnt=t(sapply(tr,function(x) colSums(st[tid==x,,drop=FALSE]))),trees=tr)}

for (amp in c("ITS","16S")) {
  b <- pool(amp); sel <- rowSums(b$cnt)>=250
  cnt <- b$cnt[sel,,drop=FALSE]; tr <- b$trees[sel]
  g <- tr%in%rownames(K0); cnt<-cnt[g,,drop=FALSE]; tg<-tr[g]
  cnt <- cnt[,colSums(cnt>0)>=5 & colSums(cnt)/sum(cnt)>=1e-4,drop=FALSE]
  n <- nrow(cnt); m <- ncol(cnt)
  Y <- clr(cnt)
  K <- K0[tg,tg]; K[is.na(K)]<-mean(K[upper.tri(K)],na.rm=TRUE); diag(K)<-max(K,na.rm=TRUE); K<-(K+t(K))/2
  mg <- man[match(tg,man$tree_id),]; j <- match(tg,ch$tree_id)
  G <- cmdscale(as.dist(max(K)-K),k=4); colnames(G)<-paste0("G",1:4)
  S4 <- pcnm(dist(cbind(ch$lon[j],ch$lat[j])))$vectors[,1:4,drop=FALSE]; colnames(S4)<-paste0("S",1:4)
  covs <- data.frame(G,S4,stand=droplevels(factor(mg$site_code)),age=factor(mg$age_class),
                     logdepth=log(rowSums(exp(Y))))
  Xc <- model.matrix(~ stand+age+S1+S2+S3+S4+logdepth, data=covs)
  Yr <- resid(lm.fit(Xc, Y))                      # covariate-residualised observed data
  cat(sprintf("\n############ %s: n=%d trees, %d ASVs ############\n", amp, n, m))

  test1 <- function(Ysim){
    d <- data.frame(covs, check.names=FALSE)
    cap <- try(dbrda(dist(Ysim) ~ G1+G2+G3+G4 +
                Condition(stand+logdepth+S1+S2+S3+S4+age), data=d), silent=TRUE)
    if (inherits(cap,"try-error")||is.null(cap$CCA)) return(NA_real_)
    a <- anova(cap, permutations=999); a$`Pr(>F)`[1]
  }
  cat(sprintf("  %-8s %10s %10s %s\n","effect w","power","mean p","(w = fraction of variance from host genetics)"))
  for (w in c(0, .02, .04, .06, .08, .12, .16, .25)) {
    ps <- unlist(mclapply(seq_len(NREP), function(r){
      Yp <- Yr[sample(n),,drop=FALSE]
      L  <- matrix(rnorm(ncol(G)*m), ncol(G), m)
      Sg <- G %*% L
      Sg <- Sg / sqrt(sum(Sg^2)); Yp <- Yp / sqrt(sum(Yp^2))
      test1(sqrt(1-w)*Yp + sqrt(w)*Sg)
    }, mc.cores=CORES))
    ps <- ps[is.finite(ps)]
    cat(sprintf("  %-8.2f %9.3f %10.3f %s\n", w, mean(ps<0.05), mean(ps),
        if (w==0) "  <- calibration: should be near 0.05" else ""))
  }
}
cat("\nThe smallest w with power >= 0.80 is the effect size this design could reliably detect.\n")
cat("Effects below it are NOT excluded by the observed nulls; effects above it ARE.\n")
sink()
