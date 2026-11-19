#!/usr/bin/env Rscript
# DOES THE MICROBIOME SHOW ISOLATION BY DISTANCE AND A COASTAL-ANDEAN SPLIT?
#
# Q5 established that stand beats climate roughly 2:1. "Stand" is a container, not a mechanism:
# it bundles climate, soil, history and the local propagule pool. This asks whether the residual
# is SPATIALLY structured, which is what dispersal limitation predicts and what environmental
# filtering on unmeasured variables does not require.
#
# Built deliberately as the analogue of Q4, which found Coastal stands MORE differentiated
# genetically than distance alone predicts (+0.044, p = 0.0042). If the microbiome shows the same
# excess across the same barrier, host genes and microbial propagules are limited by the same
# landscape. If it shows plain distance decay with no cordillera term, dispersal limitation acts
# but the Central Valley is not a special barrier for microbes.
#
# THREE TESTS:
#   1. VARIATION PARTITIONING  space (PCNM) vs climate, on tree-level composition. Answers
#      whether space explains anything AFTER climate is removed -- the question Q5's marginal
#      models could not address.
#   2. DISTANCE DECAY          stand-centroid dissimilarity ~ ln(geographic distance). The slope
#      is the decay rate; significance by permuting stand labels.
#   3. CORDILLERA EXCESS       the Q4 test. dissimilarity ~ ln(distance) + cross-cordillera.
#      A positive cross-cordillera coefficient means the split adds separation BEYOND distance.
#
# POWER WARNING, STATED BEFORE THE RESULT: there are only a handful of Coastal stands. The
# cordillera term is therefore estimated from few independent units and is confounded with stand
# identity. A null on term 3 is uninformative; only a positive result carries weight, and even
# then it rests on a small number of stands. This is declared now so it cannot be waved away
# afterwards.

suppressPackageStartupMessages({library(vegan)})
set.seed(2026)
P <- "/home/ecastron/araucaria"
sink(file.path(P, "results/Q5/microbiome_ibd.txt"), split = TRUE)

man <- read.csv(file.path(P,"metadata/manifest microbiome runs.csv"), stringsAsFactors=FALSE)
ch  <- read.delim(file.path(P,"metadata/chelsa_per_tree.tsv"), check.names=FALSE)
clr <- function(m){m<-as.matrix(m); if(any(m==0)) m[m==0]<-0.65*min(m[m>0])
  t(apply(m,1,function(x) log(x)-mean(log(x))))}
hav <- function(la1,lo1,la2,lo2){R<-6371000;p1<-la1*pi/180;p2<-la2*pi/180
  2*R*asin(sqrt(sin((p2-p1)/2)^2+cos(p1)*cos(p2)*sin((lo2-lo1)*pi/180/2)^2))}
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
  b <- pool(amp); sel <- rowSums(b$cnt) >= 250
  cnt <- b$cnt[sel,,drop=FALSE]; tr <- b$trees[sel]
  cnt <- cnt[, colSums(cnt>0)>=5 & colSums(cnt)/sum(cnt)>=1e-4, drop=FALSE]
  mg <- man[match(tr, man$tree_id),]; j <- match(tr, ch$tree_id)
  lat <- ch$lat[j]; lon <- ch$lon[j]
  stand <- droplevels(factor(mg$site_code))
  # cordillera from longitude: the Central Valley is an empty corridor at ~72-73 W
  cord <- ifelse(lon < -72.5, "Costa", "Andes")
  cat(sprintf("\n################ %s: %d trees, %d ASVs, %d stands ################\n",
              amp, nrow(cnt), ncol(cnt), nlevels(stand)))
  print(table(stand, cord))
  ns <- tapply(cord, stand, function(x) x[1])
  cat(sprintf("stands by cordillera: Andes %d, Costa %d\n", sum(ns=="Andes"), sum(ns=="Costa")))

  X <- clr(cnt); D <- dist(X)
  bio <- ch[j, grep("^bio[0-9]+$", names(ch)), drop=FALSE]
  bio <- bio[, apply(bio,2,function(x) sd(x,na.rm=TRUE)>0), drop=FALSE]
  ok <- complete.cases(bio)

  # ---- 1. variation partitioning: space vs climate ----
  cat("\n--- 1. variation partitioning: space (PCNM) vs climate ---\n")
  # NOTE ON BIAS, declared before the numbers: Gilbert & Bennett (2010) showed by simulation that
  # PCNM/MEM eigenvector approaches OVER-INFLATE the space fraction while canonical ordination
  # UNDER-FITS environment. So space > climate here is the direction the known bias produces and
  # cannot support a dispersal-limitation claim; climate >= space runs AGAINST the bias and can.
  mem <- pcnm(dist(cbind(lon[ok], lat[ok])))$vectors
  mem <- mem[, seq_len(min(8, ncol(mem))), drop=FALSE]
  # 3 climate PCs, matching Q5 -- the 19 raw bioclim variables are heavily collinear (rank 17)
  cpc <- prcomp(scale(bio[ok,,drop=FALSE]))$x[,1:3]
  vp <- try(varpart(dist(X[ok,,drop=FALSE]), mem, cpc), silent=TRUE)
  if (!inherits(vp,"try-error")) {
    ind <- vp$part$indfract
    cn <- grep("^Adj", colnames(ind), value=TRUE)[1]
    f <- ind[[cn]]
    cat(sprintf("  space alone (a)   : %+.4f\n  shared    (b)    : %+.4f\n  climate alone (c): %+.4f\n  residual  (d)    : %+.4f\n",
                f[1], f[2], f[3], f[4]))
    if (any(f[1:3] < 0)) cat("  NOTE: a negative fraction is uninterpretable (predictor suppression), not a real quantity\n")
    cat(sprintf("  -> %s\n", if (f[3] >= f[1])
        "climate >= space: runs AGAINST the MEM bias, so this ordering is trustworthy" else
        "space > climate: this is the direction MEM inflation produces; NOT evidence of dispersal limitation"))
    sp <- dbrda(dist(X[ok,,drop=FALSE]) ~ mem + Condition(cpc))
    a1 <- anova(sp, permutations=9999)
    cat(sprintf("  SPACE after removing climate: F = %.3f  p = %.4f  %s\n", a1$F[1], a1$`Pr(>F)`[1],
        if (a1$`Pr(>F)`[1] < .05) "<< spatial structure beyond climate" else "no spatial structure beyond climate"))
    cl <- dbrda(dist(X[ok,,drop=FALSE]) ~ cpc + Condition(mem))
    a2 <- anova(cl, permutations=9999)
    cat(sprintf("  CLIMATE after removing space: F = %.3f  p = %.4f  %s\n", a2$F[1], a2$`Pr(>F)`[1],
        if (a2$`Pr(>F)`[1] < .05) "<< climate signal beyond space" else "no climate signal beyond space"))
  } else cat("  varpart failed:", conditionMessage(attr(vp,"condition")), "\n")

  # ---- 2/3. stand-centroid distance decay and cordillera excess ----
  cent <- t(sapply(levels(stand), function(s) colMeans(X[stand==s,,drop=FALSE])))
  slat <- tapply(lat, stand, mean); slon <- tapply(lon, stand, mean)
  S <- nlevels(stand); iu <- which(upper.tri(matrix(0,S,S)), arr.ind=TRUE)
  cd <- as.matrix(dist(cent))[iu]
  gd <- hav(slat[iu[,1]], slon[iu[,1]], slat[iu[,2]], slon[iu[,2]])
  cr <- ns[iu[,1]] != ns[iu[,2]]
  cat(sprintf("\n--- 2. distance decay across %d stand pairs (%.0f-%.0f km) ---\n",
              length(cd), min(gd)/1000, max(gd)/1000))
  f2 <- lm(cd ~ log(gd)); s2 <- summary(f2)$coefficients
  # permutation: shuffle stand labels on the centroid matrix
  B <- 9999; obs <- s2["log(gd)","t value"]; hit <- 0
  Mc <- as.matrix(dist(cent))
  for (z in 1:B) { o <- sample(S); cdp <- Mc[o,o][iu]
    hit <- hit + (abs(summary(lm(cdp ~ log(gd)))$coefficients["log(gd)","t value"]) >= abs(obs)) }
  cat(sprintf("  slope = %+.4f per ln(m)  r = %+.3f  permutation p = %.4f  %s\n",
      s2["log(gd)","Estimate"], cor(log(gd), cd), (1+hit)/(1+B),
      if ((1+hit)/(1+B) < .05) "<< ISOLATION BY DISTANCE" else "no distance decay"))

  if (sum(ns=="Costa") >= 2 && sum(cr) >= 3) {
    cat("\n--- 3. Coastal-Andean excess beyond distance (the Q4 analogue) ---\n")
    f3 <- lm(cd ~ log(gd) + cr); s3 <- summary(f3)$coefficients
    hit <- 0
    for (z in 1:B) { np <- sample(ns); crp <- np[iu[,1]] != np[iu[,2]]
      if (length(unique(crp)) < 2) next
      tt <- try(summary(lm(cd ~ log(gd) + crp))$coefficients["crpTRUE","t value"], silent=TRUE)
      if (!inherits(tt,"try-error") && abs(tt) >= abs(s3["crTRUE","t value"])) hit <- hit + 1 }
    cat(sprintf("  cross-cordillera excess = %+.4f  t = %+.3f  permutation p = %.4f  %s\n",
        s3["crTRUE","Estimate"], s3["crTRUE","t value"], (1+hit)/(1+B),
        if ((1+hit)/(1+B) < .05) "<< barrier beyond distance" else "no excess beyond distance"))
    a4 <- adonis2(D ~ cord, permutations=9999)
    cat(sprintf("  [tree-level] cordillera R2 = %.4f  p = %.4f (unconstrained; stands not held)\n",
        a4$R2[1], a4$`Pr(>F)`[1]))
  } else cat("\n--- 3. too few Coastal stands for a cordillera term; NOT FITTED ---\n")
}
cat("\nQ4 for comparison: Coastal stands +0.044 more differentiated GENETICALLY than distance predicts (p = 0.0042).\n")
sink()
