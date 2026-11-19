#!/usr/bin/env Rscript
# conStruct: is the Coastal-Andean signal DISCRETE ancestry or a steep CLINE?
# Decision rule fixed before running:
#   spatial model wins at K=1              -> continuous IBD; no discrete layer needed
#   nonspatial K=2 beats spatial K=1       -> discrete split not reducible to distance
#   spatial K=2 wins                       -> a real layer sitting ON TOP of the cline
if (!requireNamespace("conStruct", quietly=TRUE)) { cat("conStruct NOT INSTALLED - skipping\n"); quit(status=0) }
suppressPackageStartupMessages(library(conStruct))
P <- "/home/ecastron/araucaria"; D <- file.path(P,"results/Q1/construct")
setwd(D); set.seed(2026)
fr <- as.matrix(read.csv("freqs_5k.csv", row.names=1, check.names=FALSE))
co <- read.csv("coords.csv", row.names=1)
coords <- as.matrix(co[,c("lon","lat")])
stopifnot(nrow(fr)==nrow(coords))
hav <- function(m){ R<-6371; n<-nrow(m); d<-matrix(0,n,n)
  for(i in 1:n) for(j in 1:n){ p1<-m[i,2]*pi/180; p2<-m[j,2]*pi/180
    dp<-(m[j,2]-m[i,2])*pi/180; dl<-(m[j,1]-m[i,1])*pi/180
    a<-sin(dp/2)^2+cos(p1)*cos(p2)*sin(dl/2)^2; d[i,j]<-2*R*asin(sqrt(a)) }
  d }
gd <- hav(coords)
cat(sprintf("conStruct input: %d stands x %d loci | geo dist %.1f-%.1f km\n",
            nrow(fr), ncol(fr), min(gd[upper.tri(gd)]), max(gd[upper.tri(gd)])))
res <- list()
for (sp in c(TRUE,FALSE)) for (K in 1:3) {
  tag <- sprintf("%s_K%d", ifelse(sp,"spatial","nonspatial"), K)
  cat(sprintf("\n--- %s ---\n", tag))
  r <- tryCatch(conStruct(spatial=sp, K=K, freqs=fr, geoDist=if(sp) gd else NULL,
                          coords=coords, prefix=tag, n.chains=1, n.iter=2000,
                          make.figs=FALSE, save.files=FALSE),
                error=function(e){ cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  lc <- tryCatch(calculate.layer.contribution(r, data.block=r$data.block),
                 error=function(e) NA)
  cat(sprintf("  layer contributions: %s\n", paste(sprintf("%.3f", lc), collapse=", ")))
  res[[tag]] <- list(lc=lc)
}
saveRDS(res, file.path(D,"construct_layers.rds"))
cat("\n--- cross-validation (spatial vs nonspatial, K=1:3) ---\n")
cv <- tryCatch(x.validation(train.prop=0.9, n.reps=4, K=1:3, freqs=fr, geoDist=gd,
                            coords=coords, prefix="xval", n.iter=1000,
                            make.figs=FALSE, save.files=FALSE, parallel=FALSE),
               error=function(e){ cat("x.validation FAILED:", conditionMessage(e), "\n"); NULL })
if (!is.null(cv)) { print(cv); saveRDS(cv, file.path(D,"construct_xval.rds")) }
cat("\nDONE\n")
