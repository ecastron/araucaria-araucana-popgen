#!/usr/bin/env Rscript
# conStruct cross-validation at 16 replicates to settle K = 2 vs K = 3.
# The 4-replicate run was ambiguous: the mean favoured spatial K=2 by 0.139 against
# sds > 0.5, while all four replicates individually preferred K=3.
suppressPackageStartupMessages(library(conStruct))
P <- "/home/ecastron/araucaria"; D <- file.path(P,"results/Q1/construct")
setwd(D); set.seed(2026)
fr <- as.matrix(read.csv("freqs_5k.csv", row.names=1, check.names=FALSE))
co <- read.csv("coords.csv", row.names=1); coords <- as.matrix(co[,c("lon","lat")])
hav <- function(m){ R<-6371; n<-nrow(m); d<-matrix(0,n,n)
  for(i in 1:n) for(j in 1:n){ p1<-m[i,2]*pi/180; p2<-m[j,2]*pi/180
    dp<-(m[j,2]-m[i,2])*pi/180; dl<-(m[j,1]-m[i,1])*pi/180
    a<-sin(dp/2)^2+cos(p1)*cos(p2)*sin(dl/2)^2; d[i,j]<-2*R*asin(sqrt(a)) }; d }
gd <- hav(coords)
cat(sprintf("conStruct CV: %d stands x %d loci, 16 replicates, K = 1:4\n", nrow(fr), ncol(fr)))
cv <- x.validation(train.prop=0.9, n.reps=16, K=1:4, freqs=fr, geoDist=gd, coords=coords,
                   prefix="xval16", n.iter=2000, make.figs=FALSE, save.files=FALSE,
                   parallel=TRUE, n.nodes=16)
saveRDS(cv, file.path(D,"construct_xval16.rds"))
Ks <- as.character(1:4); reps <- names(cv)
sp  <- sapply(reps, function(r) sapply(Ks, function(k) cv[[r]]$sp[[k]]))
nsp <- sapply(reps, function(r) sapply(Ks, function(k) cv[[r]]$nsp[[k]]))
cat(sprintf("\n%-16s %10s %10s %10s %10s\n","model","K=1","K=2","K=3","K=4"))
cat(sprintf("%-16s %10.3f %10.3f %10.3f %10.3f\n","SPATIAL mean",   rowMeans(sp)[1],rowMeans(sp)[2],rowMeans(sp)[3],rowMeans(sp)[4]))
cat(sprintf("%-16s %10.3f %10.3f %10.3f %10.3f\n","        sd",     apply(sp,1,sd)[1],apply(sp,1,sd)[2],apply(sp,1,sd)[3],apply(sp,1,sd)[4]))
cat(sprintf("%-16s %10.3f %10.3f %10.3f %10.3f\n","NONSPATIAL mean",rowMeans(nsp)[1],rowMeans(nsp)[2],rowMeans(nsp)[3],rowMeans(nsp)[4]))
cat(sprintf("%-16s %10.3f %10.3f %10.3f %10.3f\n","        sd",     apply(nsp,1,sd)[1],apply(nsp,1,sd)[2],apply(nsp,1,sd)[3],apply(nsp,1,sd)[4]))
win <- sapply(reps, function(r){ v<-c(sp1=cv[[r]]$sp[["1"]],sp2=cv[[r]]$sp[["2"]],sp3=cv[[r]]$sp[["3"]],sp4=cv[[r]]$sp[["4"]],
                                     ns1=cv[[r]]$nsp[["1"]],ns2=cv[[r]]$nsp[["2"]],ns3=cv[[r]]$nsp[["3"]],ns4=cv[[r]]$nsp[["4"]]); names(which.max(v)) })
cat("\nper-replicate winner:\n"); print(table(win))
# paired comparison spatial K2 vs K3 across replicates
d23 <- sp[2,]-sp[3,]
cat(sprintf("\nspatial K2 - K3, paired across %d reps: mean %+.3f, sd %.3f, K2 wins in %d/%d\n",
            length(d23), mean(d23), sd(d23), sum(d23>0), length(d23)))
cat(sprintf("paired t-test p = %.4f\n", tryCatch(t.test(d23)$p.value, error=function(e) NA)))
# layer contributions on the best spatial models
for (K in 2:3) {
  r <- tryCatch(conStruct(spatial=TRUE, K=K, freqs=fr, geoDist=gd, coords=coords,
                          prefix=sprintf("lc_K%d",K), n.chains=1, n.iter=2000,
                          make.figs=FALSE, save.files=FALSE), error=function(e) NULL)
  if (is.null(r)) { cat(sprintf("K=%d layer fit failed\n",K)); next }
  lc <- tryCatch(calculate.layer.contribution(conStruct.results=r[[1]], data.block=r$data.block),
                 error=function(e) tryCatch(calculate.layer.contribution(r, r$data.block),
                                            error=function(e2) paste("ERR:", conditionMessage(e2))))
  cat(sprintf("K=%d layer contributions: %s\n", K, paste(round(as.numeric(lc),4), collapse=", ")))
}
cat("\nDONE\n")
