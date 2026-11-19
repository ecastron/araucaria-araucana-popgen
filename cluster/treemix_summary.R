#!/usr/bin/env Rscript
# How much of the allele-frequency covariance does each TreeMix model explain, and where do the
# migration edges go? A large gain from m=0 to m=1 that then flattens indicates one real exchange
# event rather than noise being progressively fitted.
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/Q4/treemix")
setwd(OUT)
cat(sprintf("%-6s %14s %14s\n", "model", "likelihood", "var. explained"))
res <- data.frame()
for (m in 0:5) {
  ll <- file.path(OUT, sprintf("tm_m%d.llik", m))
  cv <- file.path(OUT, sprintf("tm_m%d.cov.gz", m))
  mv <- file.path(OUT, sprintf("tm_m%d.modelcov.gz", m))
  if (!file.exists(ll)) { cat(sprintf("%-6s %14s\n", paste0("m=", m), "absent")); next }
  L <- as.numeric(sub(".*: ", "", tail(readLines(ll), 1)))
  ve <- NA
  if (file.exists(cv) && file.exists(mv)) {
    a <- as.matrix(read.table(gzfile(cv), header = TRUE, row.names = 1))
    b <- as.matrix(read.table(gzfile(mv), header = TRUE, row.names = 1))
    ve <- 1 - var(as.vector(a - b)) / var(as.vector(a))
  }
  cat(sprintf("m=%-4d %14.2f %13.4f\n", m, L, ve))
  res <- rbind(res, data.frame(m = m, loglik = L, var_explained = ve))
}
write.csv(res, file.path(OUT, "treemix_fit.csv"), row.names = FALSE)

cat("\n--- migration edges ---\n")
for (m in 1:5) {
  f <- file.path(OUT, sprintf("tm_m%d.edges.gz", m))
  tf <- file.path(OUT, sprintf("tm_m%d.treeout.gz", m))
  if (!file.exists(tf)) next
  ln <- readLines(gzfile(tf))
  mig <- ln[-1]
  if (length(mig)) {
    cat(sprintf("\nm = %d:\n", m))
    for (x in mig) cat("   ", x, "\n")
  }
}
cat("\nCAUTION: TreeMix assumes discrete populations. ADMIXTURE returned K = 1, so structure here\n")
cat("is clinal, and under isolation by distance TreeMix can fit edges that merely re-express the\n")
cat("cline. An edge is suggestive, not proof. The corroborating detail is WHERE it goes: an edge\n")
cat("crossing the Central Valley between cordilleras means more than one linking neighbours.\n")
