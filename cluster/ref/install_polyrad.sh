#!/bin/bash
#SBATCH --job-name=instpolyrad
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/instpolyrad_%j.out
set -uo pipefail
module load R/R-4.4.0
export R_LIBS_USER=/home/ecastron/R/libs
export LDFLAGS="-L/usr/lib64 -L/lib64"
export LIBRARY_PATH="/usr/lib64:/lib64:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="/usr/lib64:/lib64:${LD_LIBRARY_PATH:-}"
Rscript --vanilla - <<'RS'
.libPaths(c("/home/ecastron/R/libs", .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (p in c("VariantAnnotation")) {
  if (!requireNamespace(p, quietly=TRUE))
    BiocManager::install(p, ask=FALSE, update=FALSE, lib="/home/ecastron/R/libs")
}
if (!requireNamespace("polyRAD", quietly=TRUE))
  install.packages("polyRAD", lib="/home/ecastron/R/libs")
cat("\n===== VERSIONS =====\n")
for (p in c("polyRAD","VariantAnnotation","pcaMethods","Rsamtools")) {
  ok <- requireNamespace(p, quietly=TRUE)
  cat(sprintf("%-20s %-6s %s\n", p, ok, if (ok) as.character(packageVersion(p)) else ""))
}
if (requireNamespace("polyRAD", quietly=TRUE)) {
  cat("\n===== key functions =====\n")
  for (f in c("VCF2RADdata","TestOverdispersion","IteratePopStructLD",
              "GetProbableGenotypes","GetWeightedMeanGenotypes","EstimateContaminationRate"))
    cat(sprintf("  %-28s %s\n", f, exists(f, where=asNamespace("polyRAD"))))
}
RS
