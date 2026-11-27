#!/usr/bin/env Rscript
# Revision (27 Nov 2026): export the ITS1 ASV table analyzed in the paper (forward-read run,
# dada2_its_r1.R) for archiving: read counts per library (rows) and ASV (columns), and the UNITE
# taxonomy of each ASV. ASVs are named ASV1, ASV2, ... with their sequences in the taxonomy file.
# Outputs: results/revision/its1_asv_counts.csv.gz, results/revision/its1_asv_taxonomy.csv.gz
P <- "/home/ecastron/araucaria"
st <- readRDS(file.path(P, "microbiome/ITS/seqtab_nochim_r1.rds"))
tx <- readRDS(file.path(P, "microbiome/ITS/taxonomy_r1.rds"))
stopifnot(identical(colnames(st), rownames(tx)))
ids <- paste0("ASV", seq_len(ncol(st)))
counts <- data.frame(run_accession = rownames(st), setNames(as.data.frame(st), ids), check.names = FALSE)
taxa <- data.frame(asv = ids, sequence = colnames(st), as.data.frame(tx), check.names = FALSE)
write.csv(counts, gzfile(file.path(P, "results/revision/its1_asv_counts.csv.gz")), row.names = FALSE)
write.csv(taxa, gzfile(file.path(P, "results/revision/its1_asv_taxonomy.csv.gz")), row.names = FALSE)
cat("libraries", nrow(st), "ASVs", ncol(st), "\n")
