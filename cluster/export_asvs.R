#!/usr/bin/env Rscript
# Export DADA2 ASV sequences to FASTA, one file per amplicon.
#
# DADA2 stores the ASV sequence itself as the column name of seqtab. Downstream tools need
# stable short IDs instead, so we hash each sequence to an MD5 and keep a lookup table. MD5 of
# the sequence is the same convention QIIME 2 uses, which makes the two interchangeable and
# means the ID is reproducible from the sequence alone rather than depending on column order.

suppressPackageStartupMessages({library(digest)})

PROJ <- "/home/ecastron/araucaria"

for (AMP in c("16S", "ITS")) {
  f <- file.path(PROJ, "microbiome", AMP, "seqtab_nochim.rds")
  if (!file.exists(f)) { cat(sprintf("[%s] seqtab not ready, skipping\n", AMP)); next }

  seqtab <- readRDS(f)
  seqs <- colnames(seqtab)
  ids  <- vapply(seqs, function(s) digest(s, algo = "md5", serialize = FALSE), character(1))
  stopifnot(!anyDuplicated(ids))

  fa <- file.path(PROJ, "microbiome", AMP, paste0("asv_", AMP, ".fasta"))
  writeLines(as.vector(rbind(paste0(">", ids), seqs)), fa)

  # Rewrite the count table under the same IDs so the tree tips and the table agree.
  colnames(seqtab) <- ids
  saveRDS(seqtab, file.path(PROJ, "microbiome", AMP, "seqtab_nochim_md5.rds"))
  write.csv(data.frame(asv_id = ids, sequence = seqs, nchar = nchar(seqs), row.names = NULL),
            file.path(PROJ, "microbiome", AMP, "asv_lookup.csv"), row.names = FALSE)

  cat(sprintf("[%s] %d ASVs, length min=%d median=%d max=%d -> %s\n",
              AMP, length(seqs), min(nchar(seqs)), median(nchar(seqs)), max(nchar(seqs)), fa))
}
