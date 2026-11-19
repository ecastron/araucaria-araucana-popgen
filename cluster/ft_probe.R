#!/usr/bin/env Rscript
# FEASIBILITY PROBE ONLY -- no hypothesis test. Measures how far the ITS ASV set gets
# through: genus assigned -> genus present in FungalTraits -> primary_lifestyle non-missing.
# Reported ASV-weighted AND read-weighted, because the two differ a lot and only the
# second is usually quoted.
suppressPackageStartupMessages({})
P <- "/home/ecastron/araucaria"
ft <- read.delim(file.path(P,"refs/fungaltraits_genus.tsv"), stringsAsFactors=FALSE, quote="")
ft$genus <- trimws(ft$genus)
ft$primary_lifestyle[ft$primary_lifestyle==""] <- NA
cat(sprintf("FungalTraits: %d genera, %d with primary_lifestyle\n\n", nrow(ft), sum(!is.na(ft$primary_lifestyle))))

probe <- function(sttag, txtag){
  fs <- file.path(P,"microbiome/ITS",sttag); ftx <- file.path(P,"microbiome/ITS",txtag)
  if(!file.exists(fs)||!file.exists(ftx)){cat("missing:",sttag,txtag,"\n");return(invisible())}
  st <- readRDS(fs); tax <- readRDS(ftx)
  cat(sprintf("========== %s / %s ==========\n", sttag, txtag))
  cat(sprintf("%d samples x %d ASVs, %s reads | tax ranks: %s\n",
      nrow(st), ncol(st), format(sum(st),big.mark=","), paste(colnames(tax),collapse=",")))
  if(!identical(rownames(tax), colnames(st))){cat("  !! taxonomy rownames do not match seqtab colnames\n"); return(invisible())}
  rd <- colSums(st); tot <- sum(rd)
  g <- tax[,"Genus"]
  gc <- sub("^g__","",as.character(g)); gc <- trimws(gc)
  gc[gc %in% c("","NA","unidentified")] <- NA
  cat(sprintf("  Genus assigned      : %5d/%5d ASVs (%5.1f%%)  %5.1f%% of reads\n",
      sum(!is.na(gc)), length(gc), 100*mean(!is.na(gc)), 100*sum(rd[!is.na(gc)])/tot))
  inft <- !is.na(gc) & gc %in% ft$genus
  cat(sprintf("  Genus in FungalTraits: %5d       ASVs (%5.1f%%)  %5.1f%% of reads\n",
      sum(inft), 100*mean(inft), 100*sum(rd[inft])/tot))
  pl <- rep(NA_character_,length(gc))
  pl[inft] <- ft$primary_lifestyle[match(gc[inft], ft$genus)]
  has <- !is.na(pl)
  cat(sprintf("  primary_lifestyle   : %5d       ASVs (%5.1f%%)  %5.1f%% of reads   <== GATE\n",
      sum(has), 100*mean(has), 100*sum(rd[has])/tot))
  cat(sprintf("  distinct lifestyles : %d\n", length(unique(pl[has]))))
  tb <- tapply(rd[has], pl[has], sum); tb <- sort(tb, decreasing=TRUE)
  na <- tapply(rep(1,sum(has)), pl[has], sum); na <- na[names(tb)]
  cat(sprintf("  %-28s %6s %8s\n","lifestyle","ASVs","%reads"))
  for(i in seq_along(tb)) cat(sprintf("  %-28s %6d %8.2f\n", names(tb)[i], na[i], 100*tb[i]/tot))
  cat(sprintf("  %-28s %6d %8.2f\n","(NO LIFESTYLE)", sum(!has), 100*sum(rd[!has])/tot))
  cat("\n  top 12 unassigned-at-genus lineages by reads:\n")
  fin <- rep("k__Fungi", ncol(st))
  for(r in c("Kingdom","Phylum","Class","Order","Family")) if(r %in% colnames(tax)){
    v <- as.character(tax[,r]); ok <- !is.na(v) & v!=""; fin[ok] <- paste0(substr(tolower(r),1,1),"__",v[ok]) }
  u <- !has; tu <- sort(tapply(rd[u], fin[u], sum), decreasing=TRUE)[1:12]
  for(i in seq_along(tu)) if(!is.na(tu[i])) cat(sprintf("    %-40s %6.2f%% of all reads\n", names(tu)[i], 100*tu[i]/tot))
  cat("\n")
}
probe("seqtab_nochim.rds","taxonomy.rds")
probe("seqtab_nochim_r1.rds","taxonomy_r1.rds")
