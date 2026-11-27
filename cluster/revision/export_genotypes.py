#!/usr/bin/env python
"""Revision (27 Nov 2026): export the genotype data analyzed in the paper for archiving.

Writes the 163 analyzed trees (results/Q1/feems/keep.txt) at the 78,026 biallelic SNPs of the
cr30 de novo VCF, with the same selection rules as structure.py and load_observed.py (single-base
REF and ALT). All FORMAT fields are kept; sample columns are renamed to the SRA run accession
(the ".1"/".2" library suffix removed), matching metadata/sra_runs_radseq.csv.
Output: results/revision/araucaria_163trees_78026snps.vcf.gz
"""
import gzip, re
P = "/home/ecastron/araucaria"
V = f"{P}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
site = {l.split("\t")[0]: l.split("\t")[1].strip() for l in open(f"{P}/metadata/popmap_177.tsv") if "\t" in l}
keep = {l.split()[0] for l in open(f"{P}/results/Q1/feems/keep.txt")}
n = 0
with open(V) as fi, gzip.open(f"{P}/results/revision/araucaria_163trees_78026snps.vcf.gz", "wt") as fo:
    for line in fi:
        if line.startswith("##"):
            fo.write(line); continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            cols = [i for i, s in enumerate(samples) if s in keep and s in site]
            assert len(cols) == 163, len(cols)
            fo.write("\t".join(f[:9] + [samples[i] for i in cols]) + "\n"); continue
        if len(f[3]) != 1 or len(f[4]) != 1 or "," in f[4]: continue
        fo.write("\t".join(f[:9] + [f[9 + i] for i in cols]) + "\n"); n += 1
assert n == 78026, n
print("SNPs", n)
