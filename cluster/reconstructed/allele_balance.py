#!/usr/bin/env python3
"""Allele-balance distribution of high-depth genotypes in the de novo calls (Figure S2A).

PRODUCES  results/reconstructed/allele_balance.txt
          (original: results/Q1/recal/allele_balance.txt, 24 Aug 2026, copied to
          figdata/allele_balance.txt on 26 Aug; parsed by figures/build_supp_journal.py for S2A)

INPUT     rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf   (FORMAT GT:DP:CATG, 177 libraries)

METHOD    For every SNP whose REF and first ALT are both in C/A/T/G and differ, every library cell is
          read from its CATG base counts (the genotype call is ignored, so uncalled ./. cells with
          reads are included). n = REF count + first-ALT count; a cell enters when n >= 30.
          Balance = ALT / n. Histogram in 20 bins of width 0.05, bin = int(balance * 20) with 1.0
          folded into the last bin. "Putative heterozygotes" are 0.15 < balance < 0.85; SD is the
          population SD. The overdispersion factor divides that SD by the fair-binomial SD at n = 30
          after rounding it to the 0.091 that is printed (0.1919 / 0.091 = 2.11; unrounded it would
          be 2.10). In the asymmetry line "above" includes balance exactly 0.5 (1,837 genotypes).
          These conventions are inferred: they are the ones that reproduce the original exactly
          (binning by division, balance / 0.05, misplaces float boundaries and does not).

PROVENANCE  Reconstructed 18-19 Nov 2026 to replace lost inline code (24 Aug 2026 session, written
          next to plateau_diagnosis.txt and clean_loci.py output in results/Q1/recal/).
VERIFICATION  (run 18 Nov 2026 on McCartney, srun -p basic) BYTE-IDENTICAL (cmp) to
          results/Q1/recal/allele_balance.txt; md5 80e570d5d104004013721b29fd375240, same as the
          local figdata/allele_balance.txt.
"""
import numpy as np

P = "/home/ecastron/araucaria"
VCF = f"{P}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
OUT = f"{P}/results/reconstructed/allele_balance.txt"
MINDP, B = 30, "CATG"

bal = []
for line in open(VCF):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    r, a = f[3], f[4].split(",")[0]
    if r not in B or a not in B or r == a:
        continue
    ri, ai = B.index(r), B.index(a)
    for cell in f[9:]:
        p = cell.split(":")
        try:
            g = [int(x) for x in p[2].split(",")]
        except (IndexError, ValueError):
            continue
        n = g[ri] + g[ai]
        if n >= MINDP:
            bal.append(g[ai] / n)
b = np.array(bal); N = len(b)

h = np.bincount(np.minimum((b * 20).astype(int), 19), minlength=20)
het = b[(b > 0.15) & (b < 0.85)]
below, above = int((het < 0.5).sum()), int((het >= 0.5).sum())
exp_sd = round(np.sqrt(0.25 / MINDP), 3)          # 0.091, as printed

L = [f"genotypes at depth >={MINDP}x: {N}", "",
     "allele-balance histogram (fraction ALT), bin width 0.05:"]
for k, c in enumerate(h):
    L.append(f"  {k*0.05:.2f}-{(k+1)*0.05:.2f} {c:>8} ({100*c/N:5.2f}%) {'#'*int(60*c/h.max())}")
L += ["",
      f"putative heterozygotes (0.15 < balance < 0.85): {len(het)} "
      f"({100*len(het)/N:.2f}% of high-depth genotypes)",
      f"  mean balance = {het.mean():.4f}   (0.5 expected; deviation => reference bias)",
      f"  SD           = {het.std():.4f}   (fair binomial at n={MINDP} predicts ~{exp_sd:.3f})",
      f"  overdispersion factor vs fair binomial: {het.std()/exp_sd:.2f}x",
      f"  asymmetry: {below} below 0.5 vs {above} above  "
      f"({100*below/len(het):.1f}% / {100*above/len(het):.1f}%)",
      "",
      f"homozygote spikes: balance<=0.05 = {100*(b<=0.05).mean():.2f}%, "
      f">=0.95 = {100*(b>=0.95).mean():.2f}%"]
open(OUT, "w").write("\n".join(L) + "\n")
print(f"written -> {OUT}")
