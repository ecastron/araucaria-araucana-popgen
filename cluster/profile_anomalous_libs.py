#!/usr/bin/env python3
"""What is actually wrong with the six libraries in the three non-replicate pairs?

Reclassifying them as "different trees" only makes sense if they BEHAVE like different trees.
They do not: unrelated libraries in this dataset mismatch at ~0.39, while these mismatch at
~0.69 -- worse than random, which no pair of ordinary A. araucana samples should manage. They
also carry 55k-102k shared sites against 10-25k for normal pairs: more data, worse agreement.

Per-library genotype composition separates the candidate explanations:

  excess heterozygosity      -> pooled or contaminated DNA (two individuals in one tube)
  excess homozygous-ALT      -> systematically different allele calls, e.g. divergent lineage
  ordinary composition       -> the anomaly is in the PAIRING, not the libraries

The distinction decides the action. Pooled or contaminated libraries must be dropped whatever
label they carry, because they are not one tree's genotype. A pure bookkeeping error can be
kept if the true identity is recoverable.
"""
import csv, glob, os, re, collections

PROJ = "/home/ecastron/araucaria"
VCF = sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf"))[0]

BAD = {"SRR11841736", "SRR11841792",   # CALM06
       "SRR11841723", "SRR11841740",   # CALR15
       "SRR11841688", "SRR11841786"}   # CATH04

samples, rows = [], []
with open(VCF) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            continue
        rows.append([c.split(":")[0].replace("|", "/") for c in f[9:]])

print(f"{os.path.basename(VCF)}: {len(samples)} libraries, {len(rows)} sites\n")
print(f"{'library':<15}{'flag':<8}{'called':>9}{'call%':>8}{'hom_ref':>9}{'het':>8}{'hom_alt':>9}{'het/called':>11}")

stats = {}
for i, s in enumerate(samples):
    c = collections.Counter(r[i] for r in rows)
    called = sum(v for k, v in c.items() if k not in ("./.", "."))
    if not called:
        continue
    hr = c.get("0/0", 0); ha = c.get("1/1", 0)
    het = called - hr - ha
    stats[s] = (called, hr, het, ha)
    flag = "BAD" if s in BAD else ""
    print(f"{s:<15}{flag:<8}{called:>9}{100*called/len(rows):>7.1f}%"
          f"{100*hr/called:>8.1f}%{100*het/called:>7.1f}%{100*ha/called:>8.1f}%{het/called:>11.3f}")

def summarise(sel, label):
    if not sel:
        return
    n = len(sel)
    ca = sum(stats[s][0] for s in sel) / n
    hr = sum(100*stats[s][1]/stats[s][0] for s in sel) / n
    he = sum(100*stats[s][2]/stats[s][0] for s in sel) / n
    ha = sum(100*stats[s][3]/stats[s][0] for s in sel) / n
    print(f"\n{label:<28} n={n:<3} mean called={ca:>9.0f}  hom_ref={hr:5.1f}%  het={he:5.1f}%  hom_alt={ha:5.1f}%")

bad = [s for s in stats if s in BAD]
ok = [s for s in stats if s not in BAD]
summarise(bad, "six anomalous libraries")
summarise(ok, "all other libraries")
