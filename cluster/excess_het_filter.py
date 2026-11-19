#!/usr/bin/env python3
"""Does removing excess-heterozygosity loci make individual relatedness usable?

Hall et al. (2020) recovered expected relatedness values in Pinus sylvestris GBS -- validated
against a known seed-orchard pedigree -- by removing loci with EXCESS HETEROZYGOSITY and
accounting for genotyping error. Collapsed paralogs generate persistent apparent heterozygotes,
because reads from two loci pile into one cluster, and that is the mechanism behind our residual
opposite-homozygote calls.

Sweeps the observed-heterozygosity ceiling and reports, at each level, the replicate error rate
AND its composition. Composition matters more than rate: dropout-type error (het-vs-hom)
attenuates relatedness and is correctable; non-homology (opposite homozygotes) is not.

Baselines: raw assembly 20.1% het-vs-hom / 79.8% opposite-hom at ~0.28 error.
           cr30 filter   63.3% / 36.7% at 0.0959 error.

Q2 becomes worth attempting if opposite-homozygote share keeps falling while enough loci survive
-- Gervais et al. (2019) put the floor for a usable relatedness matrix at >7,000 loci.
"""
import csv, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1] if len(sys.argv) > 1 else (
    f"{PROJ}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")

pairs = [(r["run_a"], r["run_b"]) for r in
         csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"), delimiter="\t")
         if r["status"] == "validated"]

samples, recs = [], []
for line in open(V):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
        continue
    recs.append([c.partition(":")[0].replace("|", "/") for c in f[9:]])
idx = {s: i for i, s in enumerate(samples)}

def het(g):
    a, _, b = g.partition("/")
    return a != b

Ho = np.zeros(len(recs))
for k, r in enumerate(recs):
    c = h = 0
    for g in r:
        if g in ("./.", "."):
            continue
        c += 1
        h += het(g)
    Ho[k] = h / c if c else 0.0

print(f"{V.split('/')[-1]}: {len(recs)} sites, {len(pairs)} validated pairs")
print(f"observed heterozygosity: median {np.median(Ho):.3f}, "
      f"90th pct {np.percentile(Ho, 90):.3f}, max {Ho.max():.3f}\n")

hdr = ("max Ho", "sites kept", "error rate", "het-vs-hom", "opp-hom")
print(f"  {hdr[0]:>7}{hdr[1]:>12}{hdr[2]:>12}{hdr[3]:>12}{hdr[4]:>10}")
for cut in (1.01, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10):
    keep = np.where(Ho <= cut)[0]
    n = m = hh = oh = 0
    for a, b in pairs:
        if a not in idx or b not in idx:
            continue
        ia, ib = idx[a], idx[b]
        for k in keep:
            ga, gb = recs[k][ia], recs[k][ib]
            if ga in ("./.", ".") or gb in ("./.", "."):
                continue
            n += 1
            if ga == gb:
                continue
            m += 1
            if het(ga) != het(gb):
                hh += 1
            elif not het(ga) and not het(gb):
                oh += 1
    lab = "none" if cut > 1 else f"{cut:.2f}"
    er = m / n if n else 0.0
    print(f"  {lab:>7}{len(keep):>12}{er:>12.4f}"
          f"{100 * hh / m if m else 0:>11.1f}%{100 * oh / m if m else 0:>9.1f}%")

print("\nGervais et al. (2019): a usable relatedness matrix needs >7,000 loci.")
print("Judge on opposite-homozygote share falling while loci stay above that floor.")
