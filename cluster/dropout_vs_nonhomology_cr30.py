#!/usr/bin/env python3
"""Are the opposite homozygotes caused by bidirectional dropout, or by non-homologous loci?

CORRECTION TO EARLIER REASONING. I argued that opposite homozygotes (0/0 vs 1/1) prove
non-homologous loci, because dropout can only turn a heterozygote into a homozygote. That holds
WITHIN one library. Across TWO libraries it is false: a true A/B heterozygote can drop to A/A in
one and B/B in the other, producing an opposite homozygote from dropout alone. That error was
the basis for ruling out genotype likelihoods, so it has to be tested rather than assumed.

Supporting evidence for the dropout reading: median observed heterozygosity here is 0.016, an
order of magnitude below any plausible value for a wind-pollinated outcrossing conifer. We are
barely calling heterozygotes at all.

THE DISCRIMINATING TEST is depth. Bidirectional dropout requires BOTH libraries to sample only
one allele, whose probability collapses as depth rises. Non-homology does not care about depth.

  opp-hom share falls steeply with depth -> DROPOUT. Genotype likelihoods are the right fix and
                                            Q2 is viable; my earlier diagnosis was wrong.
  opp-hom share plateaus                 -> NON-HOMOLOGY. Earlier diagnosis stands.

Also reports heterozygosity by depth: under dropout, called heterozygosity must rise steeply
toward the true value as depth increases.
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
    row = []
    for c in f[9:]:
        p = c.split(":")
        g = p[0].replace("|", "/")
        try:
            d = int(p[1])
        except (IndexError, ValueError):
            d = 0
        row.append((g, d))
    recs.append(row)
idx = {s: i for i, s in enumerate(samples)}
print(f"{V.split('/')[-1]}: {len(recs)} sites, {len(pairs)} validated pairs\n")

def het(g):
    a, _, b = g.partition("/")
    return a != b

print("A. replicate disagreements by minimum depth of BOTH libraries")
print(f"  {'min DP':>7}{'sites':>10}{'error':>9}{'het-vs-hom':>12}{'opp-hom':>10}")
for cut in (1, 3, 5, 8, 12, 20, 30, 50):
    n = m = hh = oh = 0
    for a, b in pairs:
        if a not in idx or b not in idx:
            continue
        ia, ib = idx[a], idx[b]
        for r in recs:
            ga, da = r[ia]
            gb, db = r[ib]
            if ga in ("./.", ".") or gb in ("./.", ".") or da < cut or db < cut:
                continue
            n += 1
            if ga == gb:
                continue
            m += 1
            if het(ga) != het(gb):
                hh += 1
            elif not het(ga) and not het(gb):
                oh += 1
    if n:
        print(f"  {cut:>7}{n:>10}{m/n:>9.4f}"
              f"{100*hh/m if m else 0:>11.1f}%{100*oh/m if m else 0:>9.1f}%")

print("\nB. called heterozygosity by depth (all libraries)")
print(f"  {'depth':>10}{'genotypes':>12}{'het rate':>10}")
bins = [(1, 2), (3, 4), (5, 7), (8, 11), (12, 19), (20, 34), (35, 10**9)]
tot = {b: [0, 0] for b in bins}
for r in recs:
    for g, d in r:
        if g in ("./.", ".") or d < 1:
            continue
        for lo, hi in bins:
            if lo <= d <= hi:
                tot[(lo, hi)][0] += 1
                tot[(lo, hi)][1] += het(g)
                break
for lo, hi in bins:
    c, h = tot[(lo, hi)]
    if c:
        lab = f"{lo}-{hi}" if hi < 10**9 else f"{lo}+"
        print(f"  {lab:>10}{c:>12}{h/c:>10.4f}")
print("\n  A steeply rising het rate is the dropout signature. A flat one is not.")
