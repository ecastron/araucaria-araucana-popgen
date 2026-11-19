#!/usr/bin/env python3
"""Are the "replicate pairs" the SAME tree, or two related trees?

Eduardo's hypothesis, and it fits the data better than mine. Everything attributed to genotyping
failure is equally explained by the pairs being close relatives rather than the same individual:

  - zero minority reads in both libraries -> two genuinely homozygous individuals differing
  - opposite homozygotes at 3.4% of comparisons -> normal for two related individuals
  - pair distance 0.096 vs 0.14 for random pairs -> intermediate, i.e. where RELATIVES sit
  - mutual nearest-neighbour 10/10 -> close relatives in the same stand would also do this
  - MISMATCH IS FLAT ACROSS DEPTH (0.0965 at >=1x, 0.1026 at >=30x) -> real genetic difference.
    Genotyping error would DECLINE with depth. This is the strongest evidence and I misread it.

THE TEST: estimate relatedness on high-depth sites only, where calls are reliable.

    same individual        r ~ 1.0
    parent-offspring/full-sib r ~ 0.5
    half-sib               r ~ 0.25
    unrelated              r ~ 0.0

Calibrated against the empirical distribution of r among ALL pairs in the dataset, so the scale
is anchored in this data rather than in theory.

If the pairs are relatives, there is no genotyping catastrophe: the loci and calls are fine, and
the case against Q2/Q3/Q6 collapses. If they are identical, the residual really is error.
"""
import csv, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1] if len(sys.argv) > 1 else (
    f"{PROJ}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")
MINDP = 20          # high-depth only: at 20x calls are reliable, so residual difference is real

pairs = [(r["tree_id"], r["run_a"], r["run_b"]) for r in
         csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"), delimiter="\t")
         if r["status"] in ("validated", "NOT_REPLICATE", "QUESTIONABLE")]
status = {r["tree_id"]: r["status"] for r in
          csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"), delimiter="\t")}

samples, G, D = [], [], []
for line in open(V):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
        continue
    g = np.full(len(samples), np.nan)
    d = np.zeros(len(samples), dtype=np.int32)
    for k, c in enumerate(f[9:]):
        p = c.split(":")
        gt = p[0].replace("|", "/")
        try:
            d[k] = int(p[1])
        except (IndexError, ValueError):
            d[k] = 0
        if gt in ("./.", "."):
            continue
        a, _, b = gt.partition("/")
        g[k] = (a != "0") + (b != "0")          # ALT allele dosage 0/1/2
    G.append(g); D.append(d)
G = np.vstack(G); D = np.vstack(D)
idx = {s: i for i, s in enumerate(samples)}
print(f"{V.split('/')[-1]}: {G.shape[0]} sites, {G.shape[1]} libraries")

# Allele frequencies from high-depth calls only.
hd = D >= MINDP
gg = np.where(hd, G, np.nan)
p = np.nanmean(gg, axis=1) / 2
usable = (p > 0.05) & (p < 0.95) & (np.sum(~np.isnan(gg), axis=1) >= 30)
print(f"usable polymorphic sites at >={MINDP}x: {usable.sum()}\n")
gg = gg[usable]; p = p[usable]

def ritland(i, j):
    x, y = gg[:, i], gg[:, j]
    ok = ~np.isnan(x) & ~np.isnan(y)
    if ok.sum() < 200:
        return np.nan, int(ok.sum())
    xi, yi, pi = x[ok] / 2, y[ok] / 2, p[ok]
    v = ((xi - pi) * (yi - pi)) / (pi * (1 - pi))
    return 2 * float(np.mean(v)), int(ok.sum())

# Empirical scale: relatedness among all pairs, to anchor what 0 and 1 look like here.
n = len(samples)
rng = np.random.default_rng(3)
bg = []
for _ in range(2000):
    i, j = rng.integers(0, n, 2)
    if i == j:
        continue
    r, k = ritland(i, j)
    if not np.isnan(r):
        bg.append(r)
bg = np.array(bg)
print(f"background relatedness among random library pairs (n={len(bg)}):")
print(f"  mean {bg.mean():+.4f}  median {np.median(bg):+.4f}  "
      f"5th {np.percentile(bg,5):+.4f}  95th {np.percentile(bg,95):+.4f}\n")

# Self-comparison ceiling: a library against itself must give r ~ 1 on this scale.
selfs = [ritland(idx[s], idx[s])[0] for s in samples[:20]]
print(f"library vs ITSELF (the r=1 anchor): mean {np.nanmean(selfs):.4f}\n")

print(f"  {'tree':<10}{'status':<15}{'sites':>8}{'Ritland r':>11}   interpretation")
for t, a, b in pairs:
    if a not in idx or b not in idx:
        continue
    r, k = ritland(idx[a], idx[b])
    if np.isnan(r):
        print(f"  {t:<10}{status.get(t,''):<15}{k:>8}{'':>11}   too few shared sites")
        continue
    if r > 0.80:
        lab = "SAME INDIVIDUAL"
    elif r > 0.35:
        lab = "parent-offspring / full-sib"
    elif r > 0.15:
        lab = "half-sib / cousin"
    else:
        lab = "unrelated"
    print(f"  {t:<10}{status.get(t,''):<15}{k:>8}{r:>11.4f}   {lab}")

print("\nIf these land near 0.5 rather than 1.0, the pairs are RELATIVES, not replicates:")
print("the loci and genotype calls are fine and the case against Q2/Q3/Q6 collapses.")
