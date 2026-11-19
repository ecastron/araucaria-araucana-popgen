#!/usr/bin/env python3
"""Is population structure recoverable from the full 177-library assembly?

Individual genotype calls fail on this dataset (spec, Step 2 FINDING). The question that decides
how much of the design survives is whether STAND-level signal survives anyway: allele frequencies
pooled over 10-20 trees can be usable when individual calls are not, because errors independent
across individuals average out.

Usage: pop_structure_full.py <vcf>

Vectorised with numpy. The pure-Python version used for the 30-library preliminary took ~15 min
for 435 pairs; 177 libraries is 15,576 pairs and would have taken 9-12 h before permutation
tests. Genotypes are encoded once into an int8 matrix and every distance is an array operation.

THE CONFOUND: mismatch is dominated by sequencing depth (random-pair medians 0.7030 deep-deep,
0.4030 deep-shallow, 0.3994 shallow-shallow), so every test conditions on depth. The raw version
is printed alongside only to show what the confound alone would have produced.

PREDICTIONS, fixed before running:
  P1 positive control -- validated replicate pairs closer than non-replicates. If this fails,
     the matrix is wrong and nothing below is interpretable.
  P2 same-stand pairs closer than different-stand pairs, conditioning on depth.
  P3 SPLIT-HALF STABILITY -- the stand effect computed on two DISJOINT sets of trees agrees.

P3 is the one that decides it. P2 can reach significance on noise if some stands happen to hold
similar-depth libraries; only replication across disjoint individual sets rules that out.
"""
import csv, os, re, sys, random, statistics as st
import numpy as np

PROJ = "/home/ecastron/araucaria"
random.seed(17); np.random.seed(17)
VCF = sys.argv[1]

# ---- load and encode ---------------------------------------------------------------------
codes, rows, samples = {}, [], []
with open(VCF) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            continue
        r = np.empty(len(samples), dtype=np.int8)
        for k, cell in enumerate(f[9:]):
            g = cell.partition(":")[0].replace("|", "/")
            if g in ("./.", "."):
                r[k] = -1
            else:
                c = codes.get(g)
                if c is None:
                    c = codes[g] = len(codes)
                r[k] = c
        rows.append(r)
M = np.vstack(rows)                                   # sites x samples
print(f"{os.path.basename(VCF)}: {M.shape[0]} sites, {M.shape[1]} libraries, "
      f"{len(codes)} genotype states")

site_map = {}
for r in csv.reader(open(f"{PROJ}/metadata/popmap_177.tsv"), delimiter="\t"):
    if len(r) >= 2:
        site_map[r[0]] = r[1]

idx = {s: i for i, s in enumerate(samples)}
called = (M >= 0).sum(axis=0)
use = [s for s in samples if called[idx[s]] >= 0.005 * M.shape[0] and s in site_map]
cols = np.array([idx[s] for s in use])
stands = sorted({site_map[s] for s in use})
print(f"libraries used: {len(use)} | stands: {len(stands)}")
print("  " + "  ".join(f"{x}:{sum(1 for s in use if site_map[s]==x)}" for x in stands) + "\n")

Msub = M[:, cols]
callrate = (Msub >= 0).sum(axis=0) / M.shape[0]

valid = set()
for r in csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"), delimiter="\t"):
    if r["status"] == "validated":
        valid.add(frozenset((r["run_a"], r["run_b"])))

# ---- pairwise distances, vectorised ------------------------------------------------------
n = len(use)
I, J, G, DEP, SAME, REP = [], [], [], [], [], []
for i in range(n):
    a = Msub[:, i]
    va = a >= 0
    for j in range(i + 1, n):
        b = Msub[:, j]
        both = va & (b >= 0)
        cnt = int(both.sum())
        if cnt < 500:
            continue
        mm = int((a[both] != b[both]).sum())
        I.append(i); J.append(j)
        G.append(mm / cnt)
        DEP.append(abs(callrate[i] - callrate[j]) + callrate[i] + callrate[j])
        SAME.append(1 if site_map[use[i]] == site_map[use[j]] else 0)
        REP.append(1 if frozenset((use[i], use[j])) in valid else 0)

G = np.array(G); DEP = np.array(DEP)
SAME = np.array(SAME, bool); REP = np.array(REP, bool)
I = np.array(I); J = np.array(J)
print(f"pairwise comparisons: {len(G)}  (replicates {REP.sum()}, same-stand {SAME.sum()})")
print(f"depth vs genetic distance: r = {np.corrcoef(DEP, G)[0,1]:+.3f}  (the confound)")

Gr = G - np.polyval(np.polyfit(DEP, G, 1), DEP)       # depth regressed out

def perm(y, g, label, n_perm=9999):
    a, b = y[g], y[~g]
    if len(a) < 3 or len(b) < 3:
        print(f"  {label}: too few pairs"); return
    obs = a.mean() - b.mean()
    pool = y.copy(); cnt = 0
    na = len(a)
    for _ in range(n_perm):
        np.random.shuffle(pool)
        if pool[:na].mean() - pool[na:].mean() <= obs:
            cnt += 1
    print(f"  {label}: same={a.mean():.4f} diff={b.mean():.4f} delta={obs:+.4f} "
          f"p={(cnt+1)/(n_perm+1):.4f}")

print("\nP1 positive control (validated replicate pairs)")
perm(G, REP, "raw             "); perm(Gr, REP, "depth-controlled")
print("\nP2 same stand vs different stand")
perm(G, SAME, "raw             "); perm(Gr, SAME, "depth-controlled")

# ---- P3 split-half stability -------------------------------------------------------------
print("\nP3 split-half stability of the stand effect (50 disjoint splits)")
res = []
for _ in range(50):
    perm_idx = np.random.permutation(n)
    hA = set(perm_idx[: n // 2].tolist())
    inA = np.array([I[k] in hA and J[k] in hA for k in range(len(G))])
    inB = np.array([I[k] not in hA and J[k] not in hA for k in range(len(G))])
    def eff(mask):
        y, g = Gr[mask], SAME[mask]
        if g.sum() < 3 or (~g).sum() < 3:
            return None
        return y[g].mean() - y[~g].mean()
    ea, eb = eff(inA), eff(inB)
    if ea is not None and eb is not None:
        res.append((ea, eb))
if res:
    a = np.array([x for x, _ in res]); b = np.array([y for _, y in res])
    both_neg = int(((a < 0) & (b < 0)).sum())
    print(f"  half A mean effect {a.mean():+.4f} | half B mean effect {b.mean():+.4f}")
    print(f"  correlation between halves: r = {np.corrcoef(a, b)[0,1]:+.3f}")
    print(f"  splits where BOTH halves show same-stand-closer: {both_neg}/{len(res)}")
    if both_neg >= 0.8 * len(res):
        print("\n  >>> STAND EFFECT IS STABLE across disjoint halves: population-level structure")
        print("  >>> is recoverable. Q1/Q4/Q5 and a stand-level extended-phenotype test remain.")
    else:
        print("\n  >>> STAND EFFECT DOES NOT REPLICATE: population-level structure cannot be")
        print("  >>> relied on either, and the host side needs rescoping.")
else:
    print("  insufficient data for split-half test")
