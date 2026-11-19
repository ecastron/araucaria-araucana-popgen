#!/usr/bin/env python3
"""Why did the positive control fail on the full assembly?

On the 30-library subset, validated replicate pairs separated cleanly (delta -0.102, p 0.0004).
On the full 177-library assembly the same test gives delta -0.053, p 0.144 -- a failure. Since
the design says nothing downstream is interpretable when P1 fails, the cause has to be
established before the stand result can be read at all.

Two explanations with opposite consequences:

  LOW POWER      only 11 replicate pairs against 15,206 other comparisons, and a difference-in-
                 means permutation test is weak at that imbalance. The signal would be real and
                 the test blind to it.
  SIGNAL LOSS    min_samples_locus = 4 over 177 samples retains a far larger locus set (3.19M
                 sites vs 654k), possibly dominated by the non-homologous clusters already
                 identified. The signal would genuinely be gone.

The rank test separates them, and does not depend on group-mean power: is each replicate its
partner's NEAREST NEIGHBOUR among all 174 alternatives? That test gave 10/10 mutual rank-1 on
the subset. Rank is also immune to the depth confound shifting the whole distribution, because
it asks only about ordering within a library's own row.

  still mutual rank-1  -> identity signal present; P1's failure is test design, and the stand
                          result can be reconsidered with an appropriately powered test
  ranks now scattered  -> signal genuinely lost in this locus set; nothing here is usable
"""
import csv, os, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
VCF = sys.argv[1]

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
M = np.vstack(rows)
idx = {s: i for i, s in enumerate(samples)}
print(f"{os.path.basename(VCF)}: {M.shape[0]} sites, {M.shape[1]} libraries")

pairs = []
for r in csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"), delimiter="\t"):
    if r["status"] == "validated" and r["run_a"] in idx and r["run_b"] in idx:
        pairs.append((r["tree_id"], r["run_a"], r["run_b"]))
print(f"validated pairs present in this VCF: {len(pairs)}\n")

def row_dists(a_i):
    """Distance from library a to every other library, vectorised."""
    a = M[:, a_i]
    va = a >= 0
    out = np.full(M.shape[1], np.nan)
    for j in range(M.shape[1]):
        if j == a_i:
            continue
        both = va & (M[:, j] >= 0)
        n = int(both.sum())
        if n < 500:
            continue
        out[j] = (a[both] != M[both, j]).sum() / n
    return out

print(f"  {'tree':<9}{'library':<14}{'partner rank':>13}{'pair dist':>11}{'best dist':>11}{'median':>9}")
ranks = []
for t, a, b in pairs:
    for x, y in ((a, b), (b, a)):
        d = row_dists(idx[x])
        valid = ~np.isnan(d)
        order = np.argsort(np.where(valid, d, np.inf))
        rank = int(np.where(order == idx[y])[0][0]) + 1
        ranks.append(rank)
        print(f"  {t:<9}{x:<14}{rank:>13}{d[idx[y]]:>11.4f}"
              f"{np.nanmin(d):>11.4f}{np.nanmedian(d):>9.4f}")

ranks = np.array(ranks)
n_lib = M.shape[1] - 1
print(f"\nrank-1 (partner is nearest neighbour): {(ranks == 1).sum()} / {len(ranks)}")
print(f"rank <= 3:  {(ranks <= 3).sum()} / {len(ranks)}   median rank {np.median(ranks):.0f} "
      f"of ~{n_lib} candidates")
# Under no signal, rank is uniform: P(rank == 1) = 1/n_lib per library.
exp = len(ranks) / n_lib
print(f"expected rank-1 count if there were NO identity signal: {exp:.2f}")
if (ranks == 1).sum() >= 0.8 * len(ranks):
    print("\n>>> IDENTITY SIGNAL INTACT. P1's failure is a test-power artefact, not signal loss:")
    print(">>> the group-mean permutation test is simply blind at 11 pairs vs 15,206.")
elif (ranks == 1).sum() > 3 * exp:
    print("\n>>> PARTIAL: identity signal present but degraded relative to the subset.")
else:
    print("\n>>> IDENTITY SIGNAL LOST in this locus set. Nothing downstream is usable.")
