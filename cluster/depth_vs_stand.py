#!/usr/bin/env python3
"""Is sequencing depth confounded with stand?

The stand test conditions on depth, because mismatch is dominated by it. But if depth varies
systematically BY STAND -- e.g. because a stand's samples were prepped or sequenced together --
then regressing depth out also removes the stand signal, and the negative P2/P3 result would be
an artefact of the control rather than a property of the data.

This must be ruled out before the negative result is accepted. One-way F on per-library call
rate grouped by stand, with a permutation p-value.

  p large  -> depth is spread evenly across stands; controlling for it is safe and the negative
              stand result stands.
  p small  -> depth and stand are entangled; the depth control was removing real signal and the
              stand effect must be re-estimated some other way (e.g. matched subsampling).
"""
import csv, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1]

samples, rows = [], []
for line in open(V):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
        continue
    rows.append(np.array([0 if c.partition(":")[0] in ("./.", ".") else 1
                          for c in f[9:]], dtype=np.int8))
M = np.vstack(rows)
cr = M.mean(axis=0)

site = {r[0]: r[1] for r in csv.reader(open(f"{PROJ}/metadata/popmap_177.tsv"), delimiter="\t")
        if len(r) >= 2}
idx = {s: i for i, s in enumerate(samples)}
use = [s for s in samples if s in site and cr[idx[s]] >= 0.005]
stands = sorted({site[s] for s in use})
print(f"{len(use)} libraries, {len(stands)} stands\n")

hdr = ("stand", "n", "mean call rate", "sd", "min", "max")
print(f"  {hdr[0]:<7}{hdr[1]:>4}{hdr[2]:>16}{hdr[3]:>9}{hdr[4]:>9}{hdr[5]:>9}")
groups = []
for stx in stands:
    v = np.array([cr[idx[s]] for s in use if site[s] == stx])
    groups.append(v)
    print(f"  {stx:<7}{len(v):>4}{v.mean():>16.4f}{v.std():>9.4f}{v.min():>9.4f}{v.max():>9.4f}")

allv = np.concatenate(groups)
gsz = [len(g) for g in groups]

def F(gs):
    gm = np.concatenate(gs).mean()
    k = len(gs)
    between = sum(len(g) * (g.mean() - gm) ** 2 for g in gs) / (k - 1)
    within = sum(((g - g.mean()) ** 2).sum() for g in gs) / (len(allv) - k)
    return between / within if within > 0 else 0.0

obs = F(groups)
rng = np.random.default_rng(3)
cnt = 0
for _ in range(9999):
    p = rng.permutation(allv)
    gs, k = [], 0
    for n in gsz:
        gs.append(p[k:k + n]); k += n
    if F(gs) >= obs:
        cnt += 1
pval = (cnt + 1) / 10000
print(f"\n  call rate differs by stand: F = {obs:.3f}, permutation p = {pval:.4f}")
if pval < 0.05:
    print("\n  >>> DEPTH IS CONFOUNDED WITH STAND. Regressing out depth also removes stand")
    print("  >>> signal, so the negative P2/P3 is an artefact of the control, not a result.")
    print("  >>> The stand effect must be re-estimated by depth-matched subsampling instead.")
else:
    print("\n  >>> Depth is spread evenly across stands. Controlling for it does not remove")
    print("  >>> stand signal, so the negative stand result stands as a genuine finding.")
