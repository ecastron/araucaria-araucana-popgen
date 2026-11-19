#!/usr/bin/env python3
"""For every library in a declared replicate pair, find its closest library in the whole set.

Three declared pairs (CALM06, CALR15, CATH04) mismatch at ~0.69 -- worse than the ~0.39 median
between RANDOM libraries. Two libraries of the same tree cannot be more different than two
libraries of different trees, so something in those pairs is wrong rather than merely noisy.

If a library's nearest neighbour is some OTHER library rather than its declared partner, the
labels are swapped and the true pairing is recoverable. If instead nothing in the set is close
to it, the library is an outlier in its own right (contamination, a different species, or a
failed prep) and it should be dropped rather than repaired.

This distinguishes a fixable bookkeeping error from an unfixable data problem, which the
aggregate mismatch statistic cannot do.
"""
import csv, glob, os, re

PROJ = "/home/ecastron/araucaria"
VCF = sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf"))[0]

samples, cols = [], []
with open(VCF) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            continue
        cols.append([c.split(":")[0].replace("|", "/") for c in f[9:]])
print(f"{os.path.basename(VCF)}: {len(samples)} libraries, {len(cols)} sites\n")

def mm(ia, ib):
    n = m = 0
    for r in cols:
        a, b = r[ia], r[ib]
        if a in ("./.", ".") or b in ("./.", "."):
            continue
        n += 1
        if a != b:
            m += 1
    return (m / n if n >= 500 else None), n

idx = {s: i for i, s in enumerate(samples)}
pairs = {r["tree_id"]: (r["run_a"], r["run_b"])
         for r in csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs.tsv"), delimiter="\t")}
declared = {}
for t, (a, b) in pairs.items():
    declared[a] = (b, t)
    declared[b] = (a, t)

focus = list(pairs.keys())   # all declared pairs
for t in focus:
    if t not in pairs:
        continue
    for lib in pairs[t]:
        if lib not in idx:
            continue
        scores = []
        for other in samples:
            if other == lib:
                continue
            e, n = mm(idx[lib], idx[other])
            if e is not None:
                scores.append((e, other, n))
        scores.sort()
        partner, _ = declared[lib]
        drank = next((k for k, (e, o, n) in enumerate(scores, 1) if o == partner), None)
        dval = next((e for e, o, n in scores if o == partner), None)
        print(f"{t}  {lib}")
        for e, o, n in scores[:3]:
            tag = "  <-- DECLARED PARTNER" if o == partner else ""
            print(f"      best: {o}  mismatch {e:.4f}  ({n} sites){tag}")
        if drank and drank > 3:
            print(f"      declared partner {partner}: mismatch {dval:.4f}, ranked {drank} of {len(scores)}")
        print()
