#!/usr/bin/env python3
"""Hudson F_ST per stand pair from the reference calls, weighting TREES equally.

The earlier version pooled AD across a stand and used the read total as the
chromosome count. With ~90% PCR duplicates and several-fold depth differences
between libraries, that lets the deepest few trees carry a stand's allele
frequency, and F_ST came out at 0.40 with neighbouring Andean stands at 0.32 --
not credible. Here each tree contributes one frequency estimate, stands are the
unweighted mean over trees, and n is 2 x (trees with data), which is what the
de novo comparison values use.
"""
import gzip
from itertools import combinations

VCF = "/home/ecastron/araucaria/results/polyrad_fb/fb_AD.vcf.gz"
POP = "/home/ecastron/araucaria/metadata/popmap_177.tsv"
OUT = "/home/ecastron/araucaria/results/sv/fb_fst_indiv.tsv"
MINQUAL, MINMAF, MINSTANDS = 30.0, 0.05, 8
MINDP_TREE = 4      # reads needed before a tree gives a frequency
MINTREE = 4         # trees needed before a stand gives a frequency

stand = {}
for ln in open(POP):
    p = ln.split()
    if len(p) >= 2:
        stand[p[0]] = p[1]

acc, cols, nsite = {}, None, 0
with gzip.open(VCF, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            names = [s.split(".")[0] for s in f[9:]]
            order = sorted({stand[n] for n in names if n in stand})
            cols = {st: [i + 9 for i, n in enumerate(names)
                         if n in stand and stand[n] == st] for st in order}
            for a, b in combinations(order, 2):
                acc[(a, b)] = [0.0, 0.0, 0]
            print("stands %d | libraries %d" % (len(order),
                  sum(len(v) for v in cols.values())), flush=True)
            continue
        try:
            if float(f[5]) < MINQUAL:
                continue
        except ValueError:
            continue
        fmt = f[8].split(":")
        if "AD" not in fmt:
            continue
        adi = fmt.index("AD")
        cnt, graw, ntot = {}, [0, 0], 0
        for st, idx in cols.items():
            fr = []
            for i in idx:
                cell = f[i].split(":")
                if len(cell) <= adi:
                    continue
                ad = cell[adi].split(",")
                if len(ad) < 2:
                    continue
                try:
                    r, a = int(ad[0]), int(ad[1])
                except ValueError:
                    continue
                if r + a < MINDP_TREE:
                    continue
                fr.append(a / (r + a))          # one estimate per TREE
                graw[0] += r; graw[1] += a
            if len(fr) >= MINTREE:
                cnt[st] = (sum(fr) / len(fr), 2 * len(fr))
                ntot += 1
        if ntot < MINSTANDS or graw[0] + graw[1] == 0:
            continue
        gaf = graw[1] / (graw[0] + graw[1])
        if min(gaf, 1 - gaf) < MINMAF:
            continue
        nsite += 1
        for (a, b), A in acc.items():
            if a not in cnt or b not in cnt:
                continue
            p1, n1 = cnt[a]; p2, n2 = cnt[b]
            num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
            den = p1 * (1 - p2) + p2 * (1 - p1)
            if den > 0:
                A[0] += num; A[1] += den; A[2] += 1

COAST = {"PN", "TR", "VA"}
with open(OUT, "w") as fh:
    fh.write("stand_a\tstand_b\tfst\tnsites\tbetween\n")
    for (a, b), A in sorted(acc.items()):
        if A[1] > 0:
            fh.write("%s\t%s\t%.6f\t%d\t%d\n" % (a, b, A[0] / A[1], A[2],
                     int((a in COAST) != (b in COAST))))
print("sites used: %d | written %s" % (nsite, OUT))
