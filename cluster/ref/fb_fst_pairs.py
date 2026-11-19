#!/usr/bin/env python3
"""Hudson F_ST per stand pair from the REFERENCE-based FreeBayes calls.

Read-based throughout (AD counts, never a genotype call), so it is directly
comparable with the de novo read-based values in figdata/host_fst_pairs.tsv and
is immune to the heterozygote artifact. The point is not the absolute number but
whether an independent aligner + caller recovers the same among-stand pattern.
"""
import gzip, sys
from itertools import combinations

VCF = "/home/ecastron/araucaria/results/polyrad_fb/fb_AD.vcf.gz"
POP = "/home/ecastron/araucaria/metadata/popmap_177.tsv"
OUT = "/home/ecastron/araucaria/results/sv/fb_fst_pairs_filt.tsv"
MINCHR = 6
# --pooled-continuous reports every alt allele it sees, sequencing error included:
# 54.8% of raw sites have global MAF < 0.05 and QUAL's 5th percentile is 0. The de
# novo set this is being compared against went through ipyrad's filters and a 30%
# call-rate floor, so an unfiltered comparison measures the filtering, not the
# pipeline. Match it.
MINQUAL = 30.0
MINMAF = 0.05
MINSTANDS = 8

stand = {}
for ln in open(POP):
    p = ln.split()
    if len(p) >= 2:
        stand[p[0]] = p[1]

acc = {}
cols = None
nsite = 0
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
            print("stands %d | libraries matched %d"
                  % (len(order), sum(len(v) for v in cols.values())), flush=True)
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
        cnt = {}
        for st, idx in cols.items():
            ref = alt = 0
            for i in idx:
                cell = f[i].split(":")
                if len(cell) <= adi:
                    continue
                ad = cell[adi].split(",")
                if len(ad) < 2:
                    continue
                try:
                    ref += int(ad[0]); alt += int(ad[1])
                except ValueError:
                    continue
            cnt[st] = (ref, alt)
        gr = sum(v[0] for v in cnt.values()); ga = sum(v[1] for v in cnt.values())
        if gr + ga == 0:
            continue
        gaf = ga / (gr + ga)
        if min(gaf, 1 - gaf) < MINMAF:
            continue
        if sum(1 for v in cnt.values() if v[0] + v[1] >= MINCHR) < MINSTANDS:
            continue
        nsite += 1
        for (a, b), A in acc.items():
            r1, a1 = cnt[a]; r2, a2 = cnt[b]
            n1, n2 = r1 + a1, r2 + a2
            if n1 < MINCHR or n2 < MINCHR:
                continue
            p1, p2 = a1 / n1, a2 / n2
            num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
            den = p1 * (1 - p2) + p2 * (1 - p1)
            if den > 0:
                A[0] += num; A[1] += den; A[2] += 1

COAST = {"PN", "TR", "VA"}
with open(OUT, "w") as fh:
    fh.write("stand_a\tstand_b\tfst\tnsites\tbetween\n")
    for (a, b), A in sorted(acc.items()):
        if A[1] <= 0:
            continue
        fh.write("%s\t%s\t%.6f\t%d\t%d\n"
                 % (a, b, A[0] / A[1], A[2], int((a in COAST) != (b in COAST))))
print("records scanned: %d | written %s" % (nsite, OUT))
