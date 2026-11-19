#!/usr/bin/env python3
"""Settle the F_ST estimator discrepancy on ONE identical SNP set.

Three quantities per stand pair, all ratio-of-averages over the same sites:
  hudson  Hudson (Bhatia et al. 2013) - allele frequencies only
  wc_obs  Weir & Cockerham 1984 theta - uses OBSERVED heterozygosity
  wc_hwe  the same theta with h_bar replaced by its HWE expectation 2p(1-p)

wc_obs vs wc_hwe isolates how much of the gap is caused by this dataset's
documented heterozygote deficit contaminating the W&C variance components.
If wc_hwe lands on hudson, the deficit is the whole explanation and W&C is
simply not usable here.
"""
import sys, math
from itertools import combinations

vcf, popmap, out = sys.argv[1], sys.argv[2], sys.argv[3]
MINCHR = int(sys.argv[4]) if len(sys.argv) > 4 else 6

stand = {}
for ln in open(popmap):
    p = ln.split()
    if len(p) >= 2:
        stand[p[0]] = p[1]

cols, order = None, None
acc = {}          # pair -> [hud_num, hud_den, a, abc, a_hwe, abc_hwe, nsites]

with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            names = [s.split('.')[0] for s in f[9:]]
            order = sorted(set(stand[n] for n in names if n in stand))
            cols = {st: [i + 9 for i, n in enumerate(names)
                         if n in stand and stand[n] == st] for st in order}
            for a, b in combinations(order, 2):
                acc[(a, b)] = [0.0] * 7
            continue
        # per-stand counts from GT only
        cnt = {}
        for st in order:
            nalt = nobs = nhet = 0
            for i in cols[st]:
                g = f[i].split(':', 1)[0]
                if g[0] == '.':
                    continue
                a1, a2 = g[0], g[2] if len(g) > 2 else g[0]
                if a1 not in '01' or a2 not in '01':
                    continue
                nobs += 1
                nalt += (a1 == '1') + (a2 == '1')
                nhet += (a1 != a2)
            cnt[st] = (nobs, nalt, nhet)
        for a, b in acc:
            n1, k1, h1 = cnt[a]; n2, k2, h2 = cnt[b]
            if 2 * n1 < MINCHR or 2 * n2 < MINCHR:
                continue
            p1, p2 = k1 / (2 * n1), k2 / (2 * n2)
            A = acc[(a, b)]
            # --- Hudson
            num = (p1 - p2) ** 2 - p1 * (1 - p1) / (2 * n1 - 1) - p2 * (1 - p2) / (2 * n2 - 1)
            den = p1 * (1 - p2) + p2 * (1 - p1)
            if den > 0:
                A[0] += num; A[1] += den
            # --- Weir & Cockerham theta, r = 2
            r = 2.0
            nbar = (n1 + n2) / r
            if nbar <= 1:
                continue
            nc = (r * nbar - (n1 * n1 + n2 * n2) / (r * nbar)) / (r - 1)
            pbar = (n1 * p1 + n2 * p2) / (r * nbar)
            s2 = (n1 * (p1 - pbar) ** 2 + n2 * (p2 - pbar) ** 2) / ((r - 1) * nbar)
            for hb, oa, ob in ((( n1 * (h1 / n1) + n2 * (h2 / n2)) / (r * nbar), 2, 3),
                               ((n1 * 2 * p1 * (1 - p1) + n2 * 2 * p2 * (1 - p2)) / (r * nbar), 4, 5)):
                aa = (nbar / nc) * (s2 - (1 / (nbar - 1)) *
                                    (pbar * (1 - pbar) - ((r - 1) / r) * s2 - hb / 4))
                bb = (nbar / (nbar - 1)) * (pbar * (1 - pbar) - ((r - 1) / r) * s2 -
                                            ((2 * nbar - 1) / (4 * nbar)) * hb)
                cc = hb / 2
                A[oa] += aa; A[ob] += aa + bb + cc
            A[6] += 1

COAST = {'PN', 'TR', 'VA'}
rows = []
for (a, b), A in sorted(acc.items()):
    if A[6] == 0:
        continue
    rows.append((a, b, A[0] / A[1] if A[1] else float('nan'),
                 A[2] / A[3] if A[3] else float('nan'),
                 A[4] / A[5] if A[5] else float('nan'),
                 int(A[6]), int((a in COAST) != (b in COAST))))
with open(out, 'w') as fh:
    fh.write("stand_a\tstand_b\thudson\twc_obs\twc_hwe\tnsites\tbetween\n")
    for r_ in rows:
        fh.write("%s\t%s\t%.6f\t%.6f\t%.6f\t%d\t%d\n" % r_)

def mean(v):
    v = [x for x in v if x == x]
    return sum(v) / len(v) if v else float('nan')

print("pairs: %d | identical SNP set, ratio-of-averages, min %d chromosomes/stand\n"
      % (len(rows), MINCHR))
print("%-10s %9s %9s %9s" % ("", "hudson", "wc_obs", "wc_hwe"))
for lab, sel in (("all pairs", lambda r: True),
                 ("between-region", lambda r: r[6] == 1),
                 ("within-region", lambda r: r[6] == 0)):
    s = [r for r in rows if sel(r)]
    print("%-10s %9.4f %9.4f %9.4f   (n=%d)"
          % (lab, mean([r[2] for r in s]), mean([r[3] for r in s]),
             mean([r[4] for r in s]), len(s)))
bw = [r for r in rows if r[6] == 1]; wi = [r for r in rows if r[6] == 0]
print("\nbetween/within ratio:")
for i, lab in ((2, "hudson"), (3, "wc_obs"), (4, "wc_hwe")):
    print("  %-8s %.2fx" % (lab, mean([r[i] for r in bw]) / mean([r[i] for r in wi])))
