#!/usr/bin/env python3
"""Does Ts/Tv recover under filtering, or are these SNPs random throughout?

Across all four clustering thresholds the transition/transversion ratio is 0.47-0.50. For any
reference base there is ONE possible transition and TWO possible transversions, so 0.5 is exactly
the ratio produced by random base changes. Real plant variation runs 1.5-2.5. As a set, these
SNPs are indistinguishable from noise.

That is Shafer et al.'s null expectation failing, and it is independent of the replicate and
Wahlund questions. But it sits against a real isolation-by-distance signal (Mantel r = +0.601),
which cannot arise from noise. Both are consistent only if genuine variation is present and
diluted by a large excess of spurious SNPs.

TEST: recompute Ts/Tv under increasing stringency -- call rate, per-genotype depth, and minor
allele frequency.

  Ts/Tv rises toward ~2 with filtering -> real variation is there and recoverable; the filters
                                          define the usable subset, and the r80/threshold choice
                                          should be made within it.
  Ts/Tv stays ~0.5 everywhere          -> the SNP set is dominated by artefacts at every level
                                          of stringency, and no filtering rescues it.

MAF matters here: singletons and very rare variants are where sequencing error concentrates, so
a Ts/Tv that recovers only after a MAF floor tells us the error is in rare calls rather than
spread across all loci.
"""
import re, sys
import numpy as np

V = sys.argv[1]
TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

rows = []
samples = []
for line in open(V):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples = f[9:]
        continue
    ref, alt = f[3], f[4].split(",")[0]
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        continue
    n = 0; alt_ct = 0; tot_al = 0; deep = 0
    for c in f[9:]:
        p = c.split(":")
        g = p[0].replace("|", "/")
        try:
            d = int(p[1])
        except (IndexError, ValueError):
            d = 0
        if g in ("./.", "."):
            continue
        n += 1
        if d >= 20:
            deep += 1
        a, _, b = g.partition("/")
        alt_ct += (a != "0") + (b != "0")
        tot_al += 2
    if n == 0:
        continue
    maf = min(alt_ct, tot_al - alt_ct) / tot_al if tot_al else 0
    rows.append((1 if (ref, alt) in TS else 0, n / len(samples), deep / len(samples), maf))

R = np.array(rows)
print(f"{V.split('/')[-1]}: {len(R)} biallelic SNPs, {len(samples)} libraries")
print("expected Ts/Tv: ~2 for real plant variation; 0.5 = random substitutions\n")

def tstv(mask, label):
    sub = R[mask]
    if len(sub) < 50:
        print(f"  {label:<34}{len(sub):>9}      (too few)")
        return
    ts = sub[:, 0].sum(); tv = len(sub) - ts
    r = ts / tv if tv else float("nan")
    flag = "  <-- real variation" if r > 1.2 else ("  <-- random" if r < 0.7 else "")
    print(f"  {label:<34}{len(sub):>9}{r:>9.3f}{flag}")

print(f"  {'filter':<34}{'SNPs':>9}{'Ts/Tv':>9}")
tstv(np.ones(len(R), bool), "none")
for c in (0.3, 0.5, 0.7):
    tstv(R[:, 1] >= c, f"call rate >= {c:.0%}")
for c in (0.3, 0.5):
    tstv(R[:, 2] >= c, f">= {c:.0%} of libs at >=20x")
for m in (0.02, 0.05, 0.10, 0.20):
    tstv(R[:, 3] >= m, f"MAF >= {m:.2f}")
print()
tstv((R[:, 1] >= 0.5) & (R[:, 3] >= 0.05), "call rate >=50% AND MAF >=0.05")
tstv((R[:, 2] >= 0.3) & (R[:, 3] >= 0.05), ">=30% deep AND MAF >=0.05")
tstv((R[:, 1] >= 0.5) & (R[:, 3] >= 0.10), "call rate >=50% AND MAF >=0.10")
