#!/usr/bin/env python3
"""Does the simulated error model reproduce the error we actually observe?

The relatedness simulation passed both criteria decisively (Ritland recovered 0.00/0.25/0.50
almost exactly, r = 0.996). That result is only meaningful if the simulated error process
matches the real one. It probably does not:

  - real replicate pairs (true relatedness 1.0) disagree at 9.6% of sites
  - the binomial dropout model predicts near-zero miscalls at the median depth of 19
    (P(het miscalled) ~ 2 * 0.5^19)

Those cannot both be true, so the simulation is likely far more optimistic than reality and the
PASS is not yet trustworthy.

THE CALIBRATION TEST: simulate two independent observations of the SAME individual and measure
their mismatch. That is exactly what a replicate pair is, so it must reproduce ~0.096.

  simulated mismatch ~ 0.096  -> model calibrated; the relatedness PASS stands
  simulated mismatch << 0.096 -> model too optimistic; the PASS is an artefact and the real
                                 error mechanism is something the model omits

Also reports observed heterozygosity by depth for the simulation, to compare against the real
curve (0.0139 at 5-7x rising to 0.0952 at 35+x). A model that cannot reproduce that curve is not
describing this dataset.

Candidate omitted mechanism: the assembly ran with max_Hs_consens = 0.05, a hard cap on
heterozygous sites per consensus read. That is a FILTER, not a sampling process, and no binomial
read model reproduces it.
"""
import re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1] if len(sys.argv) > 1 else (
    f"{PROJ}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")
NSITES = 20000
rng = np.random.default_rng(11)

freqs, depths, n_read = [], [], 0
for line in open(V):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    ref = alt = 0
    for c in f[9:]:
        p = c.split(":")
        g = p[0].replace("|", "/")
        try:
            d = int(p[1])
        except (IndexError, ValueError):
            d = 0
        if d > 0:
            depths.append(d)
        if g in ("./.", ".") or d < 20:
            continue
        a, _, b = g.partition("/")
        for x in (a, b):
            if x == "0":
                ref += 1
            else:
                alt += 1
    if ref + alt >= 20:
        freqs.append(alt / (ref + alt))
    n_read += 1
    if n_read >= 40000:
        break

freqs = np.array([q for q in freqs if 0.05 < q < 0.95])
depths = np.array(depths)
p = rng.choice(freqs, NSITES)
print(f"depth distribution: median {np.median(depths):.0f}, 10th {np.percentile(depths,10):.0f}, "
      f"90th {np.percentile(depths,90):.0f}")
print(f"simulating {NSITES} loci\n")

def observe(g):
    d = rng.choice(depths, len(g))
    out = g.astype(float).copy()
    het = g == 1
    if het.any():
        idx = np.where(het)[0]
        k = rng.binomial(d[het], 0.5)
        out[idx[k == 0]] = 0
        out[idx[k == d[het]]] = 2
    out[d < 5] = np.nan
    return out, d

# One individual, observed twice -- the definition of a replicate pair.
truth = rng.binomial(2, p)
o1, d1 = observe(truth)
o2, d2 = observe(truth)
ok = ~np.isnan(o1) & ~np.isnan(o2)
mism = np.mean(o1[ok] != o2[ok])

hh = np.mean(((o1[ok] == 1) & (o2[ok] != 1)) | ((o2[ok] == 1) & (o1[ok] != 1))
             ) if (o1[ok] != o2[ok]).any() else 0
diff = o1[ok] != o2[ok]
if diff.sum():
    a, b = o1[ok][diff], o2[ok][diff]
    het_hom = np.mean((a == 1) | (b == 1))
    opp = np.mean(((a == 0) & (b == 2)) | ((a == 2) & (b == 0)))
else:
    het_hom = opp = 0.0

print("CALIBRATION -- simulated replicate pair (same individual, observed twice)")
print(f"  simulated mismatch : {mism:.4f}")
print(f"  OBSERVED mismatch  : 0.0959   (10 real validated replicate pairs)")
print(f"  simulated composition: het-vs-hom {100*het_hom:.1f}%, opposite-hom {100*opp:.1f}%")
print(f"  OBSERVED composition : het-vs-hom      63.3%, opposite-hom      36.7%")

print("\nsimulated heterozygosity by depth")
print(f"  {'depth':>10}{'genotypes':>12}{'het rate':>10}")
bins = [(5, 7), (8, 11), (12, 19), (20, 34), (35, 10**9)]
for lo, hi in bins:
    m = (d1 >= lo) & (d1 <= hi) & ~np.isnan(o1)
    if m.sum():
        lab = f"{lo}-{hi}" if hi < 10**9 else f"{lo}+"
        print(f"  {lab:>10}{int(m.sum()):>12}{np.mean(o1[m] == 1):>10.4f}")
print("  REAL: 0.0139 / 0.0126 / 0.0206 / 0.0385 / 0.0952 across the same bins")

ratio = mism / 0.0959 if mism else 0.0
print(f"\nsimulated / observed mismatch = {ratio:.3f}")
if ratio < 0.5:
    print("\n  >>> MODEL NOT CALIBRATED. It reproduces far less error than the data contain, so")
    print("  >>> the relatedness PASS is optimistic and must not be acted on. The dominant real")
    print("  >>> error mechanism is missing from the model.")
elif ratio > 2:
    print("\n  >>> Model overshoots the observed error; the PASS is conservative.")
else:
    print("\n  >>> Model reproduces the observed error rate. The relatedness PASS stands.")
