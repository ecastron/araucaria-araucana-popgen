#!/usr/bin/env python3
"""Can relatedness be estimated on THIS data, given its dropout regime?

Attard et al. (2018) recommend simulating individuals of KNOWN relatedness under a dataset's own
error regime and choosing the estimator empirically, rather than assuming any estimator is safe.
That is what this does, using the dropout process measured from our own replicate pairs.

WHY THIS TEST AND NOT AN ARGUMENT: the diagnosis of this dataset has been revised twice. Both
revisions came from testing a premise that had been asserted. So nothing here rests on dropout
being "correctable in principle" -- the question is whether kinship classes are separable in
simulation at our measured dropout rate, which is a fact about our data, not about methods.

PASS/FAIL CRITERIA, FIXED BEFORE RUNNING:

  Q2 (SGS, sigma) needs the kinship-DISTANCE SLOPE, so it needs relatedness to correlate with
  truth and to separate related from unrelated:
      PASS if Pearson r(estimated, true) >= 0.7 AND the unrelated (0.0) and half-sib (0.25)
      distributions separate with < 10% overlap.

  Q3 (dispersal guild) needs the SHAPE of the kinship-distance curve, so it needs the kinship
  CLASSES resolved, not just ordered:
      PASS if half-sib (0.25) and full-sib (0.5) separate with < 10% overlap.
      Q3 is the stricter test and may fail while Q2 passes.

  Q6 (individual extended phenotype) needs per-pair relatedness within stands:
      PASS on the same criterion as Q2, since it uses relatedness as a continuous predictor.

Dropout model: true genotypes are drawn under the stated kinship, then each individual's
genotype at each site is observed at a depth drawn from the EMPIRICAL depth distribution, with
reads sampled binomially. A heterozygote whose reads all come from one allele is called
homozygous -- the actual mechanism, not an abstract error rate.
"""
import csv, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1] if len(sys.argv) > 1 else (
    f"{PROJ}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")
NSITES = 5000          # simulated loci; well above the >7,000-locus guidance once scaled
NPAIRS = 300           # pairs per kinship class
rng = np.random.default_rng(11)

# ---- learn allele frequencies and the depth distribution from the real data ---------------
# Allele frequencies are taken from HIGH-DEPTH genotypes only. Using all genotypes would embed
# the dropout bias into the frequencies themselves and make the simulation flatter than reality.
freqs, depths = [], []
n_read = 0
for line in open(V):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    alt = ref = 0
    for c in f[9:]:
        p = c.split(":")
        g = p[0].replace("|", "/")
        try:
            d = int(p[1])
        except (IndexError, ValueError):
            d = 0
        if d > 0:
            depths.append(d)
        if g in ("./.", ".") or d < 20:          # high-depth only for frequency estimation
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

freqs = np.array([p for p in freqs if 0.05 < p < 0.95])
depths = np.array(depths)
print(f"learned from {V.split('/')[-1]}")
print(f"  {len(freqs)} polymorphic sites (MAF>0.05) for allele frequencies")
print(f"  depth distribution: median {np.median(depths):.0f}, "
      f"mean {depths.mean():.1f}, 10th {np.percentile(depths,10):.0f}, "
      f"90th {np.percentile(depths,90):.0f}\n")

p = rng.choice(freqs, NSITES)

def founder(p):
    """Draw allele counts (0/1/2 copies of ALT) under HWE."""
    return rng.binomial(2, p)

def offspring_alleles(gA, gB):
    """One allele from each parent genotype (allele counts 0-2)."""
    a = (rng.random(len(gA)) < gA / 2).astype(int)
    b = (rng.random(len(gB)) < gB / 2).astype(int)
    return a + b

def make_pair(kind, p):
    """Return two genotype vectors with the stated true relationship."""
    if kind == "unrelated":
        return founder(p), founder(p)
    if kind == "half-sib":
        m = founder(p); f1 = founder(p); f2 = founder(p)
        return offspring_alleles(m, f1), offspring_alleles(m, f2)
    if kind == "full-sib":
        m = founder(p); f = founder(p)
        return offspring_alleles(m, f), offspring_alleles(m, f)
    if kind == "parent-offspring":
        m = founder(p); f = founder(p)
        return m, offspring_alleles(m, f)
    raise ValueError(kind)

def observe(g, p):
    """Apply the measured dropout process: sample reads binomially at empirical depth."""
    d = rng.choice(depths, len(g))
    out = g.copy()
    het = g == 1
    if het.any():
        k = rng.binomial(d[het], 0.5)              # ALT reads seen out of d
        allref = k == 0
        allalt = k == d[het]
        idx = np.where(het)[0]
        out[idx[allref]] = 0                       # dropped to hom-ref
        out[idx[allalt]] = 2                       # dropped to hom-alt
    miss = d < 5                                   # too shallow to call at all
    out = out.astype(float)
    out[miss] = np.nan
    return out

# ---- estimators ---------------------------------------------------------------------------
def ritland(x, y, p):
    """Ritland (1996) kinship. Attard et al. found it least affected by allelic dropout."""
    ok = ~np.isnan(x) & ~np.isnan(y)
    if ok.sum() < 100:
        return np.nan
    xi, yi, pi = x[ok] / 2, y[ok] / 2, p[ok]
    num = ((xi - pi) * (yi - pi)) / (pi * (1 - pi))
    return 2 * np.nanmean(num)

def ibs_sim(x, y, p):
    """Naive identity-by-state similarity, as a baseline that ignores dropout."""
    ok = ~np.isnan(x) & ~np.isnan(y)
    if ok.sum() < 100:
        return np.nan
    return 1 - np.mean(np.abs(x[ok] - y[ok])) / 2

TRUE = {"unrelated": 0.0, "half-sib": 0.25, "full-sib": 0.5, "parent-offspring": 0.5}
res = {k: {"ritland": [], "ibs": []} for k in TRUE}
for kind in TRUE:
    for _ in range(NPAIRS):
        gA, gB = make_pair(kind, p)
        oA, oB = observe(gA, p), observe(gB, p)
        res[kind]["ritland"].append(ritland(oA, oB, p))
        res[kind]["ibs"].append(ibs_sim(oA, oB, p))

print("estimated relatedness by TRUE relationship "
      f"({NPAIRS} pairs each, {NSITES} loci, empirical dropout)")
print(f"  {'relationship':<18}{'true':>6}{'Ritland mean':>14}{'sd':>8}{'IBS mean':>10}{'sd':>8}")
for kind in ("unrelated", "half-sib", "full-sib", "parent-offspring"):
    r = np.array(res[kind]["ritland"], float); i = np.array(res[kind]["ibs"], float)
    print(f"  {kind:<18}{TRUE[kind]:>6.2f}{np.nanmean(r):>14.4f}{np.nanstd(r):>8.4f}"
          f"{np.nanmean(i):>10.4f}{np.nanstd(i):>8.4f}")

def overlap(a, b):
    """Fraction misclassified by the midpoint threshold between two classes."""
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    t = (np.mean(a) + np.mean(b)) / 2
    if np.mean(a) < np.mean(b):
        return (np.mean(a > t) + np.mean(b < t)) / 2
    return (np.mean(a < t) + np.mean(b > t)) / 2

R = {k: np.array(res[k]["ritland"], float) for k in TRUE}
truth, est = [], []
for k in TRUE:
    v = R[k][~np.isnan(R[k])]
    truth += [TRUE[k]] * len(v); est += list(v)
corr = np.corrcoef(truth, est)[0, 1]

o_un_hs = overlap(R["unrelated"], R["half-sib"])
o_hs_fs = overlap(R["half-sib"], R["full-sib"])
print(f"\nPearson r(estimated, true) = {corr:.3f}")
print(f"overlap unrelated vs half-sib : {100*o_un_hs:.1f}%")
print(f"overlap half-sib vs full-sib  : {100*o_hs_fs:.1f}%")

print("\nVERDICT against pre-declared criteria")
q2 = corr >= 0.7 and o_un_hs < 0.10
q3 = o_hs_fs < 0.10
print(f"  Q2 / Q6 (slope, continuous relatedness): {'PASS' if q2 else 'FAIL'}"
      f"  [needs r>=0.7 and unrelated-vs-halfsib overlap <10%]")
print(f"  Q3 (kinship classes resolved):           {'PASS' if q3 else 'FAIL'}"
      f"  [needs halfsib-vs-fullsib overlap <10%]")
