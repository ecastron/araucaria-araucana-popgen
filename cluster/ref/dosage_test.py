#!/usr/bin/env python3
"""Does hemizygosity explain the flat allele balance?

Design (no genotype calls anywhere; depth and raw read counts only):
  1. normalised depth per tree per single-copy locus  = locus mean depth
     divided by that tree's median across all loci.  A tree carrying one
     haplotype at a locus sits near 0.5; a tree carrying two sits near 1.0.
  2. a locus is "dosage-polymorphic" if both classes are well populated.
  3. THE TEST: within one such locus, compare the low-dosage trees to the
     high-dosage trees.  Under hemizygosity the low-dosage trees should be
     essentially never heterozygous, and the high-dosage trees should look
     normal.  If both classes are equally flat, hemizygosity is not the cause.
"""
import gzip, glob, sys
import numpy as np

R = "/home/ecastron/araucaria/results/sv/dosage"
VCF = "/home/ecastron/araucaria/results/polyrad_fb/fb_AD.vcf.gz"
LO, HI = (0.25, 0.75), (0.75, 1.60)   # normalised-depth windows
MIN_PER_CLASS = 15                     # trees needed in each class
MIN_READS = 2                          # per allele, to call a site het from AD

# ---- 1. depth matrix -------------------------------------------------------
names, rows = [], []
for f in sorted(glob.glob(f"{R}/bedcov.part*.tsv")):
    for l in open(f):
        p = l.rstrip("\n").split("\t")
        ln = int(p[2]) - int(p[1])
        if ln <= 0:
            continue
        names.append((p[0], int(p[1]), int(p[2]), p[3]))
        rows.append([int(x) / ln for x in p[4:]])
D = np.asarray(rows, dtype=np.float32)
print(f"loci {D.shape[0]} x trees {D.shape[1]}")

samples = [l.strip().split("/")[-1].replace(".bam", "")
           for l in open("/home/ecastron/araucaria/results/sv/bamlist177.txt") if l.strip()]
assert len(samples) == D.shape[1], (len(samples), D.shape[1])

# RAD data: ~63% of (locus, tree) cells have zero coverage simply because the
# locus was not sampled in that tree. Taking a plain median across ALL loci gives
# 0 for 155 of 177 trees and silently voids the whole test. Normalise on covered
# loci only, and treat zero as MISSING rather than as a dosage of zero -- a zero
# could be a homozygous deletion or ordinary RAD dropout, and depth alone cannot
# tell those apart.
COV = np.where(D > 0, D, np.nan)
med = np.nanmedian(COV, axis=0)
print("per-tree median depth over covered loci: min %.1f med %.1f max %.1f"
      % (np.nanmin(med), np.nanmedian(med), np.nanmax(med)))
N = COV / med                          # normalised depth, NaN where uncovered
ncov = np.sum(~np.isnan(N), axis=1)
print("loci with >=60 trees covered: %d / %d" % ((ncov >= 60).sum(), len(ncov)))

# ---- 2. dosage-polymorphic loci -------------------------------------------
with np.errstate(invalid='ignore'):
    lo_mask = (N >= LO[0]) & (N < LO[1])
    hi_mask = (N >= HI[0]) & (N < HI[1])
nlo, nhi = lo_mask.sum(1), hi_mask.sum(1)
cand = np.where((nlo >= MIN_PER_CLASS) & (nhi >= MIN_PER_CLASS) & (ncov >= 60))[0]
print("normalised depth over covered cells: p5 %.2f p50 %.2f p95 %.2f"
      % tuple(np.nanpercentile(N, [5, 50, 95])))
print(f"dosage-polymorphic loci: {len(cand)} / {D.shape[0]} "
      f"({100*len(cand)/D.shape[0]:.1f}%)")
if len(cand) == 0:
    print("NO dosage polymorphism detected -> hemizygosity hypothesis fails at step 1")
    sys.exit(0)

# interval lookup for the candidate loci
byctg = {}
for i in cand:
    c, s, e, nm = names[i]
    byctg.setdefault(c, []).append((s, e, i))
for c in byctg:
    byctg[c].sort()

# ---- 3. read counts inside those loci, split by dosage class --------------
# accumulators: [n_sites, n_het, sum_balance, sum_balance_sq] per class
acc = {"lo": np.zeros(4), "hi": np.zeros(4)}
vcf_samples = None
with gzip.open(VCF, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            vcf_samples = line.rstrip("\n").split("\t")[9:]
            # column order in the VCF need not match bamlist order
            idx = {s: k for k, s in enumerate(vcf_samples)}
            colmap = np.array([idx[s] for s in samples])
            continue
        f = line.rstrip("\n").split("\t")
        ivs = byctg.get(f[0])
        if not ivs:
            continue
        pos = int(f[1]) - 1
        loc = None
        for s, e, i in ivs:
            if s <= pos < e:
                loc = i
                break
            if s > pos:
                break
        if loc is None:
            continue
        fmt = f[8].split(":")
        if "AD" not in fmt:
            continue
        adi = fmt.index("AD")
        gts = f[9:]
        for cls, mask in (("lo", lo_mask), ("hi", hi_mask)):
            sel = np.where(mask[loc])[0]
            for t in sel:
                cell = gts[colmap[t]].split(":")
                if len(cell) <= adi:
                    continue
                ad = cell[adi].split(",")
                if len(ad) < 2:
                    continue
                try:
                    r, a = int(ad[0]), int(ad[1])
                except ValueError:
                    continue
                if r + a < 8:
                    continue
                acc[cls][0] += 1
                if r >= MIN_READS and a >= MIN_READS:
                    acc[cls][1] += 1
                    b = a / (r + a)
                    acc[cls][2] += b
                    acc[cls][3] += b * b

print("\n=== THE TEST: same loci, trees split by haplotype dosage ===")
print(f"{'class':<6}{'sites':>12}{'het sites':>12}{'H_obs':>9}"
      f"{'mean bal':>10}{'sd bal':>9}")
for cls in ("lo", "hi"):
    n, h, sb, sq = acc[cls]
    if n == 0:
        print(f"{cls:<6}{'no data':>12}")
        continue
    hob = h / n
    mb = sb / h if h else float("nan")
    sd = np.sqrt(max(sq / h - mb * mb, 0)) if h else float("nan")
    lab = "low (~0.5x)" if cls == "lo" else "high (~1.0x)"
    print(f"{lab:<14}{int(n):>10}{int(h):>12}{hob:>9.4f}{mb:>10.4f}{sd:>9.4f}")

n_lo, h_lo = acc["lo"][0], acc["lo"][1]
n_hi, h_hi = acc["hi"][0], acc["hi"][1]
if n_lo and n_hi:
    ratio = (h_hi / n_hi) / (h_lo / n_lo) if h_lo else float("inf")
    print(f"\nH_obs ratio high/low = {ratio:.3f}")
    print("Hemizygosity predicts this to be LARGE (low-dosage trees cannot be "
          "heterozygous).\nA ratio near 1.0 means dosage does not explain the "
          "heterozygote deficit.")
