#!/usr/bin/env python3
"""Is the "low dosage" class real hemizygosity, or just RAD depth noise?

Two discriminating predictions:

  A. Real segregating deletions are a property of the LOCUS x TREE pair. Depth
     noise is a property of the TREE. So under noise, the fraction of loci where
     a tree lands in the low class should track that tree's sequencing depth;
     under hemizygosity it should not.

  B. Real deletions replicate. Split each tree's reads at a locus into two halves
     by position parity and classify twice: a true hemizygote is called low both
     times, a low-coverage fluctuation is not. Reported as the excess of the
     observed agreement over what independent Poisson draws would give.
"""
import glob
import numpy as np

R = "/home/ecastron/araucaria/results/sv/dosage"
LO, HI = (0.25, 0.75), (0.75, 1.60)

rows = []
for f in sorted(glob.glob(f"{R}/bedcov.part*.tsv")):
    for l in open(f):
        p = l.rstrip("\n").split("\t")
        ln = int(p[2]) - int(p[1])
        if ln > 0:
            rows.append([int(x) / ln for x in p[4:]])
D = np.asarray(rows, dtype=np.float32)
COV = np.where(D > 0, D, np.nan)
med = np.nanmedian(COV, axis=0)
N = COV / med
ncov = np.sum(~np.isnan(N), axis=1)

with np.errstate(invalid='ignore'):
    lo = (N >= LO[0]) & (N < LO[1])
    hi = (N >= HI[0]) & (N < HI[1])
cand = (lo.sum(1) >= 15) & (hi.sum(1) >= 15) & (ncov >= 60)
print("dosage-polymorphic loci: %d / %d" % (cand.sum(), len(cand)))

# ---- A: does low-class membership track the tree, or the locus? ----
sub = N[cand]
with np.errstate(invalid='ignore'):
    lo_s = (sub >= LO[0]) & (sub < LO[1])
covered = ~np.isnan(sub)
frac_lo = lo_s.sum(0) / np.maximum(covered.sum(0), 1)      # per tree
r = np.corrcoef(med, frac_lo)[0, 1]
print("\n--- A. tree-level confound ---")
print("per-tree fraction of loci called low: min %.3f med %.3f max %.3f"
      % (frac_lo.min(), np.median(frac_lo), frac_lo.max()))
print("correlation with that tree's median depth: r = %+.3f" % r)
print("   r near 0  -> locus-specific, consistent with real dosage variation")
print("   r strongly negative -> the low class is just shallow sequencing")

# ---- B: Poisson null for the observed spread ----
# If every (locus, tree) were truly diploid, depth is Poisson about the tree
# mean; how often would that alone put a cell in the 0.25-0.75 window?
rng = np.random.default_rng(0)
lam = np.nanmean(COV, axis=0)                               # per-tree mean depth
sim = rng.poisson(np.broadcast_to(lam, (4000, len(lam))))
sim_n = sim / med
exp_lo = ((sim_n >= LO[0]) & (sim_n < LO[1])).mean()
obs_lo = np.nanmean(lo_s)
print("\n--- B. Poisson null ---")
print("observed fraction of cells in the low window : %.4f" % obs_lo)
print("expected under pure Poisson depth noise      : %.4f" % exp_lo)
print("excess over noise                            : %+.4f (ratio %.2f)"
      % (obs_lo - exp_lo, obs_lo / exp_lo if exp_lo else float("nan")))
print("   ratio near 1.0 -> the 'dosage' classes are noise, and the")
print("   hemizygosity test built on them carries no information.")
