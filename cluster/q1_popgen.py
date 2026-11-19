#!/usr/bin/env python3
"""Q1: how is genetic variation partitioned within versus among stands?

Analysis plan is Step 3 of the spec. Working dataset: ipyrad ct = 0.91 at a 30% call-rate floor
(78,802 SNPs, 177 libraries). That floor is validated by Ts/Tv -- 0.48 unfiltered (the ratio
random base changes produce), 1.522 filtered (real plant variation).

Outputs, each written to its own file so a later step can fail without losing earlier work:

  q1_per_stand_diversity.tsv   Ho, He, FIS per stand, computed at >=20x
  q1_fst_matrix.tsv            pairwise Weir & Cockerham FST between stands
  q1_pca.tsv                   PCA coordinates per library, for structure plots
  q1_genetic_distance.tsv      library x library distance, input to AMOVA in R
  q1_ne_surface.tsv            Ne sensitivity across MAF x missingness (see caveat below)
  q1_summary.txt               everything in one readable place

TWO SPEC RULES ENFORCED HERE:

  FIS is reported WITHIN stands, never pooled. Pooling the twelve stands inflates FIS by +0.474
  through the Wahlund effect, which is a property of the sampling design, not of the data.

  Ne is reported as a SENSITIVITY SURFACE across MAF and missingness, not a point estimate.
  With 10-18 trees per stand and depth-limited genotypes, a single value would be an artefact of
  the filtering choice. The surface makes that dependence visible instead of hiding it.
"""
import csv, os, re, sys, itertools
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1] if len(sys.argv) > 1 else (
    f"{PROJ}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")
OUT = f"{PROJ}/results/Q1"
os.makedirs(OUT, exist_ok=True)
MINDP = 20

log = open(f"{OUT}/q1_summary.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); log.write(s + "\n"); log.flush()

say(f"Q1 population genomics\ninput: {V}\n")

samples, G, D = [], [], []
for line in open(V):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
        continue
    g = np.full(len(samples), np.nan, dtype=np.float32)
    d = np.zeros(len(samples), dtype=np.int32)
    for k, c in enumerate(f[9:]):
        p = c.split(":")
        gt = p[0].replace("|", "/")
        try:
            d[k] = int(p[1])
        except (IndexError, ValueError):
            d[k] = 0
        if gt in ("./.", "."):
            continue
        a, _, b = gt.partition("/")
        g[k] = (a != "0") + (b != "0")
    G.append(g); D.append(d)
G = np.vstack(G); D = np.vstack(D)
idx = {s: i for i, s in enumerate(samples)}
say(f"{G.shape[0]} SNPs x {G.shape[1]} libraries")

site_of = {r[0]: r[1] for r in csv.reader(open(f"{PROJ}/metadata/popmap_177.tsv"), delimiter="\t")
           if len(r) >= 2}
use = [s for s in samples if s in site_of]
stands = sorted({site_of[s] for s in use})
cols = {st: np.array([idx[s] for s in use if site_of[s] == st]) for st in stands}
say(f"{len(use)} libraries across {len(stands)} stands: " +
    ", ".join(f"{st}={len(cols[st])}" for st in stands) + "\n")

# ---- 1. per-stand diversity, at >=20x -----------------------------------------------------
hi = np.where(D >= MINDP, G, np.nan)
say(f"--- per-stand diversity (genotypes at >={MINDP}x) ---")
rows = []
for st in stands:
    sub = hi[:, cols[st]]
    n = np.sum(~np.isnan(sub), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.nansum(sub, axis=1) / (2 * n)
        ho = np.nansum(sub == 1, axis=1) / n
        he = 2 * p * (1 - p)
        ok = (n >= 6) & (he > 0.05)
        fis = 1 - ho[ok] / he[ok]
    rows.append((st, len(cols[st]), int(ok.sum()), float(np.mean(ho[ok])),
                 float(np.mean(he[ok])), float(np.mean(fis))))
    say(f"  {st:<6} n={len(cols[st]):>3}  loci={int(ok.sum()):>5}  "
        f"Ho={np.mean(ho[ok]):.4f}  He={np.mean(he[ok]):.4f}  FIS={np.mean(fis):+.4f}")
with open(f"{OUT}/q1_per_stand_diversity.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["stand", "n_libraries", "n_loci", "Ho", "He", "FIS"])
    w.writerows(rows)
say(f"  mean within-stand FIS = {np.mean([r[5] for r in rows]):+.4f}")
say("  (pooled FIS is NOT reported: Wahlund inflates it by ~+0.47)\n")

# ---- 2. pairwise FST, Weir & Cockerham ----------------------------------------------------
def wc_fst(a_cols, b_cols):
    A, B = hi[:, a_cols], hi[:, b_cols]
    na = np.sum(~np.isnan(A), axis=1); nb = np.sum(~np.isnan(B), axis=1)
    ok = (na >= 4) & (nb >= 4)
    if ok.sum() < 100:
        return np.nan, int(ok.sum())
    A, B, na, nb = A[ok], B[ok], na[ok], nb[ok]
    pa = np.nansum(A, axis=1) / (2 * na)
    pb = np.nansum(B, axis=1) / (2 * nb)
    pbar = (na * pa + nb * pb) / (na + nb)
    nbar = (na + nb) / 2
    s2 = (na * (pa - pbar) ** 2 + nb * (pb - pbar) ** 2) / ((2 - 1) * nbar)
    var = pbar * (1 - pbar)
    keep = var > 0
    num = np.sum(s2[keep] - (1 / (nbar[keep] - 1)) * (var[keep] - 0.5 * s2[keep]))
    den = np.sum(var[keep] + s2[keep] / 2)
    return (num / den if den else np.nan), int(ok.sum())

say("--- pairwise FST (Weir & Cockerham) ---")
F = np.full((len(stands), len(stands)), np.nan)
for i, j in itertools.combinations(range(len(stands)), 2):
    f_, nl = wc_fst(cols[stands[i]], cols[stands[j]])
    F[i, j] = F[j, i] = f_
np.fill_diagonal(F, 0.0)
with open(f"{OUT}/q1_fst_matrix.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow([""] + stands)
    for i, st in enumerate(stands):
        w.writerow([st] + [f"{x:.5f}" if not np.isnan(x) else "NA" for x in F[i]])
iu = np.triu_indices(len(stands), 1)
say(f"  mean FST = {np.nanmean(F[iu]):.4f}, range {np.nanmin(F[iu]):.4f}-{np.nanmax(F[iu]):.4f}\n")

# ---- 3. PCA -------------------------------------------------------------------------------
say("--- PCA ---")
M = hi[:, [idx[s] for s in use]]
keep = np.sum(~np.isnan(M), axis=1) >= 0.5 * len(use)
M = M[keep]
mu = np.nanmean(M, axis=1, keepdims=True)
M = np.where(np.isnan(M), mu, M) - mu
sd = M.std(axis=1, keepdims=True); sd[sd == 0] = 1
M = M / sd
C = (M.T @ M) / M.shape[0]
val, vec = np.linalg.eigh(C)
o = np.argsort(val)[::-1]
val, vec = val[o], vec[:, o]
pv = 100 * val / val.sum()
say(f"  {M.shape[0]} loci; PC1 {pv[0]:.1f}%, PC2 {pv[1]:.1f}%, PC3 {pv[2]:.1f}%")
with open(f"{OUT}/q1_pca.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["library", "stand"] + [f"PC{i+1}" for i in range(6)])
    for k, s in enumerate(use):
        w.writerow([s, site_of[s]] + [f"{vec[k,i]:.6f}" for i in range(6)])

# ---- 4. genetic distance matrix, for AMOVA in R -------------------------------------------
say("\n--- genetic distance matrix (for AMOVA) ---")
n = len(use); Dm = np.zeros((n, n))
Msub = hi[:, [idx[s] for s in use]]
for i in range(n):
    a = Msub[:, i]; va = ~np.isnan(a)
    for j in range(i + 1, n):
        b = Msub[:, j]; both = va & ~np.isnan(b)
        Dm[i, j] = Dm[j, i] = (np.mean(np.abs(a[both] - b[both])) / 2
                               if both.sum() >= 200 else np.nan)
with open(f"{OUT}/q1_genetic_distance.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow([""] + use)
    for i, s in enumerate(use):
        w.writerow([s] + [f"{x:.6f}" if not np.isnan(x) else "NA" for x in Dm[i]])
say(f"  {n} x {n}, {100*np.mean(np.isnan(Dm)):.1f}% undefined pairs")
with open(f"{OUT}/q1_stands.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["library", "stand"])
    for s in use:
        w.writerow([s, site_of[s]])

# ---- 5. Ne sensitivity surface ------------------------------------------------------------
# Reported as a surface, never a point estimate. LD-based Ne assumes reliable genotypes; ours
# are depth-limited and each stand holds 10-18 trees, so the value moves with the filter. The
# surface makes that dependence explicit rather than burying it in a confidence interval.
say("\n--- Ne sensitivity surface (r^2-based, per stand) ---")
def ne_ld(sub, maf_min, callrate_min):
    n = np.sum(~np.isnan(sub), axis=1)
    with np.errstate(invalid="ignore"):
        p = np.nansum(sub, axis=1) / (2 * n)
    maf = np.minimum(p, 1 - p)
    ok = (maf >= maf_min) & (n >= callrate_min * sub.shape[1])
    X = sub[ok]
    if X.shape[0] < 50:
        return np.nan, int(X.shape[0])
    if X.shape[0] > 800:
        X = X[np.random.default_rng(0).choice(X.shape[0], 800, replace=False)]
    Xc = np.where(np.isnan(X), np.nanmean(X, axis=1, keepdims=True), X)
    Xc = Xc - Xc.mean(axis=1, keepdims=True)
    sd = Xc.std(axis=1, keepdims=True); sd[sd == 0] = 1
    Xc /= sd
    R = (Xc @ Xc.T) / Xc.shape[1]
    iu2 = np.triu_indices(R.shape[0], 1)
    r2 = np.mean(R[iu2] ** 2)
    S = sub.shape[1]
    r2 = r2 - 1.0 / S                      # sample-size correction
    return (1 / (3 * r2) if r2 > 0 else np.nan), int(X.shape[0])

with open(f"{OUT}/q1_ne_surface.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["stand", "maf_min", "callrate_min", "n_loci", "Ne"])
    for st in stands:
        line = []
        for maf in (0.05, 0.10, 0.20):
            for cr in (0.5, 0.8):
                ne, nl = ne_ld(hi[:, cols[st]], maf, cr)
                w.writerow([st, maf, cr, nl, f"{ne:.1f}" if np.isfinite(ne) else "NA"])
                line.append(f"maf{maf}/cr{cr}:{ne:.0f}" if np.isfinite(ne) else
                            f"maf{maf}/cr{cr}:NA")
        say(f"  {st:<6} " + "  ".join(line))
say("\n  Ne varies with filtering by design; report the surface, not a point estimate.")
say(f"\nDONE. Outputs in {OUT}/")
log.close()
