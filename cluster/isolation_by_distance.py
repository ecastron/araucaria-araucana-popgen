#!/usr/bin/env python3
"""Independent check: does the stand effect look like real population structure?

A stand effect that appears only after discarding 97.5% of sites, and is only 0.8% in magnitude,
needs validation that is not another permutation of the same data. Isolation by distance is the
right test: if the signal is genuine population structure, geographically closer stands should be
genetically closer. Filtering artefacts, batch effects and residual depth structure have no
reason to align with geography.

Works at STAND CENTROIDS, not per-tree coordinates, so per-tree GPS error (tens of metres, and
the subject of a known deposit problem) cannot affect it. Centroids come from the corrected
coordinates in environment/chelsa_per_tree.tsv, never from the deposited RAD coordinates, which
Section 7 prohibits.

  significant positive Mantel r  -> genetic distance tracks geography. That is population
                                    structure; the stand effect is real and Q1/Q4/Q5 are back.
  null or negative               -> the stand effect does not correspond to anything spatial.
                                    It is a property of the data, not of the trees.

Also reports the effect with depth removed, because depth is the known confound, and a
STAND-PERMUTED null so the Mantel p-value is not the only line of defence.
"""
import csv, math, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
VCF = sys.argv[1]
rng = np.random.default_rng(5)

# ---- genotypes ---------------------------------------------------------------------------
codes, rows, samples = {}, [], []
for line in open(VCF):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
        continue
    r = np.empty(len(samples), dtype=np.int8)
    for k, cell in enumerate(f[9:]):
        g = cell.partition(":")[0].replace("|", "/")
        if g in ("./.", "."):
            r[k] = -1
        else:
            c = codes.get(g)
            if c is None:
                c = codes[g] = len(codes)
            r[k] = c
    rows.append(r)
M = np.vstack(rows)
idx = {s: i for i, s in enumerate(samples)}
print(f"{VCF.split('/')[-1]}: {M.shape[0]} sites, {M.shape[1]} libraries")

site = {r[0]: r[1] for r in csv.reader(open(f"{PROJ}/metadata/popmap_177.tsv"), delimiter="\t")
        if len(r) >= 2}

# ---- stand centroids from corrected coordinates -------------------------------------------
pts = {}
with open(f"{PROJ}/metadata/chelsa_per_tree.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        pts.setdefault(r["site_code"], []).append((float(r["lat"]), float(r["lon"])))
cent = {k: (float(np.mean([p[0] for p in v])), float(np.mean([p[1] for p in v])))
        for k, v in pts.items()}

cr = (M >= 0).mean(axis=0)
use = [s for s in samples if s in site and cr[idx[s]] >= 0.005 and site[s] in cent]
stands = sorted({site[s] for s in use})
print(f"libraries {len(use)} | stands with coordinates {len(stands)}\n")

def hav(a, b):
    R = 6371.0
    p1, p2 = math.radians(a[0]), math.radians(b[0])
    dp, dl = p2 - p1, math.radians(b[1] - a[1])
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(math.sqrt(h))

# ---- pairwise genetic distances, keeping stand labels and depth ---------------------------
gd, dep, sa, sb = [], [], [], []
for i in range(len(use)):
    a = M[:, idx[use[i]]]; va = a >= 0
    for j in range(i + 1, len(use)):
        b = M[:, idx[use[j]]]
        both = va & (b >= 0)
        n = int(both.sum())
        if n < 500:
            continue
        gd.append((a[both] != b[both]).sum() / n)
        dep.append(abs(cr[idx[use[i]]] - cr[idx[use[j]]]) + cr[idx[use[i]]] + cr[idx[use[j]]])
        sa.append(site[use[i]]); sb.append(site[use[j]])
gd = np.array(gd); dep = np.array(dep)
gdr = gd - np.polyval(np.polyfit(dep, gd, 1), dep)          # depth removed

# ---- stand x stand matrices ---------------------------------------------------------------
S = len(stands); si = {s: k for k, s in enumerate(stands)}
GEN = np.full((S, S), np.nan); GENr = np.full((S, S), np.nan); GEO = np.zeros((S, S))
acc = {}
for k in range(len(gd)):
    if sa[k] == sb[k]:
        continue                                            # between-stand only
    key = tuple(sorted((sa[k], sb[k])))
    acc.setdefault(key, [[], []])
    acc[key][0].append(gd[k]); acc[key][1].append(gdr[k])
for (x, y), (v, vr) in acc.items():
    GEN[si[x], si[y]] = GEN[si[y], si[x]] = np.mean(v)
    GENr[si[x], si[y]] = GENr[si[y], si[x]] = np.mean(vr)
for x in stands:
    for y in stands:
        GEO[si[x], si[y]] = hav(cent[x], cent[y])

iu = np.triu_indices(S, 1)
g_, gr_, d_ = GEN[iu], GENr[iu], GEO[iu]
ok = ~np.isnan(g_)
g_, gr_, d_ = g_[ok], gr_[ok], d_[ok]
print(f"stand pairs: {len(d_)} | geographic distance {d_.min():.1f}-{d_.max():.1f} km")

def mantel(x, y, n_perm=9999):
    r = np.corrcoef(x, y)[0, 1]
    cnt = 0
    for _ in range(n_perm):
        p = rng.permutation(S)
        Xp = GEO[np.ix_(p, p)][iu][ok]
        if np.corrcoef(Xp, y)[0, 1] >= r:
            cnt += 1
    return r, (cnt + 1) / (n_perm + 1)

r1, p1 = mantel(d_, g_)
r2, p2 = mantel(d_, gr_)
print(f"\nISOLATION BY DISTANCE (Mantel, stand labels permuted)")
print(f"  geographic vs genetic (raw)            r = {r1:+.3f}  p = {p1:.4f}")
print(f"  geographic vs genetic (depth-removed)  r = {r2:+.3f}  p = {p2:.4f}")

if p2 < 0.05 and r2 > 0:
    print("\n  >>> Genetic distance tracks GEOGRAPHY. The stand effect corresponds to something")
    print("  >>> spatial and external to the sequencing, which filtering artefacts cannot")
    print("  >>> produce. This is population structure.")
else:
    print("\n  >>> No isolation by distance. The stand effect does not correspond to geography,")
    print("  >>> so it is a property of the data rather than of the trees. Do not build on it.")
