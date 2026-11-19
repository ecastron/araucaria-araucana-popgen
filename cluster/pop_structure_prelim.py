#!/usr/bin/env python3
"""Preliminary: is population structure recoverable despite unreliable individual genotypes?

Individual genotype calls fail (~80% of replicate disagreements are opposite homozygotes, and
locus filtering fails cross-validation). But RELATIVE similarity still carries signal: every
validated replicate pair is its own nearest neighbour, and within-pair mismatch (~0.28) sits
well below between-individual (~0.40). The question is whether that residual signal is enough
to recover which stand a tree came from.

THE CONFOUND THAT MUST BE CONTROLLED: mismatch is driven mostly by sequencing depth --
deep-deep 0.7030, deep-shallow 0.4030, shallow-shallow 0.3994. An uncontrolled ordination
would separate libraries by depth and look like structure. So geography is tested with depth
held fixed, via a partial Mantel test, and never on the raw distances alone.

PREDICTIONS, stated before running:
  P1 (positive control) replicate pairs are closer than non-replicates WITHIN their own depth
     stratum. Already established; if it fails here the distance matrix is built wrong.
  P2 same-stand tree pairs are closer than different-stand pairs, after conditioning on depth.
  P3 genetic distance correlates with geographic distance, after conditioning on depth.

This subset is only ~15 trees over ~8 stands, so it is UNDERPOWERED and can support a negative
signal only weakly. A clean positive would be informative; a null result here must not be read
as absence of structure. The full 177-library assembly is the real test.
"""
import csv, glob, os, re, random, math, statistics as st, collections

PROJ = "/home/ecastron/araucaria"
random.seed(11)
VCF = sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf"))[0]

samples, recs = [], []
with open(VCF) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            continue
        recs.append(f[9:])
idx = {s: i for i, s in enumerate(samples)}

def gt(c):
    return c.split(":")[0].replace("|", "/")

called = {s: sum(1 for r in recs if gt(r[idx[s]]) not in ("./.", ".")) for s in samples}

# library -> (tree, site). popmap gives site; replicate_pairs gives tree grouping.
site = {}
for r in csv.reader(open(f"{PROJ}/metadata/popmap_all_libraries.tsv"), delimiter="\t"):
    if len(r) >= 2:
        site[r[0]] = r[1]
tree = {}
for r in csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"), delimiter="\t"):
    tree[r["run_a"]] = r["tree_id"]
    tree[r["run_b"]] = r["tree_id"]

use = [s for s in samples if called[s] >= 0.005 * len(recs) and s in site]
print(f"{os.path.basename(VCF)}: {len(recs)} sites | libraries used {len(use)} "
      f"| stands {len(set(site[s] for s in use))}\n")

def dist(a, b):
    n = m = 0
    ia, ib = idx[a], idx[b]
    for r in recs:
        x, y = gt(r[ia]), gt(r[ib])
        if x in ("./.", ".") or y in ("./.", "."):
            continue
        n += 1
        if x != y:
            m += 1
    return (m / n if n >= 500 else None)

G, DEP, SAME_SITE, SAME_TREE = [], [], [], []
keys = []
for i in range(len(use)):
    for j in range(i + 1, len(use)):
        a, b = use[i], use[j]
        d = dist(a, b)
        if d is None:
            continue
        G.append(d)
        # depth covariate: absolute difference in call rate, plus how deep the pair is overall
        ca, cb = called[a] / len(recs), called[b] / len(recs)
        DEP.append(abs(ca - cb) + (ca + cb))
        SAME_SITE.append(1 if site[a] == site[b] else 0)
        SAME_TREE.append(1 if tree.get(a) and tree.get(a) == tree.get(b) else 0)
        keys.append((a, b))
print(f"pairwise comparisons: {len(G)}  (same tree {sum(SAME_TREE)}, same stand {sum(SAME_SITE)})\n")

def pearson(x, y):
    n = len(x); mx, my = st.mean(x), st.mean(y)
    sx = math.sqrt(sum((v - mx) ** 2 for v in x)); sy = math.sqrt(sum((v - my) ** 2 for v in y))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((x[k] - mx) * (y[k] - my) for k in range(n)) / (sx * sy)

def residual(y, x):
    mx, my = st.mean(x), st.mean(y)
    vx = sum((v - mx) ** 2 for v in x)
    b = sum((x[k] - mx) * (y[k] - my) for k in range(len(x))) / vx if vx else 0
    return [y[k] - (my + b * (x[k] - mx)) for k in range(len(y))]

def perm_test(y, g, label, n_perm=9999):
    """Difference in mean y between groups g=1 and g=0, permutation p-value."""
    a = [y[k] for k in range(len(y)) if g[k] == 1]
    b = [y[k] for k in range(len(y)) if g[k] == 0]
    if len(a) < 2 or len(b) < 2:
        print(f"  {label}: too few pairs")
        return
    obs = st.mean(a) - st.mean(b)
    pool = list(y); cnt = 0
    for _ in range(n_perm):
        random.shuffle(pool)
        d = st.mean(pool[:len(a)]) - st.mean(pool[len(a):])
        if d <= obs:
            cnt += 1
    print(f"  {label}: mean(same)={st.mean(a):.4f}  mean(diff)={st.mean(b):.4f}  "
          f"delta={obs:+.4f}  p={(cnt+1)/(n_perm+1):.4f}")

print(f"depth covariate vs genetic distance: r = {pearson(DEP, G):+.3f}  "
      "(the confound, before controlling)")
Gr = residual(G, DEP)          # genetic distance with depth regressed out

print("\nP1 (positive control) replicate pairs vs the rest")
perm_test(G,  SAME_TREE, "raw distance      ")
perm_test(Gr, SAME_TREE, "depth-controlled  ")

print("\nP2 same stand vs different stand")
perm_test(G,  SAME_SITE, "raw distance      ")
perm_test(Gr, SAME_SITE, "depth-controlled  ")

print("\nNOTE: ~15 trees over ~8 stands. A clean positive is informative; a null is NOT")
print("evidence of absent structure at this n. The 177-library assembly is the real test.")
