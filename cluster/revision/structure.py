#!/usr/bin/env python
"""Revision (25 Nov 2026): robustness of population structure to linkage, sample size and
block definition.

Input   rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf, 163 trees in
        results/Q1/feems/keep.txt, biallelic SNPs only (78,026; as load_observed.py).
        CHROM in this VCF is the ipyrad RAD locus, so SNPs sharing CHROM are physically linked.

Analyses
  1. one SNP per RAD locus (random within locus, seed 1) -> thinned genotype set
  2. pairwise Hudson FST (project estimator, sumstats.hudson_fst; >= 6 called chromosomes per
     stand) on all SNPs and on the thinned set, with delete-one jackknife over 50 blocks whose
     edges fall on RAD-locus boundaries (equal numbers of loci per block)
  3. the within/between gap: largest within-range vs smallest between-range FST, its jackknife
     SE, and the same with all pairs involving VA removed
  4. equal sample size: 9 trees per stand (the smallest stand, PC), 200 draws, seed 2;
     share of draws in which every between-range FST exceeds every within-range FST
  5. f3(target; A, B) on the 270 cross-range-informative triples, all SNPs, locus-edged blocks
     (same estimator as abc/f3_targeted.py)
  6. the 12 differentiation summary statistics of the ABC (sumstats2.summarise) on the thinned set
  7. per-tree call rate by stand and range, and a stand-level permutation test of the range
     difference (all 220 assignments of 3 of 12 stands)

Outputs (results/revision/)
  structure.json            all numbers above
  fst_pairs_locusblocks.tsv     pairwise FST and locus-block jackknife CI, all SNPs and thinned
  thin1.extract                 plink variant IDs of the thinned set (for PCA/ADMIXTURE)
  keep163.plink                 FID IID for plink --keep (VCF sample names)
  thin1.geno, thin1.geno.samples  LEA input for sNMF on the thinned set
  callrate_by_tree.tsv
"""
import sys, csv, re, json, itertools as it
import numpy as np
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
from sumstats import freqs_by_deme, hudson_fst
from sumstats2 import summarise, STAT_ORDER

P = "/home/ecastron/araucaria"
V = f"{P}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
OUT = f"{P}/results/revision"
COAST = {"VA", "TR", "PN"}
NB = 50

site = {r[0]: r[1] for r in csv.reader(open(f"{P}/metadata/popmap_177.tsv"), delimiter="\t") if len(r) >= 2}
keep = {l.split()[0] for l in open(f"{P}/results/Q1/feems/keep.txt")}
lat, lon = {}, {}
for r in csv.DictReader(open(f"{P}/metadata/manifest rad runs.csv")):
    try:
        lat[r["run_accession"]] = float(r["latitude"]); lon[r["run_accession"]] = float(r["longitude"])
    except ValueError:
        pass

# ---- read genotypes (same rules as load_observed.py) plus locus and variant IDs ----
rows, chrom, vid = [], [], []
for line in open(V):
    if line.startswith("##"): continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        raw = f[9:]
        samples = [re.sub(r"\.[12]$", "", s) for s in raw]
        cols = [i for i, s in enumerate(samples) if s in keep and s in site]
        continue
    if len(f[3]) != 1 or len(f[4]) != 1 or "," in f[4]: continue
    g = []
    for i in cols:
        t = f[9 + i].split(":", 1)[0].replace("|", "/")
        if t in ("./.", "."): g.append(-1)
        else:
            a, b = t.split("/"); g.append(int(a != "0") + int(b != "0"))
    rows.append(g); chrom.append(f[0]); vid.append(f[2])
G = np.array(rows, dtype=np.int8).T.astype(float); G[G < 0] = np.nan
sel = [samples[i] for i in cols]; selraw = [raw[i] for i in cols]
stands = sorted({site[s] for s in sel})
deme = np.array([stands.index(site[s]) for s in sel])
coords = np.array([[np.mean([lon[s] for s in sel if site[s] == st]),
                    np.mean([lat[s] for s in sel if site[s] == st])] for st in stands])
coastal = np.array([st in COAST for st in stands])
chrom = np.array(chrom); L = G.shape[1]
assert G.shape == (163, 78026), G.shape
res = {"n_trees": int(G.shape[0]), "n_snps": int(L)}

# ---- locus structure and the thinned set ----
loc_ids, loc_idx = np.unique(chrom, return_inverse=True)
order = np.argsort(np.array([int(c.split("_")[1]) for c in loc_ids]))   # RAD_1, RAD_2, ... order
rank = np.empty_like(order); rank[order] = np.arange(len(order))
loc_rank = rank[loc_idx]                       # rank of each SNP's locus in VCF order
nloc = len(loc_ids)
spl = np.bincount(loc_idx)
res["n_loci"] = int(nloc)
res["snps_per_locus"] = {"mean": float(spl.mean()), "median": float(np.median(spl)),
                         "max": int(spl.max()), "share_loci_1snp": float((spl == 1).mean())}
rng = np.random.default_rng(1)
# one random SNP per locus: random key per SNP, keep the SNP with the smallest key in each locus
key = rng.random(L)
o = np.lexsort((key, loc_idx))                  # grouped by locus, random order within
first = np.r_[True, loc_idx[o][1:] != loc_idx[o][:-1]]
pick = np.zeros(L, bool); pick[o[first]] = True
assert pick.sum() == nloc
res["n_snps_thin1"] = int(pick.sum())
with open(f"{OUT}/thin1.extract", "w") as fh:
    fh.write("\n".join(np.array(vid)[pick]) + "\n")
with open(f"{OUT}/keep163.plink", "w") as fh:
    fh.write("".join(f"{s}\t{s}\n" for s in selraw))
for tag, js in (("thin1", np.flatnonzero(pick)), ("all163", np.arange(L))):
    with open(f"{OUT}/{tag}.geno", "w") as fh:      # LEA .geno: one line per SNP, 9 = missing
        for j in js:
            fh.write("".join("9" if np.isnan(x) else str(int(x)) for x in G[:, j]) + "\n")
    with open(f"{OUT}/{tag}.geno.samples", "w") as fh:
        fh.write("\n".join(sel) + "\n")

R = 6371.0
def km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, (coords[a, 1], coords[a, 0], coords[b, 1], coords[b, 0]))
    h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(h))

def pair_fst(Gs, lrank, with_jk=True):
    """pairwise FST for 66 pairs; jackknife blocks = equal numbers of loci, edges on loci"""
    Pf, Nf = freqs_by_deme(Gs, deme, 12)
    blk = (lrank * NB) // nloc
    out = []
    for a, b in it.combinations(range(12), 2):
        m = (Nf[a] >= 6) & (Nf[b] >= 6)
        f = hudson_fst(Pf[a][m], Pf[b][m], Nf[a][m], Nf[b][m])
        se = np.nan
        if with_jk:
            vals = []
            for k in range(NB):
                mm = m & (blk != k)
                vals.append(hudson_fst(Pf[a][mm], Pf[b][mm], Nf[a][mm], Nf[b][mm]))
            vals = np.array(vals); n = len(vals)
            se = float(np.sqrt((n - 1) / n * np.sum((vals - vals.mean()) ** 2)))
        out.append(dict(a=stands[a], b=stands[b], km=float(km(a, b)), fst=float(f), se=se,
                        between=bool(coastal[a] != coastal[b]), nloci=int(m.sum())))
    return out

def gap(rows, drop=None):
    rr = [r for r in rows if not (drop and (r["a"] == drop or r["b"] == drop))]
    w = [r for r in rr if not r["between"]]; b = [r for r in rr if r["between"]]
    wmax = max(w, key=lambda r: r["fst"]); bmin = min(b, key=lambda r: r["fst"])
    d = bmin["fst"] - wmax["fst"]
    se = float(np.sqrt(wmax["se"] ** 2 + bmin["se"] ** 2)) if np.isfinite(wmax["se"]) else np.nan
    ov = sum(1 for x in w for y in b
             if np.isfinite(x["se"]) and x["fst"] + 1.96 * x["se"] >= y["fst"] - 1.96 * y["se"])
    return dict(n_within=len(w), n_between=len(b),
                within_max=[wmax["a"], wmax["b"], wmax["fst"], wmax["se"]],
                between_min=[bmin["a"], bmin["b"], bmin["fst"], bmin["se"]],
                gap=float(d), gap_se=se, gap_z=float(d / se) if se and np.isfinite(se) else None,
                ci_overlapping_pairs=int(ov), total_pairs=len(w) * len(b))

def ibd(rows):
    x = np.log([r["km"] for r in rows]); y = np.array([r["fst"] for r in rows])
    return float(np.corrcoef(x, y)[0, 1])

full = pair_fst(G, loc_rank)
thin = pair_fst(G[:, pick], loc_rank[pick])
res["all_snps"] = {"gap": gap(full), "gap_without_VA": gap(full, "VA"), "ibd_r": ibd(full)}
res["thin1"] = {"gap": gap(thin), "gap_without_VA": gap(thin, "VA"), "ibd_r": ibd(thin),
                "r_fst_with_all_snps": float(np.corrcoef([r["fst"] for r in full], [r["fst"] for r in thin])[0, 1])}
with open(f"{OUT}/fst_pairs_locusblocks.tsv", "w") as fh:
    fh.write("set\tstand_a\tstand_b\tkm\tfst\tse\tlo95\thi95\tbetween\tnloci\n")
    for tag, rows_ in (("all", full), ("thin1", thin)):
        for r in rows_:
            fh.write(f"{tag}\t{r['a']}\t{r['b']}\t{r['km']:.3f}\t{r['fst']:.6f}\t{r['se']:.6f}\t"
                     f"{r['fst']-1.96*r['se']:.6f}\t{r['fst']+1.96*r['se']:.6f}\t{int(r['between'])}\t{r['nloci']}\n")

# ---- equal sample size per stand ----
nmin = int(np.bincount(deme).min())
rng2 = np.random.default_rng(2)
sep, gaps, ibds = 0, [], []
NDRAW = 200
for _ in range(NDRAW):
    idx = np.concatenate([rng2.choice(np.flatnonzero(deme == d), nmin, replace=False) for d in range(12)])
    Pf, Nf = freqs_by_deme(G[idx], deme[idx], 12)
    w, b, xs, ys = [], [], [], []
    for a, c in it.combinations(range(12), 2):
        m = (Nf[a] >= 6) & (Nf[c] >= 6)
        f = hudson_fst(Pf[a][m], Pf[c][m], Nf[a][m], Nf[c][m])
        (b if coastal[a] != coastal[c] else w).append(f)
        xs.append(np.log(km(a, c))); ys.append(f)
    gaps.append(min(b) - max(w)); sep += min(b) > max(w); ibds.append(np.corrcoef(xs, ys)[0, 1])
res["equal_n"] = {"trees_per_stand": nmin, "draws": NDRAW, "share_separated": sep / NDRAW,
                  "gap_median": float(np.median(gaps)), "gap_2.5": float(np.percentile(gaps, 2.5)),
                  "gap_97.5": float(np.percentile(gaps, 97.5)), "ibd_r_median": float(np.median(ibds))}

# ---- f3 with locus-edged blocks (estimator as abc/f3_targeted.py) ----
Pf, Nf = freqs_by_deme(G, deme, 12)
blk = (loc_rank * NB) // nloc
f3rows = []
for c in range(12):
    for a, b in it.combinations([i for i in range(12) if i != c], 2):
        cross = (coastal[c] != coastal[a]) ^ (coastal[c] != coastal[b])
        m = np.isfinite(Pf[a]) & np.isfinite(Pf[b]) & np.isfinite(Pf[c]) & (Nf[c] > 2)
        num = (Pf[c] - Pf[a]) * (Pf[c] - Pf[b]) - Pf[c] * (1 - Pf[c]) / np.maximum(Nf[c] - 1, 1)
        est = num[m].mean()
        vals = np.array([num[m & (blk != k)].mean() for k in range(NB)])
        n = len(vals)
        se = np.sqrt((n - 1) / n * np.sum((vals - vals.mean()) ** 2))
        f3rows.append((cross, est, est / se))
for lab, want in (("cross", True), ("within", False)):
    g = [r for r in f3rows if r[0] == want]
    z = np.array([r[2] for r in g]); e = np.array([r[1] for r in g])
    res[f"f3_{lab}"] = {"n": len(g), "negative": int((e < 0).sum()), "z_lt_-3": int(((e < 0) & (z < -3)).sum()),
                        "min_z": float(z.min())}

# ---- ABC differentiation statistics on the thinned set ----
d = summarise(G[:, pick], deme, coords, coastal, nd=12)
res["abc_stats_thin1"] = {k: (float(d[k]) if d and k in d and np.isfinite(d[k]) else None) for k in STAT_ORDER}
d0 = summarise(G, deme, coords, coastal, nd=12)
res["abc_stats_all"] = {k: (float(d0[k]) if d0 and k in d0 and np.isfinite(d0[k]) else None) for k in STAT_ORDER}

# ---- call rate by tree, stand and range ----
cr = 1 - np.isnan(G).mean(axis=1)
with open(f"{OUT}/callrate_by_tree.tsv", "w") as fh:
    fh.write("tree\tstand\trange\tcall_rate\n")
    for s, d_, c in zip(sel, deme, cr):
        fh.write(f"{s}\t{stands[d_]}\t{'Coastal' if coastal[d_] else 'Andean'}\t{c:.4f}\n")
smean = np.array([cr[deme == d_].mean() for d_ in range(12)])
obs = smean[coastal].mean() - smean[~coastal].mean()
null = [smean[list(k)].mean() - np.delete(smean, list(k)).mean() for k in it.combinations(range(12), 3)]
res["callrate"] = {"by_stand": {stands[i]: float(smean[i]) for i in range(12)},
                   "coastal_mean": float(smean[coastal].mean()), "andean_mean": float(smean[~coastal].mean()),
                   "tree_min": float(cr.min()), "tree_median": float(np.median(cr)), "tree_max": float(cr.max()),
                   "perm_p_two_sided": float(np.mean(np.abs(null) >= abs(obs) - 1e-12))}

json.dump(res, open(f"{OUT}/structure.json", "w"), indent=1)
print(json.dumps(res, indent=1))
