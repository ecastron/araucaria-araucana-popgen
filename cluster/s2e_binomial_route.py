#!/usr/bin/env python3
"""Figure S2E on ONE genotype route: every filter and its unfiltered baseline on the same SNPs and bins.

PRODUCES  results/derived/s2e/binomial_route.tsv   one row per intervention (spread + per-bin F, n)
          stdout: the same table, human readable

ROUTE     "binomial all-SNP route": the ipyrad cr30 VCF, every biallelic SNP; each genotype re-called
          from its CATG read counts with the binomial rule used by scripts/ref/{rarefy_and_cv,
          locus_depth_filter, cnv_filter, F_by_copynumber}.py (het if L(0.5) > L(err=0.001));
          genotypes with < 6 reads ignored; F within stand, He from the stand's pooled read
          frequency (stand needs >= 6 genotypes and both alleles); a SNP enters only with >= 20
          genotypes at >= 6 reads. This is exactly the universe and estimator of those four
          scripts, so the unfiltered row must reproduce their 1.421 (7 bins) and 1.242 (4 bins >= 20x).

BINS      7 genotype-depth bins 6-10 ... 100+ (the S2B/S2D bins). A bin counts toward the spread
          only with >= MIN_N genotypes (500, the rule of locus_depth_filter.py and rarefy_and_cv.py).
          spread_ge20 = the same over the 4 bins >= 20x (the rarefy_and_cv.py bin set), reported
          because rarefaction to 30 reads cannot change bins below 30x.

ROWS      every S2E intervention that can be expressed on this route, each against the baseline
          row computed here. H_obs and the structural filter are re-expressed on binomial calls
          (their published values are polyRAD-route; see s2e_polyrad_route.R). The polyRAD model
          and the FreeBayes pipeline are route changes, not filters, and cannot be done here.

USAGE     s2e_binomial_route.py <cr30.vcf> <popmap_177.tsv> <loci_copies.tsv> <loci_sd.tsv> <out.tsv>
"""
import sys, math, random, collections, statistics as stt

IDX = {'C': 0, 'A': 1, 'T': 2, 'G': 3}
vcf, popmap, copies, sdf, out = sys.argv[1:6]
MIN_N, K = 500, 30
E = 0.001
EDG = [6, 10, 15, 20, 30, 50, 100, 10**9]
NM = ['6-10', '10-15', '15-20', '20-30', '30-50', '50-100', '100+']
GE20 = [3, 4, 5, 6]

stand = {}
for ln in open(popmap):
    p = ln.split()
    if len(p) >= 2: stand[p[0]] = p[1]
cn = {}
for ln in open(copies):
    p = ln.split()
    if len(p) >= 2: cn[p[0]] = int(p[1])
sd = {}
for ln in open(sdf):
    p = ln.split()
    if len(p) >= 2 and p[1] in ('0', '1'): sd[p[0]] = int(p[1])


def lnb(k, n, p):
    if p <= 0: return 0.0 if k == 0 else -1e18
    if p >= 1: return 0.0 if k == n else -1e18
    return k * math.log(p) + (n - k) * math.log(1 - p)


def call(ac, d):
    lhet = lnb(ac, d, 0.5); lhom = max(lnb(ac, d, E), lnb(ac, d, 1 - E))
    return 1 if lhet > lhom else (0 if (d - ac) >= ac else 2)


def dbin(d):
    for j in range(len(EDG) - 1):
        if EDG[j] <= d < EDG[j + 1]: return j


# ---- one pass: keep, per SNP, the per-stand observations and the locus-level metrics ----------
# rarefied calls are drawn here, in the same order and with the same seed as rarefy_and_cv.py
# (seed 42, a draw for every genotype > 30 reads, stands in header order), so its 1.093 must reproduce.
rng = random.Random(42)
recs = []   # (chrom, metrics dict, per=[(he1, [(d, g, g_rarefied), ...]), ...])
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols = [(i, s.split('.')[0]) for i, s in enumerate(f[9:], 9)]
            cols = [(i, s) for i, s in cols if s in stand]
            byst = collections.defaultdict(list)
            for i, s in cols: byst[stand[s]].append(i)
            continue
        if len(f[4]) != 1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt = f[8].split(':')
        try: ci = fmt.index('CATG')
        except ValueError: continue
        ri, ai = IDX[f[3]], IDX[f[4]]
        alld, present, nhet, ncall = [], [], 0, 0
        raw = {}
        for st, idxs in byst.items():
            for i in idxs:
                p = f[i].split(':')
                if len(p) <= ci: continue
                try: catg = [int(x) for x in p[ci].split(',')]
                except ValueError: continue
                rc, ac = catg[ri], catg[ai]; d = rc + ac
                if d > 0: present.append(d)
                if d >= 6:
                    alld.append(d); raw.setdefault(st, []).append((d, rc, ac))
        if len(alld) < 20: continue
        mu = sum(alld) / len(alld)
        cv = math.sqrt(sum((x - mu) ** 2 for x in alld) / len(alld)) / mu
        lmed = stt.median(alld)
        per = []
        for st in byst:                    # header order, as in the source scripts
            obs = raw.get(st, [])
            if len(obs) < 6: continue
            sref = sum(o[1] for o in obs); salt = sum(o[2] for o in obs)
            q = salt / (sref + salt)
            if min(q, 1 - q) <= 0: continue
            rows = []
            for d, rc, ac in obs:
                g = call(ac, d)
                if d > K:
                    s = rng.sample([1] * ac + [0] * rc, K); g2 = call(sum(s), K)
                else:
                    g2 = g
                rows.append((d, g, g2))
            per.append((2 * q * (1 - q), rows))
        # locus-level H_obs on binomial calls, over all genotypes >= 6 reads
        for st, obs in raw.items():
            for d, rc, ac in obs:
                ncall += 1; nhet += (call(ac, d) == 1)
        if not per: continue
        pm = sum(present) / len(present)
        pcv = (math.sqrt(sum((x - pm) ** 2 for x in present) / (len(present) - 1)) / pm
               if len(present) > 1 else 0.0)   # sd() in R is n-1, as locus_structural_filter.R
        met = dict(mean=mu, cv=cv, lmed=lmed, hobs=nhet / ncall,
                   pmean=pm, pcv=pcv, miss=1 - len(present) / len(cols),
                   cn=cn.get(f[0], 0), sd=sd.get(f[0]))
        recs.append((f[0], met, per))

med = sorted(m['mean'] for _, m, _ in recs)[len(recs) // 2]


def F_by_bin(keep_snp=lambda m: True, keep_geno=lambda m, d: True, rarefied=False):
    acc = [[0, 0.0, 0] for _ in NM]; nsnp = 0; nloc = set()
    for chrom, m, per in recs:
        if not keep_snp(m): continue
        nsnp += 1; nloc.add(chrom)
        for he1, rows in per:
            for d, g, g2 in rows:
                if not keep_geno(m, d): continue
                b = dbin(d); a = acc[b]
                a[0] += ((g2 if rarefied else g) == 1); a[1] += he1; a[2] += 1
    F = [1 - a[0] / a[1] if a[2] >= MIN_N and a[1] > 0 else None for a in acc]
    return F, [a[2] for a in acc], nsnp, len(nloc)


def spread(F, bins=range(7)):
    v = [F[b] for b in bins if F[b] is not None]
    return (max(v) - min(v)) if len(v) > 1 else float('nan'), len(v)


TESTS = [
    ("none - unfiltered baseline", {}),
    ("rarefy every genotype to 30 reads", dict(rarefied=True)),
    ("locus mean depth <= 1.25x median", dict(keep_snp=lambda m: m['mean'] <= 1.25 * med)),
    ("locus mean depth <= 1.50x median", dict(keep_snp=lambda m: m['mean'] <= 1.5 * med)),
    ("locus mean depth <= 3.00x median", dict(keep_snp=lambda m: m['mean'] <= 3.0 * med)),
    ("H_obs ceiling 0.45 (binomial calls)", dict(keep_snp=lambda m: m['hobs'] <= 0.45)),
    ("H_obs ceiling 0.20 (binomial calls)", dict(keep_snp=lambda m: m['hobs'] <= 0.20)),
    ("depth-CV filter (CV <= 1.0)", dict(keep_snp=lambda m: m['cv'] <= 1.0)),
    ("locus-relative depth (CNV) filter <= 1.5", dict(keep_geno=lambda m, d: d / m['lmed'] <= 1.5)),
    ("structural filter cv<=0.80 depth<=95 miss<=0.70",
     dict(keep_snp=lambda m: m['pcv'] <= 0.80 and m['pmean'] <= 95 and m['miss'] <= 0.70)),
    ("structural filter cv<=0.80 depth<=64 miss<=0.65",
     dict(keep_snp=lambda m: m['pcv'] <= 0.80 and m['pmean'] <= 64 and m['miss'] <= 0.65)),
    ("single-copy loci only (genome mapping)", dict(keep_snp=lambda m: m['cn'] == 1)),
    ("2-copy loci only", dict(keep_snp=lambda m: m['cn'] == 2)),
    ("3-5-copy loci only", dict(keep_snp=lambda m: 3 <= m['cn'] <= 5)),
    ("6+-copy loci only", dict(keep_snp=lambda m: m['cn'] >= 6)),
    ("exclude segmental-duplication loci", dict(keep_snp=lambda m: m['sd'] == 0)),
]

cols_out = ["mechanism", "route", "n_snps_kept", "n_loci_kept", "pct_snps_kept", "n_bins", "baseline_spread",
            "filtered_spread", "reduction", "n_bins_ge20", "baseline_spread_ge20", "filtered_spread_ge20",
            "reduction_ge20", "genotypes_used"] + [f"F_{b}" for b in NM] + [f"n_{b}" for b in NM]
res = []
base = None
for lab, kw in TESTS:
    F, n, nsnp, nloc = F_by_bin(**kw)
    s7, b7 = spread(F); s4, b4 = spread(F, GE20)
    if base is None: base = (s7, s4, nsnp)
    res.append([lab, "binomial all-SNP (ipyrad cr30 VCF)", nsnp, nloc, round(100 * nsnp / base[2], 1), b7,
                round(base[0], 4), round(s7, 4), round(base[0] - s7, 4), b4, round(base[1], 4), round(s4, 4),
                round(base[1] - s4, 4), sum(n)]
               + [("NA" if x is None else round(x, 4)) for x in F] + n)

with open(out, "w") as fo:
    fo.write("\t".join(cols_out) + "\n")
    for r in res: fo.write("\t".join(map(str, r)) + "\n")
print(f"SNPs in universe: {len(recs)} | median SNP mean depth {med:.1f} | min genotypes per bin {MIN_N}")
print("%-50s %7s %3s %8s %8s %8s | %8s %8s" % ("mechanism", "SNPs", "bins", "base", "filt", "reduct",
                                                "base>=20", "filt>=20"))
for r in res:
    print("%-50s %7d %3d %8.3f %8.3f %8.3f | %8.3f %8.3f" % (r[0], r[2], r[5], r[6], r[7], r[8], r[10], r[11]))
print(f"written -> {out}")
