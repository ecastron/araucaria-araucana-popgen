#!/usr/bin/env python3
"""Simulate genotypes under HARD-WIRED HWE (F = 0 by construction) and ask what F comes back out.

THE ONLY THING THAT CHANGES is the genotype truth. Everything that could plausibly manufacture a
heterozygote deficit is held at its observed value:
  * per-individual, per-site READ DEPTH  -- taken from the observed CATG counts
  * the MISSINGNESS PATTERN              -- a site with zero reads stays missing
  * the ALLELE FREQUENCY at each site    -- estimated from the observed reads
  * the GL MODEL and error rate          -- identical code path to catg_to_beagle.py, e = 0.01

So the simulated data have the same sparsity, the same depth heterogeneity and the same
uncertainty as the real data, and differ only in that no individual is inbred.

  F_sim ~= 0.70  ->  the pipeline/estimator manufactures the deficit from missingness and depth
  F_sim ~= 0     ->  the deficit is a property of the data, not of the estimator

Optional second arm: --dropout D makes each ALLELE independently invisible with probability D,
mimicking restriction-site polymorphism. A heterozygote losing one allele is read as a homozygote.
This is how much null-allele dropout would be needed to explain the observed F, if depth and
missingness alone cannot.
"""
import math, re, sys, random, argparse

ap = argparse.ArgumentParser()
ap.add_argument("vcf"); ap.add_argument("out")
ap.add_argument("--err", type=float, default=0.01)
ap.add_argument("--dropout", type=float, default=0.0)
ap.add_argument("--seed", type=int, default=2026)
ap.add_argument("--popmap", default=None,
                help="two-column lib<TAB>stand; draw HWE WITHIN stand at per-stand allele "
                     "frequencies, so the pooled estimator sees a real Wahlund effect")
a = ap.parse_args()
random.seed(a.seed)
BASES = "CATG"; CODE = {"A":0,"C":1,"G":2,"T":3}

def binom_logpmf(k, n, p):
    if p <= 0: return 0.0 if k == 0 else -math.inf
    if p >= 1: return 0.0 if k == n else -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
            + k*math.log(p) + (n-k)*math.log(1-p))

n_out = n_skip = 0
het_true = het_obs = 0
with open(a.vcf) as fh, open(a.out, "w") as out:
    for line in fh:
        if line.startswith("##"): continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$","",s) for s in f[9:]]
            if a.popmap:
                pm = {}
                for L in open(a.popmap):
                    q = L.split()
                    if len(q) >= 2: pm[re.sub(r"\.[12]$","",q[0])] = q[1]
                grp = [pm.get(s, "NA") for s in samples]
                groups = sorted(set(grp))
                gidx = {g: [i for i, x in enumerate(grp) if x == g] for g in groups}
                sys.stderr.write(f"popmap: {len(groups)} stands over {len(samples)} libs\n")
            hdr = ["marker","allele1","allele2"]
            for s in samples: hdr += [s,s,s]
            out.write("\t".join(hdr)+"\n"); continue
        ref, alt = f[3], f[4].split(",")[0]
        if ref not in BASES or alt not in BASES or ref == alt:
            n_skip += 1; continue
        ri, ai = BASES.index(ref), BASES.index(alt)

        depths = []
        tot_r = tot_a = 0
        for c in f[9:]:
            p = c.split(":")
            try: catg = [int(x) for x in p[2].split(",")]
            except (IndexError, ValueError): depths.append(0); continue
            nr, na = catg[ri], catg[ai]
            depths.append(nr+na); tot_r += nr; tot_a += na
        if tot_r + tot_a == 0: n_skip += 1; continue
        p_alt = tot_a / (tot_r + tot_a)                    # observed allele frequency
        pvec = [p_alt] * len(depths)
        if a.popmap:
            # per-stand frequency from that stand's own reads; falls back to pooled when a
            # stand has no reads at this site
            for g, idx in gidx.items():
                sr = sa = 0
                for i in idx:
                    c = f[9 + i].split(":")
                    try: catg = [int(x) for x in c[2].split(",")]
                    except (IndexError, ValueError): continue
                    sr += catg[ri]; sa += catg[ai]
                if sr + sa > 0:
                    for i in idx: pvec[i] = sa / (sr + sa)

        row = [f"{f[0]}_{f[1]}", str(CODE[ref]), str(CODE[alt])]
        for _i, n in enumerate(depths):
            p_alt = pvec[_i]
            if n == 0:
                row += ["0.333333","0.333333","0.333333"]; continue
            u = random.random()                            # HWE draw: F = 0 by construction
            if   u < (1-p_alt)**2:            g = 0
            elif u < (1-p_alt)**2 + 2*p_alt*(1-p_alt): g = 1
            else:                              g = 2
            if g == 1: het_true += 1
            if a.dropout > 0 and g == 1:                   # allele dropout -> looks homozygous
                lostR = random.random() < a.dropout
                lostA = random.random() < a.dropout
                if lostR and not lostA: g = 2
                elif lostA and not lostR: g = 0
            if g == 1: het_obs += 1
            theta = a.err if g == 0 else (0.5 if g == 1 else 1 - a.err)
            na = sum(1 for _ in range(n) if random.random() < theta)
            ls = [binom_logpmf(na,n,a.err), binom_logpmf(na,n,0.5), binom_logpmf(na,n,1-a.err)]
            m = max(ls); ex = [math.exp(x-m) for x in ls]; t = sum(ex)
            row += [f"{v/t:.6f}" for v in ex]
        out.write("\t".join(row)+"\n"); n_out += 1

print(f"simulated {n_out} sites x {len(samples)} individuals (skipped {n_skip})")
print(f"  dropout rate: {a.dropout}")
print(f"  true heterozygous calls: {het_true}; surviving dropout: {het_obs}"
      + (f" ({100*(1-het_obs/het_true):.1f}% lost)" if het_true else ""))
