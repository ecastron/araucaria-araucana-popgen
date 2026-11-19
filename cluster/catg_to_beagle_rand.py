#!/usr/bin/env python3
"""Random-thinning control for the depth test.

Genotypes whose informative depth (n = nRef + nAlt from CATG) falls below MINDP are emitted as
a FLAT likelihood (0.333/0.333/0.333) -- ngsRelate then draws no information from them, which is
the correct semantics for 'not confidently callable' as opposed to 'confidently homozygous'.

Purpose: the GL model assumes a heterozygote's reads are a fair binomial sample at p = 0.5. At
n = 6 a 6:0 split yields a ~60:1 likelihood ratio favouring homozygote, so systematic allelic
imbalance (PCR duplication, restriction-site nulls, capture bias) is indistinguishable from a
true homozygote. If F is driven by that, F must fall as MINDP rises. If F is real, it must not.
Usage: catg_to_beagle_mindp.py <vcf> <out> <err> <mindp>
"""
import math, re, sys, random
VCF, OUT, ERR, MINDP = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
KEEP = float(sys.argv[5]); random.seed(int(sys.argv[6]) if len(sys.argv)>6 else 1)
BASES = "CATG"; CODE = {"A":0,"C":1,"G":2,"T":3}
FLAT = ["0.333333","0.333333","0.333333"]

def binom_logpmf(k,n,p):
    if p<=0: return 0.0 if k==0 else -math.inf
    if p>=1: return 0.0 if k==n else -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
            +k*math.log(p)+(n-k)*math.log(1-p))

n_out=n_skip=0; kept=0; dropped=0; per_ind=None
with open(VCF) as fh, open(OUT,"w") as out:
    for line in fh:
        if line.startswith("##"): continue
        f=line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples=[re.sub(r"\.[12]$","",s) for s in f[9:]]
            per_ind=[0]*len(samples)
            hdr=["marker","allele1","allele2"]
            for s in samples: hdr+=[s,s,s]
            out.write("\t".join(hdr)+"\n"); continue
        ref,alt=f[3],f[4].split(",")[0]
        if ref not in BASES or alt not in BASES or ref==alt:
            n_skip+=1; continue
        ri,ai=BASES.index(ref),BASES.index(alt)
        row=[f"{f[0]}_{f[1]}",str(CODE[ref]),str(CODE[alt])]
        for j,c in enumerate(f[9:]):
            p=c.split(":")
            try: catg=[int(x) for x in p[2].split(",")]
            except (IndexError,ValueError): row+=FLAT; continue
            nr,na=catg[ri],catg[ai]; n=nr+na
            if n==0: row+=FLAT; continue
            if random.random() > KEEP: row+=FLAT; dropped+=1; continue
            kept+=1; per_ind[j]+=1
            ls=[binom_logpmf(na,n,ERR),binom_logpmf(na,n,0.5),binom_logpmf(na,n,1-ERR)]
            m=max(ls); ex=[math.exp(x-m) for x in ls]; t=sum(ex)
            row+=[f"{v/t:.6f}" for v in ex]
        out.write("\t".join(row)+"\n"); n_out+=1
import statistics as st
print(f"KEEP={KEEP}: {n_out} sites, {kept} informative genotypes kept, {dropped} zeroed")
print(f"KEEP={KEEP}: per-individual informative sites  median={st.median(per_ind)} "
      f"min={min(per_ind)} max={max(per_ind)} | inds with >=2000: {sum(1 for x in per_ind if x>=2000)}/{len(per_ind)}")
