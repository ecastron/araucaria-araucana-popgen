#!/usr/bin/env python3
"""Rebuild genotype likelihoods from HIGH-DEPTH genotypes only, to test the F = 0.70 estimate.

ngsRelate inferred median inbreeding F = 0.70 from the cr30 likelihoods. That is impossible in a
dioecious species, and the suspected cause is heterozygote under-calling at low coverage: at 6-11x
observed heterozygosity is 0.013, rising to 0.115 at 60x+. Within stands at >=20x it is a normal
0.17-0.27.

If low-depth sites are the cause, restricting the likelihoods to well-covered genotypes should
pull F down sharply. If F stays high, the explanation is something else and the kinship caveat is
worse than described.

Genotypes below the depth floor are emitted as flat (1/3, 1/3, 1/3) rather than dropped, so the
site stays in the matrix and only the uninformative cells are neutralised.
"""
import math, re, sys
VCF, OUT, MINDP = sys.argv[1], sys.argv[2], int(sys.argv[3])
ERR = 0.01; BASES = "CATG"; CODE = {"A":0,"C":1,"G":2,"T":3}
def lp(k, n, p):
    if p <= 0: return 0.0 if k == 0 else -math.inf
    if p >= 1: return 0.0 if k == n else -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
            + k*math.log(p) + (n-k)*math.log(1-p))
kept = flat = tot = 0
with open(VCF) as fh, open(OUT,"w") as out:
    for line in fh:
        if line.startswith("##"): continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples=[re.sub(r"\.[12]$","",s) for s in f[9:]]
            out.write("\t".join(["marker","allele1","allele2"]+[s for s in samples for _ in range(3)])+"\n")
            continue
        ref,alt = f[3], f[4].split(",")[0]
        if ref not in BASES or alt not in BASES or ref==alt: continue
        ri,ai = BASES.index(ref), BASES.index(alt)
        row=[f"{f[0]}_{f[1]}", str(CODE[ref]), str(CODE[alt])]
        for c in f[9:]:
            p=c.split(":")
            try: catg=[int(x) for x in p[2].split(",")]
            except (IndexError,ValueError):
                row += ["0.333333"]*3; flat+=1; tot+=1; continue
            nr,na = catg[ri], catg[ai]; n = nr+na; tot+=1
            if n < MINDP:
                row += ["0.333333"]*3; flat+=1; continue
            ls=[lp(na,n,ERR), lp(na,n,0.5), lp(na,n,1-ERR)]
            mx=max(ls); ex=[math.exp(x-mx) for x in ls]; t=sum(ex)
            row += [f"{v/t:.6f}" for v in ex]
        out.write("\t".join(row)+"\n"); kept+=1
print(f"depth floor {MINDP}x: {kept} sites, {100*flat/tot:.1f}% of genotypes flattened as uninformative")
