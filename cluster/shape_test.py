#!/usr/bin/env python3
"""Depth-corrected SHAPE test for real heterozygotes, with split-sample validation.

For each genotype at >=MINDP with intermediate allele balance, z = (nAlt - n/2)/sqrt(n/4).
A true diploid heterozygote has z ~ N(0,1) at ANY depth. Mixed templates (collapsed paralogs,
contamination) give SD(z) >> 1. Per locus we take SD(z) over its intermediate genotypes.

NON-CIRCULAR BY CONSTRUCTION: loci are SELECTED using half the trees (A) and the balance
distribution is then examined in the OTHER half (B). Selecting on A cannot manufacture a peak
at 0.5 in B.
"""
import sys, random, statistics as st
from collections import defaultdict, Counter
V="rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
MINDP=20; LO,HI=0.10,0.90; MINGT=8
random.seed(1)

samples=[]; A=set(); B="CATG"
locZ=defaultdict(list)      # locus -> [z...] from half A
balB=defaultdict(list)      # locus -> [balance...] from half B
nsnp=0
for line in open(V):
    if line.startswith("##"): continue
    f=line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples=f[9:]
        idx=list(range(len(samples))); random.shuffle(idx)
        A=set(idx[:len(idx)//2])
        print(f"{len(samples)} libraries -> half A {len(A)}, half B {len(samples)-len(A)}")
        continue
    r,a=f[3],f[4].split(",")[0]
    if r not in B or a not in B or r==a: continue
    ri,ai=B.index(r),B.index(a); nsnp+=1
    loc=f[0]
    for k,c in enumerate(f[9:]):
        p=c.split(":")
        try: g=[int(x) for x in p[2].split(",")]
        except (IndexError,ValueError): continue
        nr,na=g[ri],g[ai]; n=nr+na
        if n<MINDP: continue
        b=na/n
        if not (LO<b<HI): continue
        if k in A: locZ[loc].append((na-n/2)/((n/4)**0.5))
        else:      balB[loc].append(b)

sd={L:st.pstdev(z) for L,z in locZ.items() if len(z)>=MINGT}
print(f"{nsnp} SNPs | loci scorable in half A (>= {MINGT} intermediate genotypes): {len(sd)}")
if not sd: sys.exit("no scorable loci")

print(f"\nSD(z) across loci: min {min(sd.values()):.2f} median {st.median(sd.values()):.2f} "
      f"max {max(sd.values()):.2f}   [real hets -> ~1.0]")
h=Counter()
for v in sd.values(): h[min(int(v/0.5),11)]+=1
print("\nSD(z) distribution (bin 0.5):")
for k in range(12):
    lab=f"{k*0.5:.1f}-{(k+1)*0.5:.1f}" if k<11 else "5.5+"
    v=h.get(k,0)
    print(f"  {lab:>9} {v:>6} {'#'*int(55*v/max(h.values()))}")

def histB(loci,lab):
    hh=Counter(); tot=0
    for L in loci:
        for b in balB.get(L,[]):
            hh[min(int(b*20),19)]+=1; tot+=1
    if tot<200: print(f"\n{lab}: only {tot} genotypes in half B -- too few"); return
    peak=sum(hh.get(k,0) for k in (9,10))/tot
    print(f"\n--- HALF B: {lab} ({len(loci)} loci, {tot} genotypes) ---")
    for k in range(20):
        v=hh.get(k,0)
        print(f"  {k*0.05:.2f}-{(k+1)*0.05:.2f} {100*v/tot:6.2f}% {'#'*int(45*v/max(hh.values()))}")
    print(f"  mass in 0.45-0.55 = {100*peak:.2f}%   (uniform over 0.10-0.90 would give ~12.5%)")

tight=[L for L,v in sd.items() if v<=1.5]
loose=[L for L,v in sd.items() if v>=3.0]
print(f"\nselected in half A: tight SD(z)<=1.5 -> {len(tight)} loci | loose SD(z)>=3.0 -> {len(loose)} loci")
histB(tight,"loci selected as TIGHT in half A")
histB(loose,"loci selected as LOOSE in half A (contrast)")
