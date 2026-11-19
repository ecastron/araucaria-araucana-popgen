#!/usr/bin/env python3
"""Split the VCF into 'clean' and 'plateau' loci.

A locus is scored by the fraction of its >=20x genotypes whose allele balance falls in the flat
plateau (0.15-0.85). Real heterozygotes cluster at 0.5; the plateau is where collapsed paralogs
and mixed templates live. Loci needing >=MINGT genotypes at >=20x to be scorable.
Usage: clean_loci.py <in.vcf> <out_clean.vcf> <out_plateau.vcf> <maxrate> <mingt>
"""
import sys
from collections import defaultdict
IN,OC,OP,MAXRATE,MINGT = sys.argv[1],sys.argv[2],sys.argv[3],float(sys.argv[4]),int(sys.argv[5])
B="CATG"
stat=defaultdict(lambda:[0,0])
hdr=[]
for line in open(IN):
    if line.startswith("#"):
        hdr.append(line); continue
    f=line.rstrip("\n").split("\t")
    r,a=f[3],f[4].split(",")[0]
    if r not in B or a not in B or r==a: continue
    ri,ai=B.index(r),B.index(a)
    s=stat[f[0]]
    for c in f[9:]:
        p=c.split(":")
        try: g=[int(x) for x in p[2].split(",")]
        except (IndexError,ValueError): continue
        n=g[ri]+g[ai]
        if n<20: continue
        s[0]+=1
        if 0.15 < g[ai]/n < 0.85: s[1]+=1

clean=set(); plateau=set(); unscored=0
for loc,(tot,pl) in stat.items():
    if tot<MINGT: unscored+=1; continue
    (clean if pl/tot<=MAXRATE else plateau).add(loc)
print(f"loci: {len(clean)} clean (<= {MAXRATE}), {len(plateau)} plateau, {unscored} unscored (<{MINGT} gt at 20x)")

nc=np=0
with open(OC,"w") as oc, open(OP,"w") as op:
    for h in hdr: oc.write(h); op.write(h)
    for line in open(IN):
        if line.startswith("#"): continue
        loc=line.split("\t",1)[0]
        if loc in clean: oc.write(line); nc+=1
        elif loc in plateau: op.write(line); np+=1
print(f"SNPs written: clean={nc}  plateau={np}")
