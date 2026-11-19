#!/usr/bin/env python3
"""Two-group F_ST: Coastal vs Andes, plus the hierarchical partition.

Read-based Hudson estimator (ratio of averages, Bhatia et al. 2013), computed from
CATG read counts so it never passes through a genotype call -- appropriate given the
documented heterozygote-calling artifact in this dataset.
Block jackknife over loci for the SE.
"""
import sys, math, collections, random
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
COAST={'PN','TR','VA'}
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
num=[]; den=[]; loci=[]
cols=None
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in stand]
            grp=collections.defaultdict(list)
            for i,s in cols: grp['C' if stand[s] in COAST else 'A'].append(i)
            continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        p={} ; n={}
        for g,idxs in grp.items():
            a=b=0
            for i in idxs:
                q=f[i].split(':')
                if len(q)<=ci: continue
                try: catg=[int(x) for x in q[ci].split(',')]
                except ValueError: continue
                a+=catg[ri]; b+=catg[ai]
            if a+b>=40: p[g]=b/(a+b); n[g]=a+b
        if len(p)<2: continue
        x,y=p['C'],p['A']; nx,ny=n['C'],n['A']
        nn=(x-y)**2 - x*(1-x)/(nx-1) - y*(1-y)/(ny-1)
        dd=x*(1-y)+y*(1-x)
        if dd>0:
            num.append(nn); den.append(dd); loci.append(f[0])
N=len(num)
fst=sum(num)/sum(den)
# block jackknife over 100 blocks of contiguous loci
B=100; blk=[[] for _ in range(B)]
for k in range(N): blk[k*B//N].append(k)
vals=[]
for b in range(B):
    keep=set(range(N))-set(blk[b])
    if not keep: continue
    vals.append(sum(num[k] for k in keep)/sum(den[k] for k in keep))
m=sum(vals)/len(vals)
se=math.sqrt((len(vals)-1)/len(vals)*sum((v-m)**2 for v in vals))
print("COASTAL vs ANDES, two-group Hudson F_ST (read-based)")
print("  SNPs used   : %d"%N)
print("  F_ST        : %.4f"%fst)
print("  jackknife SE: %.4f   95%% CI %.4f - %.4f"%(se,fst-1.96*se,fst+1.96*se))
