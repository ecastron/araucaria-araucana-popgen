#!/usr/bin/env python3
"""NULL MODEL for the F-vs-depth gradient.

Simulates a population in PERFECT Hardy-Weinberg -- no inbreeding, no paralogs,
no artifacts -- then applies our exact depth distribution, read sampling and
genotype-calling rule, and computes F by depth bin.

Whatever gradient appears here is MECHANICAL: it is what any dataset produces.
Our observed gradient is only meaningful to the extent it EXCEEDS this.

Allele frequencies and the depth distribution are taken from the real VCF so the
null is matched to our data rather than idealised.
"""
import sys, math, random, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
rng=random.Random(42)
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
BINS=[(6,10),(10,15),(15,20),(20,30),(30,50),(50,100),(100,10**9)]
NAMES=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
E=0.001
# harvest observed allele frequencies and the depth pool from the real data
freqs=[]; depths=[]
cols=None; nl=0
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in stand]; continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        aa=0; tot=0; dl=[]
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; d=rc+ac
            if d<6: continue
            dl.append(d); aa+=ac; tot+=d
        if len(dl)<20: continue
        nl+=1
        freqs.append(aa/tot); depths.extend(dl[:40])
        if nl>=8000: break
print('loci sampled for the null: %d ; depth pool: %d'%(len(freqs),len(depths)))
NIND=177
def lnb(k,n,p):
    if p<=0: return -1e18 if k>0 else 0.0
    if p>=1: return -1e18 if k<n else 0.0
    return k*math.log(p)+(n-k)*math.log(1-p)
acc=collections.defaultdict(lambda:[0,0.0])
for q in freqs:
    if q<=0 or q>=1: continue
    per=collections.defaultdict(lambda:[0,0,0])
    for _ in range(NIND):
        # TRUE genotype under HWE -- no inbreeding whatsoever
        u=rng.random()
        g = 0 if u<(1-q)**2 else (2 if u>1-q*q else 1)
        d=depths[rng.randrange(len(depths))]
        pr = E if g==0 else (1-E if g==2 else 0.5)
        ac=sum(1 for _ in range(d) if rng.random()<pr)
        k=next((j for j,(lo,hi) in enumerate(BINS) if lo<=d<hi),None)
        if k is None: continue
        lhet=lnb(ac,d,0.5); lhom=max(lnb(ac,d,E),lnb(ac,d,1-E))
        call = 1 if 2*(lhet-lhom)>0 else (0 if (d-ac)>=ac else 2)
        per[k][call]+=1
    for k,(a0,a1,a2) in per.items():
        n=a0+a1+a2
        if n<10: continue
        p_=(a1+2*a2)/(2*n)
        if min(p_,1-p_)<=0: continue
        acc[k][0]+=a1; acc[k][1]+=2*p_*(1-p_)*n
print()
print('%-10s %12s %14s'%('depth','NULL F (HWE)','OBSERVED F'))
OBS={'6-10':0.792,'10-15':0.814,'15-20':0.747,'20-30':0.654,'30-50':0.322,'50-100':-0.260,'100+':-0.191}
print('-'*40)
for k,nm in enumerate(NAMES):
    a=acc.get(k)
    if not a or a[1]<=0: continue
    print('%-10s %12.3f %14.3f'%(nm,1-a[0]/a[1],OBS.get(nm,float('nan'))))
print('\nIf NULL reproduces OBSERVED, the gradient is mechanical and means nothing.')
