#!/usr/bin/env python3
"""Filter LOCI (not genotypes) by mean depth relative to the median locus.

A multi-copy / CNV locus attracts reads from several genomic copies, so its MEAN
depth sits above the median locus and its reads are a mixture regardless of how
many you keep -- which is why rarefaction did not remove the effect.

For each threshold, F is recomputed by genotype depth. If these loci are the
cause, removing them should flatten the depth gradient.
"""
import sys, math, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
E=0.001
def lnb(k,n,p):
    if p<=0: return 0.0 if k==0 else -1e18
    if p>=1: return 0.0 if k==n else -1e18
    return k*math.log(p)+(n-k)*math.log(1-p)
def call(ac,d):
    lhet=lnb(ac,d,0.5); lhom=max(lnb(ac,d,E),lnb(ac,d,1-E))
    return 1 if lhet>lhom else (0 if (d-ac)>=ac else 2)
EDG=[6,10,15,20,30,50,100,10**9]
NM=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
def dbin(d):
    for j in range(len(EDG)-1):
        if EDG[j]<=d<EDG[j+1]: return j

# pass 1: mean depth per locus, and the median across loci
recs=[]; means=[]
cols=None
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in stand]
            byst=collections.defaultdict(list)
            for i,s in cols: byst[stand[s]].append(i)
            continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        per=[]; alld=[]
        for st,idxs in byst.items():
            obs=[]; sref=salt=0
            for i in idxs:
                p=f[i].split(':')
                if len(p)<=ci: continue
                try: catg=[int(x) for x in p[ci].split(',')]
                except ValueError: continue
                rc,ac=catg[ri],catg[ai]; d=rc+ac
                if d<6: continue
                obs.append((d,rc,ac)); sref+=rc; salt+=ac; alld.append(d)
            if len(obs)<6: continue
            tot=sref+salt
            if tot==0: continue
            q=salt/tot
            if min(q,1-q)<=0: continue
            per.append((2*q*(1-q), obs))
        if not per or len(alld)<20: continue
        m=sum(alld)/len(alld)
        recs.append((m, per)); means.append(m)
means.sort(); med=means[len(means)//2]
print('loci: %d | median locus mean depth: %.1f'%(len(recs),med))
print('depth distribution of loci: p50 %.0f  p75 %.0f  p90 %.0f  p95 %.0f  p99 %.0f  max %.0f'%(
    med, means[int(len(means)*.75)], means[int(len(means)*.90)],
    means[int(len(means)*.95)], means[int(len(means)*.99)], means[-1]))

for mult in (None, 3.0, 2.0, 1.5, 1.25):
    acc=collections.defaultdict(lambda:[0,0.0,0]); kept=0
    for m,per in recs:
        if mult is not None and m > mult*med: continue
        kept+=1
        for he1,obs in per:
            for d,rc,ac in obs:
                b=dbin(d)
                if b is None: continue
                g=call(ac,d)
                a=acc[b]; a[0]+=(g==1); a[1]+=he1; a[2]+=1
    lab='no locus filter' if mult is None else 'mean depth <= %.2fx median (%.0fx)'%(mult,mult*med)
    vals=[]; row=[]
    for b in range(len(NM)):
        a=acc.get(b)
        if not a or a[1]<=0 or a[2]<500: row.append('%8s'%'-'); continue
        F=1-a[0]/a[1]; vals.append(F); row.append('%8.3f'%F)
    print('\n%-42s loci kept %6d (%.0f%%)'%(lab,kept,100*kept/len(recs)))
    print('  %-8s'%'depth', ' '.join('%8s'%n for n in NM))
    print('  %-8s'%'F_IS ', ' '.join(row), '   spread %.3f'%(max(vals)-min(vals) if len(vals)>1 else float('nan')))
