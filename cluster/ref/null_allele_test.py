#!/usr/bin/env python3
"""Null alleles from restriction-site polymorphism?

If a SNP disrupts the PstI site, that allele is never cut and never sequenced.
Trees heterozygous for it look homozygous; trees homozygous for the uncut allele
have NO DATA at the locus. So null alleles predict a POSITIVE association between
a locus's missing-data rate and its heterozygote deficit.

Real biology predicts NO association: a locus being absent in some individuals
says nothing about whether the individuals that DO have data are in HWE.

Control for depth: low-depth loci have both more missingness and more undercalling,
which would create the association spuriously. So the association is also computed
WITHIN depth strata, holding mean locus depth constant.
"""
import sys, collections, math
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]

rows=[]
cols=None
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in stand]
            ntot=len(cols); print('libraries: %d'%ntot); continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: gi=fmt.index('GT'); di=fmt.index('DP')
        except ValueError: continue
        naa=nab=nbb=miss=0; dpsum=0; ndp=0
        for i,s in cols:
            p=f[i].split(':')
            gt=p[gi]
            if gt in ('./.','.|.'): miss+=1; continue
            try: dp=int(p[di])
            except (ValueError,IndexError): dp=0
            dpsum+=dp; ndp+=1
            if gt in ('0/0','0|0'): naa+=1
            elif gt in ('1/1','1|1'): nbb+=1
            else: nab+=1
        n=naa+nab+nbb
        if n<30: continue
        p_=(2*naa+nab)/(2*n); q=1-p_
        if min(p_,q)<=0.02: continue          # need real polymorphism to measure F
        hexp=2*p_*q*n
        if hexp<1: continue
        F=1-nab/hexp
        rows.append((miss/ntot, F, dpsum/ndp if ndp else 0))
print('loci usable: %d'%len(rows))

def pearson(xs,ys):
    n=len(xs)
    if n<10: return float('nan')
    mx=sum(xs)/n; my=sum(ys)/n
    sx=math.sqrt(sum((x-mx)**2 for x in xs)); sy=math.sqrt(sum((y-my)**2 for y in ys))
    if sx==0 or sy==0: return float('nan')
    return sum((x-mx)*(y-my) for x,y in zip(xs,ys))/(sx*sy)

miss=[r[0] for r in rows]; F=[r[1] for r in rows]
print('\n=== OVERALL ===')
print('  correlation( missing-data rate , F ) = %+.4f'%pearson(miss,F))
print('\n=== F by missing-data decile ===')
rows.sort(key=lambda r:r[0])
k=len(rows)//10
print('%-16s %8s %10s %10s'%('missing rate','n_loci','mean F','mean depth'))
for d in range(10):
    ch=rows[d*k:(d+1)*k] if d<9 else rows[9*k:]
    if not ch: continue
    print('%-16s %8d %10.4f %10.1f'%('%.3f-%.3f'%(ch[0][0],ch[-1][0]),len(ch),
          sum(r[1] for r in ch)/len(ch), sum(r[2] for r in ch)/len(ch)))
print('\n=== CONTROL: same association WITHIN depth strata ===')
print('(if the overall correlation is just depth, it vanishes here)')
strata=[(0,10),(10,15),(15,25),(25,40),(40,10**9)]
print('%-12s %8s %12s'%('mean depth','n_loci','corr(miss,F)'))
for lo,hi in strata:
    sub=[r for r in rows if lo<=r[2]<hi]
    if len(sub)<200: continue
    print('%-12s %8d %+12.4f'%('%d-%d'%(lo,min(hi,999)),len(sub),
          pearson([r[0] for r in sub],[r[1] for r in sub])))
