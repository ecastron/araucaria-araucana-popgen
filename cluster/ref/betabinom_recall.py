#!/usr/bin/env python3
"""Re-call genotypes under a BETA-BINOMIAL heterozygote model and test whether
F becomes depth-INDEPENDENT.

Binomial het model: reads ~ Bin(n, 0.5).  Measured overdispersion is 2-4.7x, so
that is misspecified. Beta-binomial het model: p ~ Beta(s,s), reads ~ Bin(n,p).
s -> infinity recovers the binomial; small s allows the observed scatter.

CRITERION: not "does F reach 0" (a free beta-binomial can always achieve that by
making balance uninformative -- that is fitting the artifact). The test is whether
F becomes FLAT ACROSS DEPTH. A correct model must give the same answer at every depth.
"""
import sys, math, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
EDGES=[6,10,15,20,30,50,100,10**9]
NAMES=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
def binof(d):
    for j in range(len(EDGES)-1):
        if EDGES[j]<=d<EDGES[j+1]: return j
    return None
E=0.001
lg=math.lgamma
def lbeta(a,b): return lg(a)+lg(b)-lg(a+b)
def lbb(k,n,a,b):   # log beta-binomial pmf, ignoring the constant C(n,k)
    return lbeta(k+a,n-k+b)-lbeta(a,b)
def lbin(k,n,p):
    if p<=0: return 0.0 if k==0 else -1e18
    if p>=1: return 0.0 if k==n else -1e18
    return k*math.log(p)+(n-k)*math.log(1-p)
SVALS=[1.0,2.0,3.0,5.0,10.0,1e6]   # 1e6 ~ binomial
res={s:collections.defaultdict(lambda:[0,0.0,0]) for s in SVALS}
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
        for st,idxs in byst.items():
            obs=[]; sref=salt=0
            for i in idxs:
                p=f[i].split(':')
                if len(p)<=ci: continue
                try: catg=[int(x) for x in p[ci].split(',')]
                except ValueError: continue
                rc,ac=catg[ri],catg[ai]; d=rc+ac
                if d<6: continue
                obs.append((d,ac)); sref+=rc; salt+=ac
            if len(obs)<6: continue
            tot=sref+salt
            if tot==0: continue
            q=salt/tot
            if min(q,1-q)<=0: continue
            for s in SVALS:
                for d,ac in obs:
                    b=binof(d)
                    if b is None: continue
                    lhet = lbin(ac,d,0.5) if s>1e5 else lbb(ac,d,s,s)
                    lhom = max(lbin(ac,d,E), lbin(ac,d,1-E))
                    ishet = 1 if lhet>lhom else 0
                    a=res[s][b]; a[0]+=ishet; a[1]+=2*q*(1-q); a[2]+=1
print('%-10s'%'Beta(s,s)', ' '.join('%8s'%n for n in NAMES), '%9s %8s'%('spread','mean F'))
print('-'*104)
for s in SVALS:
    row=[]; vals=[]
    for b in range(len(NAMES)):
        a=res[s].get(b)
        if not a or a[1]<=0 or a[2]<200: row.append('%8s'%'-'); continue
        F=1-a[0]/a[1]; row.append('%8.3f'%F); vals.append(F)
    sp=(max(vals)-min(vals)) if len(vals)>1 else float('nan')
    mn=sum(vals)/len(vals) if vals else float('nan')
    lab='binomial' if s>1e5 else 's=%.0f'%s
    print('%-10s'%lab, ' '.join(row), '%9.3f %8.3f'%(sp,mn))
print('\nFLATTEST row = the self-consistent model. Its value is the F this data supports.')
print('If no row is flat, a single overdispersion parameter cannot describe the data.')
