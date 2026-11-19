#!/usr/bin/env python3
"""Re-call genotypes from raw CATG counts at several stringencies.

Strict calling misses real heterozygotes at low depth (F too high).
Relaxed calling turns sequencing errors into fake ones (F too low).
A CORRECT rule gives the SAME F at every depth -- so tune until the depth
gradient flattens, then read off the value.

Rule: likelihood ratio of het vs best homozygote, using a binomial with
per-read error e. Call het when 2*ln(LR) exceeds a cutoff.
"""
import sys, math, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
BINS=[(6,10),(10,15),(15,20),(20,30),(30,50),(50,100),(100,10**9)]
NAMES=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
E=0.001
def lnb(k,n,p):
    if p<=0: return -1e18 if k>0 else 0.0
    if p>=1: return -1e18 if k<n else 0.0
    return k*math.log(p)+(n-k)*math.log(1-p)

CUTS=[-6.0,-2.0,0.0,2.0,6.0,12.0]   # 2lnLR cutoff; lower = more permissive
res={c:collections.defaultdict(lambda:[0,0,0,0.0]) for c in CUTS}  # cut->bin->[aa,ab,bb,hexp_acc]
cols=None
recs=[]
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
        obs=[]
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; d=rc+ac
            if d<6: continue
            k=next((j for j,(lo,hi) in enumerate(BINS) if lo<=d<hi),None)
            if k is None: continue
            obs.append((k,rc,ac,d))
        if len(obs)<20: continue
        recs.append(obs)
print('loci: %d'%len(recs))
for cut in CUTS:
    for obs in recs:
        per=collections.defaultdict(lambda:[0,0,0])
        for k,rc,ac,d in obs:
            lhet=lnb(ac,d,0.5)
            lhom=max(lnb(ac,d,E), lnb(ac,d,1-E))
            g = 1 if 2*(lhet-lhom)>cut else (0 if rc>=ac else 2)
            per[k][g]+=1
        for k,(aa,ab,bb) in per.items():
            n=aa+ab+bb
            if n<10: continue
            p_=(ab+2*bb)/(2*n)
            if min(p_,1-p_)<=0: continue
            a=res[cut][k]; a[0]+=aa; a[1]+=ab; a[2]+=bb; a[3]+=2*p_*(1-p_)*n
print()
print('%-9s'%'2lnLR cut', ' '.join('%8s'%n for n in NAMES), '   spread')
print('-'*88)
for cut in CUTS:
    row=[]; vals=[]
    for k in range(len(NAMES)):
        a=res[cut].get(k)
        if not a or a[3]<=0: row.append('%8s'%'-'); continue
        F=1-a[1]/a[3]; row.append('%8.3f'%F); vals.append(F)
    sp = (max(vals)-min(vals)) if len(vals)>1 else float('nan')
    print('%-9.1f'%cut, ' '.join(row), '   %6.3f'%sp)
print('\nThe rule whose row is FLATTEST across depth is the self-consistent one.')
print('Its common value is the best estimate of F that this data supports.')
