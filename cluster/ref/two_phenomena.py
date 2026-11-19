#!/usr/bin/env python3
"""Are the heterozygote DEFICIT and the flat allele BALANCE one phenomenon or two?

Both computed per GENOTYPE-depth stratum on the same data, within stands
(so Wahlund cannot contribute to F).

  F_IS   = 1 - Hobs/Hexp   -- are heterozygotes missing?
  spread = SD(balance) / sqrt(0.25/depth) among balance 0.15-0.85
           -- do the heterozygotes we DO see look real?

If one phenomenon, both should move the same way with depth.
Undercalling predicts F falls with depth. It says nothing about spread rising.
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
# per depth bin, per stand: genotype counts for F ; and balance moments for spread
gc=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0,0]))  # bin->stand->[aa,ab,bb]
bal=collections.defaultdict(lambda:[0,0.0,0.0,0.0])                          # bin->[n,sum,sumsq,sumdepth]
cols=None
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in stand]; continue
        ref,alt=f[3],f[4]
        if len(alt)!=1 or ref not in IDX or alt not in IDX: continue
        ri,ai=IDX[ref],IDX[alt]
        fmt=f[8].split(':')
        try: gi=fmt.index('GT'); ci=fmt.index('CATG')
        except ValueError: continue
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            gt=p[gi]
            if gt in ('./.','.|.'): continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; dp=rc+ac
            if dp<6: continue
            k=next((j for j,(lo,hi) in enumerate(BINS) if lo<=dp<hi),None)
            if k is None: continue
            g=gc[k][stand[s]]
            if gt in ('0/0','0|0'): g[0]+=1
            elif gt in ('1/1','1|1'): g[2]+=1
            else: g[1]+=1
            b=ac/dp
            if 0.15<b<0.85:
                r=bal[k]; r[0]+=1; r[1]+=b; r[2]+=b*b; r[3]+=dp
print('%-8s %10s %10s %8s %12s %10s'%('depth','n_genos','F_IS(within)','','n_balance','spread'))
print('-'*66)
for k,nm in enumerate(NAMES):
    if k not in gc: continue
    ho=he=0.0; tot=0
    for st,(aa,ab,bb) in gc[k].items():
        n=aa+ab+bb
        if n<20: continue
        p_=(2*aa+ab)/(2*n); q=1-p_
        if min(p_,q)<=0: continue
        ho+=ab; he+=2*p_*q*n; tot+=n
    F = 1-ho/he if he>0 else float('nan')
    r=bal[k]
    if r[0]>=50:
        m=r[1]/r[0]; sd=math.sqrt(max(r[2]/r[0]-m*m,0)); dbar=r[3]/r[0]
        sp=sd/math.sqrt(0.25/dbar)
    else: sp=float('nan')
    print('%-8s %10d %10.3f %8s %12d %10.2f'%(nm,tot,F,'',r[0],sp))
print('\nF_IS falling with depth  -> heterozygote UNDERCALLING (a depth artifact)')
print('spread rising with depth -> a SEPARATE phenomenon; undercalling does not predict it')
