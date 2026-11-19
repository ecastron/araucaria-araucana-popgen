#!/usr/bin/env python3
"""De novo control: depth-stratified allele balance from the ipyrad cr30 VCF.

Same statistic and same depth bins as balance_from_pileup.py, on the SAME
libraries, so the reference-mapped arm has a like-for-like comparison instead
of being held against the old single-threshold (>=30x) number.

Two views:
  ALL  - every genotype at a variant site
  HET  - only genotypes ipyrad called 0/1 (the original framing)
"""
import sys, math, argparse
ap = argparse.ArgumentParser()
ap.add_argument('vcf'); ap.add_argument('--libs', required=True)
a = ap.parse_args()
want = [l.split()[0] for l in open(a.libs) if l.strip()]
stand = {l.split()[0]: l.split()[1] for l in open(a.libs) if l.strip()}

BINS=[(6,10),(10,15),(15,20),(20,30),(30,50),(50,100),(100,10**9)]
NAMES=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
IDX={'C':0,'A':1,'T':2,'G':3}

def newacc(): return dict(n=[0]*7, hn=[0]*7, hs=[0.0]*7, hq=[0.0]*7, hd=[0.0]*7, c=[0]*7)
acc={s:{'ALL':newacc(),'HET':newacc()} for s in want}

col={}
with open(a.vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            # ipyrad appends a '.1' assembly suffix to sample names, so the VCF
            # carries SRR11841620.1 rather than the bare accession.
            for i,s in enumerate(f[9:],9):
                base=s.split('.')[0]
                if base in acc: col[base]=i
            miss=[s for s in want if s not in col]
            if miss: sys.exit('libraries not in VCF: %s (VCF has e.g. %s)'%(miss,f[9:12]))
            continue
        ref,alt=f[3],f[4]
        if len(alt)!=1 or ref not in IDX or alt not in IDX: continue
        ri,ai=IDX[ref],IDX[alt]
        fmt=f[8].split(':')
        try: gi,ci=fmt.index('GT'),fmt.index('CATG')
        except ValueError: continue
        for s,i in col.items():
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
            bal=ac/dp
            views=['ALL']+(['HET'] if gt in ('0/1','1/0','0|1','1|0') else [])
            for v in views:
                A=acc[s][v]; A['n'][k]+=1
                if 0.15<bal<0.85:
                    A['hn'][k]+=1; A['hs'][k]+=bal; A['hq'][k]+=bal*bal; A['hd'][k]+=dp
                if 0.45<=bal<0.55: A['c'][k]+=1

for s in want:
    for v in ('ALL','HET'):
        A=acc[s][v]
        print('\n=== %s (%s)  DE NOVO ipyrad  view=%s ==='%(s,stand[s],v))
        print('%-7s %9s %9s %7s %7s %7s %8s %8s %7s'%(
            'depth','sites','het band','het%','mean','SD','binom SD','overdisp','c0.5%'))
        print('-'*80)
        for k,nm in enumerate(NAMES):
            if not A['n'][k]: continue
            row='%-7s %9d %9d %6.2f%%'%(nm,A['n'][k],A['hn'][k],100*A['hn'][k]/A['n'][k])
            if A['hn'][k]>=50:
                m=A['hs'][k]/A['hn'][k]
                sd=math.sqrt(max(A['hq'][k]/A['hn'][k]-m*m,0))
                bsd=math.sqrt(0.25/(A['hd'][k]/A['hn'][k]))
                row+=' %7.4f %7.4f %8.4f %7.2fx %6.2f%%'%(
                    m,sd,bsd,sd/bsd,100*A['c'][k]/A['hn'][k])
            print(row)
