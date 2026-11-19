#!/usr/bin/env python3
"""Estimate allele frequency from READS instead of from genotype calls.

p_GT   = (2*nAA + nAB)/(2n)          -- derived from the calls we suspect are wrong
p_read = mean over individuals of (alt reads / total reads)  -- каждый individual
         weighted equally, independent of how the caller resolved the genotype

If undercalling is converting heterozygotes into homozygotes at random, p_GT stays
roughly unbiased but noisy. If it is biased toward one allele (e.g. reference bias),
p_GT and p_read diverge systematically. Either way H_exp should be built from p_read.
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
agg=collections.defaultdict(lambda:[0.0,0.0,0.0,0,0.0,0.0])  # hobs, hexpGT, hexpRD, n, sum pGT, sum pRD
cols=None
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
        try: gi=fmt.index('GT'); ci=fmt.index('CATG')
        except ValueError: continue
        byb=collections.defaultdict(lambda:[0,0,0,[]])   # bin -> aa,ab,bb,[read dosage]
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            gt=p[gi]
            if gt in ('./.','.|.'): continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; d=rc+ac
            if d<6: continue
            k=next((j for j,(lo,hi) in enumerate(BINS) if lo<=d<hi),None)
            if k is None: continue
            b=byb[k]
            g = 0 if gt in ('0/0','0|0') else (2 if gt in ('1/1','1|1') else 1)
            b[g]+=1; b[3].append(ac/d)
        for k,(aa,ab,bb,dos) in byb.items():
            n=aa+ab+bb
            if n<20: continue
            pGT=(ab+2*bb)/(2*n)                 # ALT frequency from calls
            pRD=sum(dos)/len(dos)               # ALT frequency from reads
            if min(pGT,1-pGT)<=0 and min(pRD,1-pRD)<=0: continue
            a=agg[k]
            a[0]+=ab
            a[1]+=2*pGT*(1-pGT)*n
            a[2]+=2*pRD*(1-pRD)*n
            a[3]+=1; a[4]+=pGT; a[5]+=pRD
print('%-8s %8s %10s %10s %12s %12s'%('depth','n_loci','mean p_GT','mean p_read','F (p from GT)','F (p from reads)'))
print('-'*72)
for k,nm in enumerate(NAMES):
    if k not in agg: continue
    ho,heG,heR,nl,sG,sR=agg[k]
    if heG<=0 or heR<=0: continue
    print('%-8s %8d %10.4f %10.4f %12.3f %12.3f'%(nm,nl,sG/nl,sR/nl,1-ho/heG,1-ho/heR))
print('\np_read is independent of genotype calling. If the two F columns diverge,')
print('the genotype-based allele frequency was itself corrupted by miscalling.')
