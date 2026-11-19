#!/usr/bin/env python3
"""Observed vs HWE-EXPECTED heterozygotes, stratified by allele frequency.

Raw counts of heterozygotes are uninformative because rare variants legitimately
have almost none. The question is whether each locus has FEWER heterozygotes than
its own allele frequency predicts, i.e. per-locus F = 1 - Hobs/Hexp.
"""
import sys, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf=sys.argv[1]
bins=[(0.0,0.02),(0.02,0.05),(0.05,0.10),(0.10,0.20),(0.20,0.30),(0.30,0.40),(0.40,0.501)]
names=['0-0.02','0.02-0.05','0.05-0.10','0.10-0.20','0.20-0.30','0.30-0.40','0.40-0.50']
acc={i:[0,0.0,0.0,0] for i in range(len(bins))}   # nloci, sumHobs, sumHexp, ninds
with open(vcf) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        ref,alt=f[3],f[4]
        if len(alt)!=1 or ref not in IDX or alt not in IDX: continue
        fmt=f[8].split(':')
        try: gi=fmt.index('GT')
        except ValueError: continue
        naa=nab=nbb=0
        for cell in f[9:]:
            gt=cell.split(':')[gi]
            if gt in ('0/0','0|0'): naa+=1
            elif gt in ('0/1','1/0','0|1','1|0'): nab+=1
            elif gt in ('1/1','1|1'): nbb+=1
        n=naa+nab+nbb
        if n<30: continue
        p=(2*naa+nab)/(2*n); q=1-p
        maf=min(p,q)
        if maf<=0: continue
        hexp=2*p*q*n; hobs=nab
        for i,(lo,hi) in enumerate(bins):
            if lo<=maf<hi:
                a=acc[i]; a[0]+=1; a[1]+=hobs; a[2]+=hexp; a[3]+=n
                break
print('%-12s %9s %12s %12s %10s %8s'%('MAF bin','n_loci','obs het','exp het (HWE)','obs/exp','F'))
print('-'*70)
TO=TE=0
for i,nm in enumerate(names):
    nl,ho,he,ni=acc[i]
    if nl<10: continue
    TO+=ho; TE+=he
    r=ho/he if he else float('nan')
    print('%-12s %9d %12.0f %12.0f %10.3f %8.3f'%(nm,nl,ho,he,r,1-r))
print('-'*70)
print('%-12s %9s %12.0f %12.0f %10.3f %8.3f'%('ALL','',TO,TE,TO/TE,1-TO/TE))
print('\nF = 1 - obs/exp.  F~0 means HWE.  F=0.5 means half the expected heterozygotes are missing.')
print('Araucaria araucana is DIOECIOUS, so true inbreeding cannot produce a large positive F.')
