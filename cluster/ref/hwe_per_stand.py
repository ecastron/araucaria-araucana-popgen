#!/usr/bin/env python3
"""Same HWE test, computed WITHIN each stand.

The pooled version confounds the deficit with the Wahlund effect: pooling
differentiated stands creates apparent heterozygote deficit, worst at high MAF.
Within a single stand there is no subdivision, so Wahlund cannot contribute.
Allele frequencies are recomputed per stand.
"""
import sys, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
acc=collections.defaultdict(lambda:[0,0.0,0.0])
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
            print('stands: %d, libraries: %d'%(len(byst),len(cols)))
            continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: gi=fmt.index('GT')
        except ValueError: continue
        for st,idxs in byst.items():
            naa=nab=nbb=0
            for i in idxs:
                gt=f[i].split(':')[gi]
                if gt in ('0/0','0|0'): naa+=1
                elif gt in ('0/1','1/0','0|1','1|0'): nab+=1
                elif gt in ('1/1','1|1'): nbb+=1
            n=naa+nab+nbb
            if n<6: continue
            p=(2*naa+nab)/(2*n); q=1-p
            if min(p,q)<=0: continue
            a=acc[st]; a[0]+=1; a[1]+=nab; a[2]+=2*p*q*n
print()
print('%-6s %10s %10s %12s %10s %8s'%('stand','n_loci','obs het','exp het','obs/exp','F_IS'))
print('-'*62)
TO=TE=0
for st in sorted(acc):
    nl,ho,he=acc[st]
    if he<=0: continue
    TO+=ho; TE+=he
    print('%-6s %10d %10.0f %12.0f %10.3f %8.3f'%(st,nl,ho,he,ho/he,1-ho/he))
print('-'*62)
print('%-6s %10s %10.0f %12.0f %10.3f %8.3f'%('MEAN','',TO,TE,TO/TE,1-TO/TE))
print('\nWithin a stand there is no subdivision, so Wahlund cannot contribute.')
print('Whatever F_IS remains here is NOT population structure.')
