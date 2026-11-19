#!/usr/bin/env python3
"""Between-individual distance from READ-BASED allele dosage.

Genotype calls are unreliable in this dataset (documented), so distance is built
from the ALT read fraction per individual per locus -- a continuous dosage that
never passes through a genotype call. Manhattan distance over shared loci,
normalised by the number of loci compared (missing data differs between pairs).
"""
import sys, collections, math
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap,out=sys.argv[1],sys.argv[2],sys.argv[3]
keep={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: keep[p[0]]=p[1]
cols=None; dos=None; names=None
NL=0
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cc=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cc=[(i,s) for i,s in cc if s in keep]
            cols=cc; names=[s for _,s in cc]
            n=len(cols)
            S=[[0.0]*n for _ in range(n)]      # summed |d_i - d_j|
            C=[[0]*n   for _ in range(n)]      # count of shared loci
            continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        d=[None]*len(cols)
        for k,(i,s) in enumerate(cols):
            p=f[i].split(':')
            if len(p)<=ci: continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; t=rc+ac
            if t>=8: d[k]=ac/t
        NL+=1
        idxs=[k for k in range(len(cols)) if d[k] is not None]
        for a_ in range(len(idxs)):
            ka=idxs[a_]; da=d[ka]
            row=S[ka]; crow=C[ka]
            for b_ in range(a_+1,len(idxs)):
                kb=idxs[b_]
                row[kb]+=abs(da-d[kb]); crow[kb]+=1
n=len(cols)
with open(out,'w') as fo:
    fo.write('\t'+'\t'.join(names)+'\n')
    for i in range(n):
        vals=[]
        for j in range(n):
            if i==j: vals.append('0')
            else:
                a,b=(i,j) if i<j else (j,i)
                vals.append('%.6f'%(S[a][b]/C[a][b]) if C[a][b]>=500 else 'NA')
        fo.write(names[i]+'\t'+'\t'.join(vals)+'\n')
print('loci scanned %d | individuals %d'%(NL,n))
