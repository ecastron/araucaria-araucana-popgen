#!/usr/bin/env python3
"""Are high-depth loci a fair sample, or is F~0 an artifact of WHICH loci get deep?

PAIRED design: use only loci that carry BOTH shallow (6-15x) and deep (>=50x)
genotypes, and compute F separately from each subset OF THE SAME LOCI. Locus
identity is then held constant -- any difference is depth, not locus choice.

Also reports whether deep loci differ from the rest in allele frequency and in
how many loci contribute at all.
"""
import sys, math, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
SH=(6,16); DP=(50,10**9)
nloc_any=0; nloc_deep=0; nloc_paired=0
sh_ho=sh_he=dp_ho=dp_he=0.0
maf_deep=[]; maf_all=[]
ndeep_geno=0
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
        sh=[0,0,0]; dp=[0,0,0]
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            gt=p[gi]
            if gt in ('./.','.|.'): continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; d=rc+ac
            g = 0 if gt in ('0/0','0|0') else (2 if gt in ('1/1','1|1') else 1)
            if SH[0]<=d<SH[1]: sh[g]+=1
            elif d>=DP[0]:     dp[g]+=1; ndeep_geno+=1
        nloc_any+=1
        nS=sum(sh); nD=sum(dp)
        if nD>=10: nloc_deep+=1
        if nS<10 or nD<10: continue
        nloc_paired+=1
        for arr,n,acc in ((sh,nS,'s'),(dp,nD,'d')):
            p_=(2*arr[0]+arr[1])/(2*n); q=1-p_
            if min(p_,q)<=0: continue
            he=2*p_*q*n
            if acc=='s': sh_ho+=arr[1]; sh_he+=he
            else:        dp_ho+=arr[1]; dp_he+=he
        p_=(2*dp[0]+dp[1])/(2*nD); maf_deep.append(min(p_,1-p_))
print('loci in VCF (biallelic)              : %d'%nloc_any)
print('loci with >=10 deep (>=50x) genotypes: %d  (%.1f%%)'%(nloc_deep,100*nloc_deep/nloc_any))
print('loci with >=10 BOTH shallow and deep : %d  (%.1f%%)  <- the paired set'%(nloc_paired,100*nloc_paired/nloc_any))
print('deep genotypes total                 : %d'%ndeep_geno)
print()
print('=== PAIRED COMPARISON on the SAME %d loci ==='%nloc_paired)
if sh_he>0 and dp_he>0:
    Fs=1-sh_ho/sh_he; Fd=1-dp_ho/dp_he
    print('  shallow genotypes (6-15x) : obs %.0f / exp %.0f  ->  F_IS = %.3f'%(sh_ho,sh_he,Fs))
    print('  deep genotypes  (>=50x)   : obs %.0f / exp %.0f  ->  F_IS = %.3f'%(dp_ho,dp_he,Fd))
    print('  difference                : %.3f'%(Fs-Fd))
    print()
    print('  Same loci in both rows. A large drop means DEPTH, not locus identity.')
if maf_deep:
    maf_deep.sort()
    print('\n=== allele frequency of the deep-locus set ===')
    print('  median MAF %.3f   quartiles %.3f / %.3f'%(
        maf_deep[len(maf_deep)//2], maf_deep[len(maf_deep)//4], maf_deep[3*len(maf_deep)//4]))
