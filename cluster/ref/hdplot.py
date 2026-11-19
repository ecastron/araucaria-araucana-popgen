#!/usr/bin/env python3
"""HDplot paralog detection (McKinney et al. 2017 Mol Ecol Resour).

Per SNP:
  H = fraction of genotyped individuals called heterozygous
  D = (A-B)/sqrt((A+B)/2), the read-ratio deviation summed over het calls
Collapsed paralogs: H high (every individual 'het' because two loci are stacked)
and |D| large (read ratio not 1:1). Real SNPs: H follows HWE, D ~ 0.

Reference-free and assembly-agnostic, so the collapsed haploid genome is irrelevant.
Also reports whether the FLAT-BALANCE genotypes are concentrated in paralog-like loci,
which is the actual question.
"""
import sys, math, argparse
ap=argparse.ArgumentParser()
ap.add_argument('--vcf',required=True); ap.add_argument('--minind',type=int,default=30)
ap.add_argument('--out',required=True)
a=ap.parse_args()
IDX={'C':0,'A':1,'T':2,'G':3}
out=open(a.out,'w'); out.write('chrom\tpos\tn\tH\tD\tflatfrac\n')
nloc=0; hist={}
import collections
Hb=collections.Counter(); para=0; tot=0; flat_in_para=0; flat_tot=0
with open(a.vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            ncol=len(line.split('\t')); continue
        f=line.rstrip('\n').split('\t')
        ref,alt=f[3],f[4]
        if len(alt)!=1 or ref not in IDX or alt not in IDX: continue
        ri,ai=IDX[ref],IDX[alt]
        fmt=f[8].split(':')
        try: gi=fmt.index('GT'); ci=fmt.index('CATG')
        except ValueError: continue
        n=0; nhet=0; A=0; B=0; flat=0; withdp=0
        for cell in f[9:]:
            p=cell.split(':')
            if len(p)<=ci: continue
            gt=p[gi]
            if gt in ('./.','.|.'): continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; dp=rc+ac
            if dp==0: continue
            n+=1
            if gt in ('0/1','1/0','0|1','1|0'):
                nhet+=1; A+=rc; B+=ac
            if dp>=20:
                withdp+=1
                if 0.15<ac/dp<0.85: flat+=1
        if n<a.minind: continue
        H=nhet/n
        D=(A-B)/math.sqrt((A+B)/2) if (A+B)>0 else 0.0
        ff=flat/withdp if withdp else 0.0
        nloc+=1; tot+=1
        Hb[round(H,1)]+=1
        isp = (H>0.55 and abs(D)>7)
        if isp: para+=1
        flat_tot+=flat
        if isp: flat_in_para+=flat
        out.write('%s\t%s\t%d\t%.4f\t%.3f\t%.4f\n'%(f[0],f[1],n,H,D,ff))
out.close()
print('SNPs evaluated (>=%d individuals): %d'%(a.minind,nloc))
print('\n=== distribution of H (fraction of individuals called heterozygous) ===')
for k in sorted(Hb): print('  H=%.1f  %8d  (%5.2f%%)'%(k,Hb[k],100*Hb[k]/nloc))
print('\n=== paralog-like loci (H>0.55 AND |D|>7, McKinney criterion) ===')
print('  %d / %d = %.2f%% of SNPs'%(para,tot,100*para/tot))
print('\n=== are the FLAT-BALANCE genotypes concentrated in those loci? ===')
print('  flat-balance genotypes total          : %d'%flat_tot)
print('  of which in paralog-like loci         : %d (%.2f%%)'%(flat_in_para,100*flat_in_para/flat_tot if flat_tot else 0))
print('  paralog-like loci are %.2f%% of all loci'%(100*para/tot))
print('\n  If flat-balance genotypes were CAUSED by paralogs, the second')
print('  percentage should be far larger than the third.')
