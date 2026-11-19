#!/usr/bin/env python3
"""Add a proper AD field to the ipyrad VCF so polyRAD's VCF2RADdata() can read it.

ipyrad writes GT:DP:CATG, where CATG is per-base counts in C,A,T,G order.
polyRAD needs AD with Number=R -- one integer per allele, REF first then each ALT.

CRITICAL: the ##FORMAT header line must be declared, otherwise VariantAnnotation
parses the VCF without the field and VCF2RADdata() silently sees no read depth.
"""
import sys
IDX={'C':0,'A':1,'T':2,'G':3}
inp,out=sys.argv[1],sys.argv[2]
HDR='##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n'
n_var=0; n_gt=0; n_skip=0; wrote_hdr=False
with open(inp) as fh, open(out,'w') as fo:
    for line in fh:
        if line.startswith('##'):
            # insert the AD declaration immediately after the last FORMAT line
            if line.startswith('##FORMAT') :
                fo.write(line); last_format=True; continue
            if not wrote_hdr and locals().get('last_format'):
                fo.write(HDR); wrote_hdr=True
            fo.write(line); continue
        if line.startswith('#CHROM'):
            if not wrote_hdr: fo.write(HDR); wrote_hdr=True
            fo.write(line); continue
        f=line.rstrip('\n').split('\t')
        ref=f[3]; alts=f[4].split(',')
        alleles=[ref]+alts
        fmt=f[8].split(':')
        if 'CATG' not in fmt or any(a not in IDX for a in alleles):
            n_skip+=1; fo.write(line); continue
        ci=fmt.index('CATG')
        pos=[IDX[a] for a in alleles]
        newfmt=fmt+['AD'] if 'AD' not in fmt else fmt
        cells=[]
        for cell in f[9:]:
            p=cell.split(':')
            if len(p)<=ci:
                cells.append(cell+':'+','.join(['0']*len(alleles))); continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError:
                cells.append(cell+':'+','.join(['0']*len(alleles))); continue
            ad=[catg[j] for j in pos]
            cells.append(cell+':'+','.join(str(x) for x in ad))
            n_gt+=1
        f[8]=':'.join(newfmt)
        fo.write('\t'.join(f[:9]+cells)+'\n')
        n_var+=1
print('variants written with AD : %d'%n_var)
print('variants skipped (indel/non-ACGT/no CATG): %d'%n_skip)
print('genotype cells given AD  : %d'%n_gt)
