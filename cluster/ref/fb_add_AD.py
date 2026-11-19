#!/usr/bin/env python3
"""FreeBayes -> polyRAD input: keep biallelic SNPs, add an AD field from RO/AO.

FreeBayes emits GT:DP:RO:QR:AO:QA:GL. polyRAD's VCF2RADdata needs AD (Number=R).
The ##FORMAT header line must be declared or VariantAnnotation drops the field.
"""
import sys
inp,out=sys.argv[1],sys.argv[2]
HDR='##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n'
nin=nout=0; wrote=False; last=False
with open(inp) as fh, open(out,'w') as fo:
    for line in fh:
        if line.startswith('##'):
            if line.startswith('##FORMAT'): fo.write(line); last=True; continue
            if last and not wrote: fo.write(HDR); wrote=True
            fo.write(line); continue
        if line.startswith('#CHROM'):
            if not wrote: fo.write(HDR); wrote=True
            fo.write(line); continue
        f=line.rstrip('\n').split('\t')
        nin+=1
        if len(f[3])!=1 or len(f[4])!=1 or f[4]=='.': continue      # biallelic SNPs only
        if f[3] not in 'ACGT' or f[4] not in 'ACGT': continue
        fmt=f[8].split(':')
        try: ro=fmt.index('RO'); ao=fmt.index('AO')
        except ValueError: continue
        cells=[]
        ok=True
        for c in f[9:]:
            p=c.split(':')
            # a no-data cell is "." -- appending AD to it would put "0,0" in the DP
            # slot and lose AD entirely. VCF allows trailing fields to be dropped,
            # so leave fully-missing genotypes as a bare ".".
            if c=='.' or len(p)<=max(ro,ao): cells.append('.'); continue
            try:
                r=int(p[ro]) if p[ro] not in ('.','') else 0
                a=int(p[ao].split(',')[0]) if p[ao] not in ('.','') else 0
            except ValueError:
                r=a=0
            cells.append(c+':%d,%d'%(r,a))
        f[8]=f[8]+':AD'
        fo.write('\t'.join(f[:9]+cells)+'\n'); nout+=1
print('records in: %d | biallelic SNPs written with AD: %d'%(nin,nout))
