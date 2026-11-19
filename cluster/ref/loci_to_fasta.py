#!/usr/bin/env python3
"""One consensus sequence per ipyrad locus, named RAD_<n> to match the VCF CHROM.

.loci format: sample<TAB>SEQUENCE lines, each locus block terminated by a line
starting with '//'. Blocks are in the same order as RAD_1..RAD_N in the VCF.
Consensus = most common non-gap base per column (ties -> first seen).
"""
import sys, collections
inp,out=sys.argv[1],sys.argv[2]
n=0; kept=0
with open(inp) as fh, open(out,'w') as fo:
    seqs=[]
    for line in fh:
        if line.startswith('//'):
            n+=1
            if seqs:
                L=max(len(s) for s in seqs)
                cons=[]
                for i in range(L):
                    c=collections.Counter(s[i] for s in seqs if i<len(s) and s[i] in 'ACGT')
                    cons.append(c.most_common(1)[0][0] if c else 'N')
                seq=''.join(cons).replace('N','')
                if len(seq)>=60:
                    fo.write('>RAD_%d\n%s\n'%(n,seq)); kept+=1
            seqs=[]
            continue
        p=line.rstrip('\n').split()
        if len(p)>=2 and not line.startswith(' '):
            seqs.append(p[1].upper())
print('locus blocks: %d | written: %d'%(n,kept))
