#!/usr/bin/env python3
"""H and F with a MAXIMUM depth cap at 1.5x and 2x mean global coverage.

Rationale: high-depth genotypes show heterozygote EXCESS (paralog-like), so
capping removes the collapsed-repeat tail while retaining the bulk of the data.
Reports Ho, He and F_IS computed WITHIN stands (Wahlund cannot contribute).
"""
import sys, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]

# pass 1: mean global coverage over all called genotypes
tot=0; n=0
cols=None
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in stand]; continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            d=catg[ri]+catg[ai]
            if d>0: tot+=d; n+=1
mean=tot/n
print('mean global coverage over %d called genotypes: %.2fx'%(n,mean))
CAPS=[('no cap',10**9),('1.5x mean',1.5*mean),('2x mean',2.0*mean)]
for nm,cap in CAPS: print('  %-10s -> max depth %.1f'%(nm,cap))

for nm,cap in CAPS:
    acc=collections.defaultdict(lambda:[0,0.0,0])   # stand -> hobs, hexp, ngeno
    cols=None
    with open(vcf) as fh:
        for line in fh:
            if line.startswith('##'): continue
            f=line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'):
                cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
                cols=[(i,s) for i,s in cols if s in stand]
                byst=collections.defaultdict(list)
                for i,s in cols: byst[stand[s]].append((i,s))
                continue
            if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
            fmt=f[8].split(':')
            try: gi=fmt.index('GT'); ci=fmt.index('CATG')
            except ValueError: continue
            ri,ai=IDX[f[3]],IDX[f[4]]
            for st,members in byst.items():
                aa=ab=bb=0
                for i,s in members:
                    p=f[i].split(':')
                    if len(p)<=ci: continue
                    gt=p[gi]
                    if gt in ('./.','.|.'): continue
                    try: catg=[int(x) for x in p[ci].split(',')]
                    except ValueError: continue
                    d=catg[ri]+catg[ai]
                    if d<6 or d>cap: continue          # <-- the cap
                    if gt in ('0/0','0|0'): aa+=1
                    elif gt in ('1/1','1|1'): bb+=1
                    else: ab+=1
                m=aa+ab+bb
                if m<6: continue
                q=(ab+2*bb)/(2*m)
                if min(q,1-q)<=0: continue
                a=acc[st]; a[0]+=ab; a[1]+=2*q*(1-q)*m; a[2]+=m
    print('\n================ %s ================'%nm)
    print('%-6s %12s %10s %10s %9s'%('stand','n_genotypes','Ho','He','F_IS'))
    print('-'*54)
    TO=TE=TN=0
    for st in sorted(acc):
        ho,he,ng=acc[st]
        if he<=0 or ng==0: continue
        TO+=ho; TE+=he; TN+=ng
        print('%-6s %12d %10.4f %10.4f %9.3f'%(st,ng,ho/ng,he/ng,1-ho/he))
    print('-'*54)
    if TE>0:
        print('%-6s %12d %10.4f %10.4f %9.3f'%('ALL',TN,TO/TN,TE/TN,1-TO/TE))
