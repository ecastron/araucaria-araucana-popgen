#!/usr/bin/env python3
"""Two filters, then H and F.

1. MAX depth cap at 1.5x / 2x mean global coverage (drops the paralog-like tail)
2. Drop SITES whose mean alt-allele balance across heterozygous genotypes falls
   outside 0.30-0.70 (drops loci whose 'heterozygotes' are not 1:1)

F_IS computed WITHIN stands so Wahlund cannot contribute.
NOTE: filter 2 selects on allele balance, so the balance/spread statistic would be
circular afterwards. F (heterozygote COUNT) is a different measurement and is not.
"""
import sys, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]

recs=[]; tot=0; n=0; cols=None
with open(vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in stand]; continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: gi=fmt.index('GT'); ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        row=[]
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            gt=p[gi]
            if gt in ('./.','.|.'): continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]; d=rc+ac
            if d<=0: continue
            tot+=d; n+=1
            g = 0 if gt in ('0/0','0|0') else (2 if gt in ('1/1','1|1') else 1)
            row.append((s,g,rc,ac,d))
        if row: recs.append(row)
mean=tot/n
print('mean global coverage: %.2fx  over %d genotypes, %d sites'%(mean,n,len(recs)))
print('caps: 1.5x = %.1f   2x = %.1f\n'%(1.5*mean,2*mean))

def run(cap, balfilter, label):
    acc=collections.defaultdict(lambda:[0,0.0,0])
    nsite=0; ndrop=0
    for row in recs:
        keep=[(s,g,rc,ac,d) for (s,g,rc,ac,d) in row if 6<=d<=cap]
        if not keep: continue
        if balfilter:
            hb=[ac/d for (s,g,rc,ac,d) in keep if g==1]
            if len(hb)>=3:
                mb=sum(hb)/len(hb)
                if not (0.30<=mb<=0.70): ndrop+=1; continue
        nsite+=1
        byst=collections.defaultdict(lambda:[0,0,0])
        for (s,g,rc,ac,d) in keep: byst[stand[s]][g]+=1
        for st,(aa,ab,bb) in byst.items():
            m=aa+ab+bb
            if m<6: continue
            q=(ab+2*bb)/(2*m)
            if min(q,1-q)<=0: continue
            a=acc[st]; a[0]+=ab; a[1]+=2*q*(1-q)*m; a[2]+=m
    TO=TE=TN=0
    for st,(ho,he,ng) in acc.items(): TO+=ho; TE+=he; TN+=ng
    if TE<=0: print('%-34s  no data'%label); return
    print('%-34s sites %6d (dropped %5d)  Ho %.4f  He %.4f  F_IS %.3f'%(
        label,nsite,ndrop,TO/TN,TE/TN,1-TO/TE))
    return acc

print('%-34s %s'%('filter','result'))
print('-'*104)
run(10**9, False, 'no cap, no balance filter')
run(1.5*mean, False, 'cap 1.5x mean')
run(2.0*mean, False, 'cap 2x mean')
run(10**9, True,  'balance filter only (0.30-0.70)')
a15=run(1.5*mean, True,  'cap 1.5x + balance filter')
a20=run(2.0*mean, True,  'cap 2x   + balance filter')
if a15:
    print('\n=== per stand, cap 1.5x + balance filter ===')
    print('%-6s %12s %9s %9s %8s'%('stand','n_genotypes','Ho','He','F_IS'))
    for st in sorted(a15):
        ho,he,ng=a15[st]
        if he<=0: continue
        print('%-6s %12d %9.4f %9.4f %8.3f'%(st,ng,ho/ng,he/ng,1-ho/he))
