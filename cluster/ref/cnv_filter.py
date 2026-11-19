#!/usr/bin/env python3
"""Segregating copy-number filter.

A duplication present in SOME trees gives those trees ~2x the depth AT THAT LOCUS
plus paralogous variants that read as heterozygous. The diagnostic is therefore a
genotype's depth relative to ITS OWN locus median -- not absolute depth, and not
the locus average (which is diluted because only carriers are deep).

Test A: is F better explained by the depth RATIO than by absolute depth?
Test B: does removing high-ratio genotypes flatten the depth gradient?
"""
import sys, math, collections, statistics as stt
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
E=0.001
def lnb(k,n,p):
    if p<=0: return 0.0 if k==0 else -1e18
    if p>=1: return 0.0 if k==n else -1e18
    return k*math.log(p)+(n-k)*math.log(1-p)
def call(ac,d):
    lhet=lnb(ac,d,0.5); lhom=max(lnb(ac,d,E),lnb(ac,d,1-E))
    return 1 if lhet>lhom else (0 if (d-ac)>=ac else 2)
DEDG=[6,10,15,20,30,50,100,10**9]; DNM=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
REDG=[0,0.6,0.8,1.0,1.25,1.5,2.0,3.0,10**9]
RNM=['<0.6','0.6-0.8','0.8-1.0','1.0-1.25','1.25-1.5','1.5-2','2-3','3+']
def ib(v,e):
    for j in range(len(e)-1):
        if e[j]<=v<e[j+1]: return j
byratio=collections.defaultdict(lambda:[0,0.0,0])
bydepth=collections.defaultdict(lambda:[0,0.0,0])
filt={t:collections.defaultdict(lambda:[0,0.0,0]) for t in (1.5,2.0,3.0)}
nflag=collections.Counter(); ntot=0
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
            continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        alld=[]; per=[]
        for st,idxs in byst.items():
            obs=[]; sref=salt=0
            for i in idxs:
                p=f[i].split(':')
                if len(p)<=ci: continue
                try: catg=[int(x) for x in p[ci].split(',')]
                except ValueError: continue
                rc,ac=catg[ri],catg[ai]; d=rc+ac
                if d<6: continue
                obs.append((d,rc,ac)); sref+=rc; salt+=ac; alld.append(d)
            if len(obs)<6: continue
            tot=sref+salt
            if tot==0: continue
            q=salt/tot
            if min(q,1-q)<=0: continue
            per.append((2*q*(1-q),obs))
        if len(alld)<20 or not per: continue
        lmed=stt.median(alld)
        if lmed<=0: continue
        for he1,obs in per:
            for d,rc,ac in obs:
                g=call(ac,d); r=d/lmed; ntot+=1
                rb=ib(r,REDG); db=ib(d,DEDG)
                if rb is not None:
                    a=byratio[rb]; a[0]+=(g==1); a[1]+=he1; a[2]+=1
                if db is not None:
                    a=bydepth[db]; a[0]+=(g==1); a[1]+=he1; a[2]+=1
                    for t in filt:
                        if r<=t:
                            a2=filt[t][db]; a2[0]+=(g==1); a2[1]+=he1; a2[2]+=1
                        else: nflag[t]+=1
def show(nm,acc,names,minn=2000):
    vals=[]; rows=[]
    for b in range(len(names)):
        a=acc.get(b)
        if not a or a[1]<=0 or a[2]<minn: rows.append((names[b],None,None,None)); continue
        F=1-a[0]/a[1]; vals.append(F); rows.append((names[b],a[2],a[0]/a[2],F))
    print('\n%s'%nm)
    print('  %-10s %10s %8s %9s'%('bin','n','Ho','F_IS'))
    for n_,c,h,F in rows:
        if c is None: continue
        print('  %-10s %10d %8.4f %9.4f'%(n_,c,h,F))
    if len(vals)>1: print('  spread: %.3f'%(max(vals)-min(vals)))
    return vals
print('genotypes: %d'%ntot)
show('=== A: F by DEPTH RATIO (genotype depth / its locus median) ===', byratio, RNM)
show('=== reference: F by ABSOLUTE depth ===', bydepth, DNM)
for t in (3.0,2.0,1.5):
    show('=== B: F by absolute depth, keeping only ratio <= %.1f  (%d genotypes dropped, %.1f%%) ==='
         %(t,nflag[t],100*nflag[t]/ntot), filt[t], DNM)
