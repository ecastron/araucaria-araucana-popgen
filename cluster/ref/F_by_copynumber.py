#!/usr/bin/env python3
"""F by depth, split by GENOMIC COPY NUMBER of each RAD locus.

Copy number comes from mapping each locus consensus against the 19.3 Gb reference
(>=80% query coverage, >=90% identity, distinct placements). It uses no genotype
information, so this is not circular.

If widespread structural variation drives the depth gradient, SINGLE-COPY loci
should be flat and near F=0, and multi-copy loci should carry the swing.
"""
import sys, math, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap,copies=sys.argv[1],sys.argv[2],sys.argv[3]
cn={}
for ln in open(copies):
    p=ln.split()
    if len(p)>=2: cn[p[0]]=int(p[1])
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
    return 1 if lnb(ac,d,0.5)>max(lnb(ac,d,E),lnb(ac,d,1-E)) else (0 if (d-ac)>=ac else 2)
EDG=[6,10,15,20,30,50,100,10**9]; NM=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
def dbin(d):
    for j in range(len(EDG)-1):
        if EDG[j]<=d<EDG[j+1]: return j
GRP=[('1 copy',lambda c:c==1),('2 copies',lambda c:c==2),
     ('3-5 copies',lambda c:3<=c<=5),('6+ copies',lambda c:c>=6),
     ('unmapped',lambda c:c==0)]
acc={g:collections.defaultdict(lambda:[0,0.0,0]) for g,_ in GRP}
nloc=collections.Counter()
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
        c=cn.get(f[0],0)
        grp=next((g for g,t in GRP if t(c)), None)
        if grp is None: continue
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        nloc[grp]+=1
        for st,idxs in byst.items():
            obs=[]; sref=salt=0
            for i in idxs:
                p=f[i].split(':')
                if len(p)<=ci: continue
                try: catg=[int(x) for x in p[ci].split(',')]
                except ValueError: continue
                rc,ac=catg[ri],catg[ai]; d=rc+ac
                if d<6: continue
                obs.append((d,rc,ac)); sref+=rc; salt+=ac
            if len(obs)<6: continue
            tot=sref+salt
            if tot==0: continue
            q=salt/tot
            if min(q,1-q)<=0: continue
            he1=2*q*(1-q)
            for d,rc,ac in obs:
                b=dbin(d)
                if b is None: continue
                a=acc[grp][b]; a[0]+=(call(ac,d)==1); a[1]+=he1; a[2]+=1
print('%-12s %8s  %s'%('locus class','SNPs',''.join('%8s'%n for n in NM)))
print('-'*88)
for g,_ in GRP:
    row=[]; vals=[]
    for b in range(len(NM)):
        a=acc[g].get(b)
        if not a or a[1]<=0 or a[2]<1000: row.append('%8s'%'-'); continue
        F=1-a[0]/a[1]; vals.append(F); row.append('%8.3f'%F)
    sp='%.3f'%(max(vals)-min(vals)) if len(vals)>1 else '-'
    print('%-12s %8d  %s   spread %s'%(g,nloc[g],''.join(row),sp))
