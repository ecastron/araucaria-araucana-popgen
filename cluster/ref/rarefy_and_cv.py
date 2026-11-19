#!/usr/bin/env python3
"""Two tests.

(1) RAREFACTION: subsample every genotype's reads down to a fixed ceiling K.
    Equalises information without discarding genotypes. Then stratify F by the
    genotype's ORIGINAL depth. If the gradient vanishes, it was read count. If
    genotypes that were 100x still differ from natively-30x ones after both are
    cut to K, the difference is in the loci, not the depth.

(2) DEPTH-VARIANCE FILTER: drop loci whose depth varies wildly across individuals
    (CV of depth), the signature of segregating copy number.
"""
import sys, math, random, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap=sys.argv[1],sys.argv[2]
K=int(sys.argv[3]) if len(sys.argv)>3 else 30
rng=random.Random(42)
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
ORIG=[(20,30),(30,50),(50,100),(100,10**9)]
ONAME=['orig 20-30','orig 30-50','orig 50-100','orig 100+']
def obin(d):
    for j,(lo,hi) in enumerate(ORIG):
        if lo<=d<hi: return j
raw=collections.defaultdict(lambda:[0,0.0,0])
rar=collections.defaultdict(lambda:[0,0.0,0])
cvkeep=collections.defaultdict(lambda:[0,0.0,0])
cvdrop=collections.defaultdict(lambda:[0,0.0,0])
cols=None; ncv=0; ndrop=0
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
        # locus-level depth CV across individuals
        alld=[]
        for st,idxs in byst.items():
            for i in idxs:
                p=f[i].split(':')
                if len(p)<=ci: continue
                try: catg=[int(x) for x in p[ci].split(',')]
                except ValueError: continue
                d=catg[ri]+catg[ai]
                if d>=6: alld.append(d)
        if len(alld)<20: continue
        mu=sum(alld)/len(alld)
        sd=math.sqrt(sum((x-mu)**2 for x in alld)/len(alld))
        cv=sd/mu if mu>0 else 0
        ncv+=1
        highcv = cv>1.0
        if highcv: ndrop+=1
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
                b=obin(d)
                if b is None: continue
                g=call(ac,d)
                a=raw[b]; a[0]+= (g==1); a[1]+=he1; a[2]+=1
                tgt = cvdrop[b] if highcv else cvkeep[b]
                tgt[0]+= (g==1); tgt[1]+=he1; tgt[2]+=1
                # rarefy to K reads without replacement
                if d>K:
                    pool=[1]*ac+[0]*rc
                    s=rng.sample(pool,K); ac2=sum(s); d2=K
                else:
                    ac2,d2=ac,d
                g2=call(ac2,d2)
                r2=rar[b]; r2[0]+= (g2==1); r2[1]+=he1; r2[2]+=1
def show(title,acc):
    print('\n%s'%title)
    print('%-12s %10s %8s %9s'%('orig depth','n','Ho','F_IS'))
    vals=[]
    for b,nm in enumerate(ONAME):
        a=acc.get(b)
        if not a or a[1]<=0 or a[2]<500: continue
        F=1-a[0]/a[1]; vals.append(F)
        print('%-12s %10d %8.4f %9.4f'%(nm,a[2],a[0]/a[2],F))
    if len(vals)>1: print('  spread: %.3f'%(max(vals)-min(vals)))
print('loci evaluated: %d | high depth-CV (>1.0) dropped: %d (%.1f%%)'%(ncv,ndrop,100*ndrop/ncv))
show('=== (0) NO FILTER, native depth ===', raw)
show('=== (1) RAREFIED to %dx, stratified by ORIGINAL depth ==='%K, rar)
show('=== (2) DEPTH-CV FILTER: loci KEPT (CV<=1.0) ===', cvkeep)
show('=== (2) DEPTH-CV FILTER: loci DROPPED (CV>1.0) ===', cvdrop)
