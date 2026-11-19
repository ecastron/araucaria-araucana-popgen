#!/usr/bin/env python3
"""F_IS vs read depth, fine resolution, with a matched HWE null.

Allele frequency is estimated ONCE per locus per stand from READ counts across
all individuals (independent of genotype calling and of depth). Observed
heterozygosity is then measured within each depth bin against that frequency.
Estimating q inside a narrow depth slice, as I did before, starves the bins.
"""
import sys, math, random, collections
IDX={'C':0,'A':1,'T':2,'G':3}
vcf,popmap,out=sys.argv[1],sys.argv[2],sys.argv[3]
rng=random.Random(42)
stand={}
for ln in open(popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]
EDGES=[6,8,10,12,14,16,18,20,23,26,30,35,40,50,60,75,100,150,250]
def binof(d):
    for j in range(len(EDGES)-1):
        if EDGES[j]<=d<EDGES[j+1]: return j
    return None
obs=collections.defaultdict(lambda:[0,0.0,0])
freqs=[]; depths=[]
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
        try: gi=fmt.index('GT'); ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        allref=0; allalt=0; dl=[]
        for st,idxs in byst.items():
            # q from READS over the whole locus within this stand
            sref=salt=0; g=[]
            for i in idxs:
                p=f[i].split(':')
                if len(p)<=ci: continue
                gt=p[gi]
                if gt in ('./.','.|.'): continue
                try: catg=[int(x) for x in p[ci].split(',')]
                except ValueError: continue
                rc,ac=catg[ri],catg[ai]; d=rc+ac
                if d<6: continue
                sref+=rc; salt+=ac; dl.append(d)
                g.append((d, 1 if gt in ('0/1','1/0','0|1','1|0') else 0))
            tot=sref+salt
            if tot==0 or len(g)<6: continue
            q=salt/tot
            if min(q,1-q)<=0: continue
            allref+=sref; allalt+=salt
            for d,ishet in g:
                b=binof(d)
                if b is None: continue
                a=obs[b]; a[0]+=ishet; a[1]+=2*q*(1-q); a[2]+=1
        if allref+allalt>0 and len(dl)>=20:
            freqs.append(allalt/(allref+allalt)); depths.extend(dl[:40])
E=0.001
def lnb(k,n,p):
    if p<=0: return -1e18 if k>0 else 0.0
    if p>=1: return -1e18 if k<n else 0.0
    return k*math.log(p)+(n-k)*math.log(1-p)
nul=collections.defaultdict(lambda:[0,0.0,0])
for q in freqs[:6000]:
    if q<=0 or q>=1: continue
    for _ in range(60):
        u=rng.random()
        g = 0 if u<(1-q)**2 else (2 if u>1-q*q else 1)
        d=depths[rng.randrange(len(depths))]
        pr = E if g==0 else (1-E if g==2 else 0.5)
        ac=sum(1 for _ in range(d) if rng.random()<pr)
        b=binof(d)
        if b is None: continue
        lhet=lnb(ac,d,0.5); lhom=max(lnb(ac,d,E),lnb(ac,d,1-E))
        call = 1 if 2*(lhet-lhom)>0 else 0
        a=nul[b]; a[0]+=call; a[1]+=2*q*(1-q); a[2]+=1
with open(out,'w') as fo:
    fo.write('lo\thi\tmid\tn_geno\tHo\tHe\tF_obs\tF_null\n')
    for b in range(len(EDGES)-1):
        a=obs.get(b)
        if not a or a[1]<=0 or a[2]<200: continue
        lo,hi=EDGES[b],EDGES[b+1]
        nb=nul.get(b)
        fn = (1-nb[0]/nb[1]) if (nb and nb[1]>0 and nb[2]>=200) else float('nan')
        fo.write('%d\t%d\t%.1f\t%d\t%.5f\t%.5f\t%.5f\t%.5f\n'%(
            lo,hi,(lo+hi)/2,a[2],a[0]/a[2],a[1]/a[2],1-a[0]/a[1],fn))
print(open(out).read())
