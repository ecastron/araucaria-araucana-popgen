#!/usr/bin/env python3
"""Do multi-copy loci change the paper's structure inferences?

Hudson FST (ratio of averages, as used in the manuscript) per stand pair, computed on
  (a) ALL loci   (b) verified SINGLE-COPY loci only  (c) multi-copy only
Copy number from mapping locus consensus sequences against the reference genome.
Also reports the among-stand variance component (AMOVA-like) for each set.
"""
import sys, math, collections, itertools
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
SETS={'all':lambda c:True,'single-copy':lambda c:c==1,'multi-copy':lambda c:c>=2}
# accumulate Hudson numerator/denominator per stand pair per set
num=collections.defaultdict(float); den=collections.defaultdict(float)
# and per-set global heterozygosity partition
HT=collections.defaultdict(float); HS=collections.defaultdict(float); NL=collections.Counter()
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
            STANDS=sorted(byst)
            continue
        if len(f[4])!=1 or f[3] not in IDX or f[4] not in IDX: continue
        c=cn.get(f[0],0)
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        ri,ai=IDX[f[3]],IDX[f[4]]
        # allele counts per stand from READS (independent of genotype calling)
        p_=dict(); n_=dict()
        for st,idxs in byst.items():
            a=b=0
            for i in idxs:
                q=f[i].split(':')
                if len(q)<=ci: continue
                try: catg=[int(x) for x in q[ci].split(',')]
                except ValueError: continue
                a+=catg[ri]; b+=catg[ai]
            if a+b>=20: p_[st]=b/(a+b); n_[st]=a+b
        if len(p_)<8: continue
        for name,test in SETS.items():
            if not test(c): continue
            NL[name]+=1
            for s1,s2 in itertools.combinations(sorted(p_),2):
                x,y=p_[s1],p_[s2]
                nx,ny=n_[s1],n_[s2]
                # Hudson 1992 / Bhatia 2013 numerator and denominator
                nn=(x-y)**2 - x*(1-x)/(nx-1) - y*(1-y)/(ny-1)
                dd=x*(1-y)+y*(1-x)
                num[(name,s1,s2)]+=nn; den[(name,s1,s2)]+=dd
            pb=sum(p_.values())/len(p_)
            HT[name]+=2*pb*(1-pb)
            HS[name]+=sum(2*v*(1-v) for v in p_.values())/len(p_)
print('%-13s %8s %10s %10s %10s %12s'%('locus set','loci','mean FST','median','max','among-stand'))
print('-'*68)
res={}
for name in ('all','single-copy','multi-copy'):
    vals=[]
    for k in num:
        if k[0]!=name: continue
        if den[k]>0: vals.append(max(num[k]/den[k],0.0))
    if not vals: continue
    vals.sort(); res[name]=dict()
    for k in num:
        if k[0]==name and den[k]>0: res[name][(k[1],k[2])]=max(num[k]/den[k],0.0)
    among = 1-HS[name]/HT[name] if HT[name]>0 else float('nan')
    print('%-13s %8d %10.4f %10.4f %10.4f %11.1f%%'%(
        name,NL[name],sum(vals)/len(vals),vals[len(vals)//2],vals[-1],100*among))
print()
if 'single-copy' in res and 'all' in res:
    keys=sorted(set(res['all'])&set(res['single-copy']))
    a=[res['all'][k] for k in keys]; b=[res['single-copy'][k] for k in keys]
    n=len(a); ma=sum(a)/n; mb=sum(b)/n
    sa=math.sqrt(sum((x-ma)**2 for x in a)); sb=math.sqrt(sum((x-mb)**2 for x in b))
    r=sum((x-ma)*(y-mb) for x,y in zip(a,b))/(sa*sb)
    print('ALL vs SINGLE-COPY across %d stand pairs:'%n)
    print('  correlation r = %.4f'%r)
    print('  mean FST      %.4f vs %.4f   (ratio %.3f)'%(ma,mb,mb/ma))
    print('  max abs diff  %.4f'%max(abs(x-y) for x,y in zip(a,b)))
