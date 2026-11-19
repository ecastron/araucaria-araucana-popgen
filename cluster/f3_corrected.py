#!/usr/bin/env python
"""H2 test: BIAS-CORRECTED f3 (Patterson et al. 2012, eq. for f3-hat).

  f3(C; A,B) = mean_loci [ (c-a)(c-b) - c(1-c)/(n_C - 1) ]

The subtracted term removes the sampling variance of the TARGET population C; without it the
estimator is biased upward and cannot detect admixture. n_C = called chromosomes in C at that
locus. Frequencies only -- no heterozygosity, no SFS. Block jackknife (200 blocks) for SE.
Decision rule fixed in advance: Z <= -3 is an admixture signal.
"""
import numpy as np, itertools, csv, re
V="rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
site={r[0]:r[1] for r in csv.reader(open("metadata/popmap_177.tsv"),delimiter="\t") if len(r)>=2}
keep={l.split()[0] for l in open("results/Q1/feems/keep.txt")}
samples=[]; A_alt=[]; A_n=[]
for line in open(V):
    if line.startswith("##"): continue
    f=line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples=[re.sub(r"\.[12]$","",s) for s in f[9:]]
        stands=sorted({site[s] for s in samples if s in site and s in keep})
        idx={s:[k for k,x in enumerate(samples) if site.get(x)==s and x in keep] for s in stands}
        continue
    if len(f[3])!=1 or len(f[4])!=1 or "," in f[4]: continue
    alt=[];nn=[]
    for s in stands:
        a=n=0
        for k in idx[s]:
            t=f[9+k].split(":",1)[0].replace("|","/")
            if t in ("./.","."): continue
            x,_,y=t.partition("/"); a+=(x!="0")+(y!="0"); n+=2
        alt.append(a); nn.append(n)
    A_alt.append(alt); A_n.append(nn)
ALT=np.array(A_alt,float).T; NN=np.array(A_n,float).T      # stands x loci
good=(NN>=6).all(axis=0)                                    # >=3 individuals per stand
ALT=ALT[:,good]; NN=NN[:,good]
F=ALT/NN
n,L=F.shape
print(f"{n} stands x {L} loci (>=3 genotyped individuals in every stand), {len(keep)} libraries\n")
CORR=F*(1-F)/(NN-1)                                         # per-locus bias term for each stand
NB=200; edges=np.linspace(0,L,NB+1).astype(int)
def f3_all(mask):
    X=F[:,mask]; Cc=CORR[:,mask]; out={}
    for ci in range(n):
        for ai,bi in itertools.combinations([i for i in range(n) if i!=ci],2):
            out[(ci,ai,bi)]=np.mean((X[ci]-X[ai])*(X[ci]-X[bi]) - Cc[ci])
    return out
full=f3_all(np.ones(L,bool))
jack=[f3_all(np.r_[np.ones(edges[b],bool),np.zeros(edges[b+1]-edges[b],bool),np.ones(L-edges[b+1],bool)])
      for b in range(NB)]
res=[]
for k,v in full.items():
    js=np.array([j[k] for j in jack])
    se=np.sqrt((NB-1)/NB*np.sum((js-js.mean())**2))
    res.append((v/se if se>0 else 0.0, v, se, k))
res.sort()
print(f"{'target C':>9} {'source A':>9} {'source B':>9} {'f3':>12} {'SE':>10} {'Z':>8}")
for z,v,se,(c,a,b) in res[:15]:
    print(f"{stands[c]:>9} {stands[a]:>9} {stands[b]:>9} {v:>12.6f} {se:>10.6f} {z:>8.2f}"
          + ("  <- ADMIXED" if z<=-3 else ""))
sig=[r for r in res if r[0]<=-3]
print(f"\ntriples tested: {len(res)} | Z <= -3: {len(sig)} | most negative Z = {res[0][0]:.2f}")
if sig:
    tc={}
    for z,v,se,(c,a,b) in sig: tc[stands[c]]=tc.get(stands[c],0)+1
    print("  stands appearing as ADMIXED target:")
    for s,cnt in sorted(tc.items(), key=lambda kv:-kv[1]): print(f"    {s}: {cnt} triples")
else:
    print("  no admixture signal at Z<=-3")
