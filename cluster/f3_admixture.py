#!/usr/bin/env python
"""H2 test: f3(C; A,B) admixture statistics from per-stand allele frequencies.

f3(C; A,B) = mean_loci (c - a)(c - b).  Significantly NEGATIVE means C's frequencies sit
BETWEEN A and B more often than drift from a tree allows -> C is admixed.  Positive or zero
is consistent with a tree (H1 or H3).  Uses allele frequencies only: no heterozygosity, no SFS.
SE by block jackknife over loci.  Decision rule fixed in advance: Z <= -3 is a signal.
"""
import numpy as np, itertools, csv, sys
D="/home/ecastron/araucaria/results/Q1"
rows=list(csv.reader(open(f"{D}/construct/freqs_all.csv")))
stands=[r[0] for r in rows[1:]]
F=np.array([[np.nan if x=="NA" else float(x) for x in r[1:]] for r in rows[1:]])
ok=~np.isnan(F).any(axis=0)
F=F[:,ok]
n,L=F.shape
print(f"{n} stands x {L} loci complete in all stands\n")
NB=200
edges=np.linspace(0,L,NB+1).astype(int)
def f3_all(mask):
    X=F[:,mask]
    out={}
    for ci in range(n):
        for ai,bi in itertools.combinations([i for i in range(n) if i!=ci],2):
            out[(ci,ai,bi)]=np.mean((X[ci]-X[ai])*(X[ci]-X[bi]))
    return out
full=f3_all(np.ones(L,bool))
jack=[]
for b in range(NB):
    m=np.ones(L,bool); m[edges[b]:edges[b+1]]=False
    jack.append(f3_all(m))
res=[]
for k,v in full.items():
    js=np.array([j[k] for j in jack])
    se=np.sqrt((NB-1)/NB*np.sum((js-js.mean())**2))
    res.append((v/se if se>0 else 0.0, v, se, k))
res.sort()
print(f"{'target C':>9} {'source A':>9} {'source B':>9} {'f3':>12} {'SE':>10} {'Z':>8}")
sig=[r for r in res if r[0]<=-3]
for z,v,se,(c,a,b) in res[:15]:
    flag="  <- ADMIXED" if z<=-3 else ""
    print(f"{stands[c]:>9} {stands[a]:>9} {stands[b]:>9} {v:>12.6f} {se:>10.6f} {z:>8.2f}{flag}")
print(f"\ntriples tested: {len(res)} | significantly negative (Z <= -3): {len(sig)}")
if sig:
    tc={}
    for z,v,se,(c,a,b) in sig: tc[stands[c]]=tc.get(stands[c],0)+1
    print("  stands appearing as the ADMIXED target:")
    for s,cnt in sorted(tc.items(), key=lambda kv:-kv[1]): print(f"    {s}: {cnt} triples")
else:
    print("  NO stand shows an admixture signal -> H2 (secondary contact) not supported;")
    print("  the data are consistent with a tree/cline (H1 or H3).")
