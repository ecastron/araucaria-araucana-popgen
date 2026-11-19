#!/usr/bin/env python
"""f3 decomposed by whether the triple can actually detect CROSS-GAP admixture.

f3(C; A, B) = mean over SNPs of (p_C - p_A)(p_C - p_B), with Patterson's bias correction
-c(1-c)/(n_C-1). Negative => the target C is a mixture of lineages related to A and B.

Only triples whose target could have received migrants ACROSS the gap test H3'. A triple with all
three stands in the same range tests within-range admixture, which is a different question.
Standard errors by block jackknife over contiguous SNP blocks.
"""
import sys, numpy as np, itertools as it
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
from load_observed import load
from sumstats import freqs_by_deme
G,DEME,COORDS,COASTAL,stands = load()
nd=12; P,N = freqs_by_deme(G,DEME,nd)
NBLK=50
L=P.shape[1]; edges=np.linspace(0,L,NBLK+1).astype(int)

def f3_blocks(c,a,b):
    m=np.isfinite(P[a])&np.isfinite(P[b])&np.isfinite(P[c])&(N[c]>2)
    if m.sum()<200: return None
    num=(P[c]-P[a])*(P[c]-P[b]) - P[c]*(1-P[c])/np.maximum(N[c]-1,1)
    vals=[]
    for k in range(NBLK):
        sl=slice(edges[k],edges[k+1]); mm=m[sl]
        if mm.sum()<5: continue
        vals.append(num[sl][mm].mean())
    vals=np.array(vals)
    if len(vals)<10: return None
    est=num[m].mean()
    # delete-one jackknife
    n=len(vals); pseudo=np.array([ (n*est-(n-1)*np.mean(np.delete(vals,i))) for i in range(n)])
    se=np.sqrt(np.var(pseudo,ddof=1)/n)
    return est, est/se if se>0 else np.nan

CO=[i for i in range(nd) if COASTAL[i]]; AN=[i for i in range(nd) if not COASTAL[i]]
rows=[]
for c in range(nd):
    for a,b in it.combinations([i for i in range(nd) if i!=c],2):
        r=f3_blocks(c,a,b)
        if r is None: continue
        # informative for cross-gap contact: exactly one source on the far side of the target
        cross = (COASTAL[c]!=COASTAL[a]) ^ (COASTAL[c]!=COASTAL[b])
        rows.append((stands[c],stands[a],stands[b],r[0],r[1],cross))
import collections
tot=len(rows)
cross=[r for r in rows if r[5]]; within=[r for r in rows if not r[5]]
print(f"total triples computed : {tot}")
for lab,grp in (("CROSS-GAP informative (target could have received migrants across the gap)",cross),
                ("within-range only (irrelevant to H3')",within)):
    if not grp: continue
    est=np.array([g[3] for g in grp]); z=np.array([g[4] for g in grp])
    neg=(est<0).sum(); sig=((est<0)&(z<-3)).sum()
    print(f"\n{lab}")
    print(f"  n = {len(grp)}")
    print(f"  negative f3            : {neg}   ({100*neg/len(grp):.1f}%)")
    print(f"  significantly negative (Z < -3): {sig}")
    print(f"  most negative f3       : {est.min():+.6f}  (Z = {z[np.argmin(est)]:+.2f})")
    print(f"  minimum Z              : {np.nanmin(z):+.2f}")
    print(f"  median f3              : {np.median(est):+.6f}   median Z = {np.nanmedian(z):+.2f}")
print("\nthe five lowest-Z CROSS-GAP triples (where contact would show first):")
cs=sorted(cross,key=lambda r:r[4])[:5]
for c_,a_,b_,e_,z_,_ in cs:
    print(f"  f3({c_}; {a_}, {b_}) = {e_:+.6f}   Z = {z_:+.2f}")
