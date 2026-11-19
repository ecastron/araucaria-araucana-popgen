#!/usr/bin/env python
"""MC1: block-jackknife CIs on every pairwise FST. The non-overlap claim rests on a gap of
0.0065 between the largest within-range value and the smallest between-range value; without
interval estimates it is not testable."""
import sys, numpy as np, itertools as it
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
from load_observed import load
from sumstats import freqs_by_deme, hudson_fst
from chain import _km
G,DEME,COORDS,COASTAL,stands=load()
P,N=freqs_by_deme(G,DEME,12)
L=P.shape[1]; NB=50; edges=np.linspace(0,L,NB+1).astype(int)
rows=[]
for a,b in it.combinations(range(12),2):
    m=(N[a]>=6)&(N[b]>=6)
    if m.sum()<50: continue
    full=hudson_fst(P[a][m],P[b][m],N[a][m],N[b][m])
    vals=[]
    for k in range(NB):                       # delete-one-block jackknife
        sl=np.ones(L,bool); sl[edges[k]:edges[k+1]]=False
        mm=m&sl
        if mm.sum()<50: continue
        v=hudson_fst(P[a][mm],P[b][mm],N[a][mm],N[b][mm])
        if np.isfinite(v): vals.append(v)
    vals=np.array(vals); n=len(vals)
    se=np.sqrt((n-1)/n*np.sum((vals-vals.mean())**2)) if n>5 else np.nan
    rows.append((stands[a],stands[b],_km(COORDS[a],COORDS[b]),full,se,
                 int(COASTAL[a]!=COASTAL[b]),int(m.sum())))
w=[r for r in rows if not r[5]]; b=[r for r in rows if r[5]]
wmax=max(w,key=lambda r:r[3]); bmin=min(b,key=lambda r:r[3])
print(f"pairs: {len(w)} within, {len(b)} between   ({NB}-block delete-one jackknife)\n")
print(f"largest WITHIN-range  : {wmax[0]}-{wmax[1]}  FST = {wmax[3]:.4f} +/- {wmax[4]:.4f}  "
      f"95% CI [{wmax[3]-1.96*wmax[4]:.4f}, {wmax[3]+1.96*wmax[4]:.4f}]")
print(f"smallest BETWEEN-range: {bmin[0]}-{bmin[1]}  FST = {bmin[3]:.4f} +/- {bmin[4]:.4f}  "
      f"95% CI [{bmin[3]-1.96*bmin[4]:.4f}, {bmin[3]+1.96*bmin[4]:.4f}]")
gap=bmin[3]-wmax[3]; se=np.sqrt(wmax[4]**2+bmin[4]**2)
print(f"\ngap = {gap:+.4f}   SE of the difference = {se:.4f}   z = {gap/se:.2f}   "
      f"{'SEPARATED' if gap-1.96*se>0 else 'CIs OVERLAP'}")
ov=sum(1 for x in w for yb in b if x[3]+1.96*x[4] >= yb[3]-1.96*yb[4])
print(f"within/between pair comparisons whose 95% CIs overlap: {ov} of {len(w)*len(b)}")
print(f"\nSE range across all pairs: {min(r[4] for r in rows):.5f} - {max(r[4] for r in rows):.5f}")
print(f"loci per pair: {min(r[6] for r in rows):,} - {max(r[6] for r in rows):,}")
with open("/home/ecastron/araucaria/results/figdata/host_fst_pairs_ci.tsv","w") as f:
    f.write("stand_a\tstand_b\tkm\tfst\tse\tlo95\thi95\tbetween\tnloci\n")
    for r in sorted(rows,key=lambda x:x[3]):
        f.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.6f}\t{r[4]:.6f}\t"
                f"{r[3]-1.96*r[4]:.6f}\t{r[3]+1.96*r[4]:.6f}\t{r[5]}\t{r[6]}\n")
print("\nwritten -> results/figdata/host_fst_pairs_ci.tsv")
