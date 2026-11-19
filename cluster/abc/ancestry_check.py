#!/usr/bin/env python
"""Does a NO-GENE-FLOW model reproduce the K=2 ancestry barplot?

sNMF/ADMIXTURE decompose variation into K components; every individual is forced to be a mixture,
so intermediate coefficients arise from drift, unbalanced sampling and isolation by distance with
no gene flow whatsoever (Lawson et al. 2018). If H2' (two refugia, ZERO cross-gap migration)
reproduces the observed per-stand ancestry proportions, then Fig 2E is not evidence of admixture
and does not conflict with f3.

Decomposition here is NMF with 2 components on the mean-imputed dosage matrix -- the same model
family as sNMF, applied identically to observed and simulated data.
"""
import sys, numpy as np, pickle
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
from sklearn.decomposition import NMF
import sim5 as S
from load_observed import load

def ancestry(Gd, deme, nd=12):
    X = Gd.copy().astype(float)
    cm = np.nanmean(X, axis=0); ix = np.where(np.isnan(X))
    X[ix] = np.take(cm, ix[1])
    X = X[:, np.isfinite(X).all(0) & (X.std(0) > 0)]
    X = X - X.min()
    W = NMF(n_components=2, init="nndsvda", max_iter=400, random_state=1).fit_transform(X)
    q = W[:, 0] / np.maximum(W.sum(1), 1e-12)
    # orient so component 1 is the one enriched in the Coastal demes
    co = np.isin(deme, [i for i in range(nd) if S.COASTAL[i]])
    if q[co].mean() < q[~co].mean(): q = 1 - q
    return np.array([q[deme == i].mean() for i in range(nd)]), q

G,DEME,COORDS,COASTAL,stands = load()
obs_stand, obs_q = ancestry(G, DEME)
print("OBSERVED per-stand ancestry (component enriched in Coastal):")
order = np.argsort(obs_stand)[::-1]
print("  " + "  ".join(f"{stands[i]}={obs_stand[i]:.3f}" for i in order))
an = [i for i in range(12) if not COASTAL[i]]
print(f"  Andean stands: min {obs_stand[an].min():.3f}  max {obs_stand[an].max():.3f}  "
      f"spread {obs_stand[an].max()-obs_stand[an].min():.3f}")

print("\nSIMULATED under H2' (two refugia, ZERO cross-gap migration):")
p = dict(Ne=3000, m=0.003, T1=10000, Nanc=150000, Tsep=10000)
sims = []
for k, seed in enumerate((81001, 81002, 81003)):
    Gd = S.simulate("H2", p, seed)
    if Gd is None: continue
    st, _ = ancestry(Gd, S.SIM_DEME)
    sims.append(st)
    o = np.argsort(st)[::-1]
    print(f"  rep {k+1}: " + "  ".join(f"{stands[i]}={st[i]:.3f}" for i in o))
    print(f"          Andean stands: min {st[an].min():.3f}  max {st[an].max():.3f}  "
          f"spread {st[an].max()-st[an].min():.3f}")
if sims:
    M = np.vstack(sims)
    print("\nmean over replicates, Andean stands only:")
    print(f"  simulated  {M[:,an].mean():.3f}  (range {M[:,an].min():.3f}-{M[:,an].max():.3f})")
    print(f"  observed   {obs_stand[an].mean():.3f}  (range {obs_stand[an].min():.3f}-{obs_stand[an].max():.3f})")
    print("\nCoastal stands:")
    co = [i for i in range(12) if COASTAL[i]]
    for i in co:
        print(f"  {stands[i]}: observed {obs_stand[i]:.3f}   simulated {M[:,i].mean():.3f}")
    print("\n-> if the simulated Andean stands carry a comparable non-zero 'Coastal' component,")
    print("   the observed barplot is reproduced by a model with NO gene flow at all.")
