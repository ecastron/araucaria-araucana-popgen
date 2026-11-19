#!/usr/bin/env python
"""Three controls a referee will demand.

MC3  Calibrate the migration-rate regressor against ZERO-migration data. The regressor is trained
     on H3' only, so it structurally cannot return zero. Applying it to H2' simulations (true
     Nm = 0) establishes its floor: if H2' data yield estimates around 0.154, the observed value
     is the estimator's floor rather than a measurement.
MC2  Posterior-predictive check. The prior-predictive check asks whether the model CLASS can reach
     the data; this asks whether the SELECTED model at its estimated parameters reproduces it.
MC1  Block-jackknife confidence intervals on every pairwise FST, to test whether the
     between/within non-overlap survives uncertainty (the gap is only 0.0065).
"""
import sys, numpy as np, pickle, itertools as it
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
from sklearn.ensemble import RandomForestRegressor
P="/home/ecastron/araucaria"
d=np.load(f"{P}/results/figdata/abc_for_fig.npz",allow_pickle=True)
X,y,PP,obs=d["X"],d["y"],d["P"],d["obs"]
STAT=[str(s) for s in d["stat"]]; PN=[str(s) for s in d["param"]]
imc,iNe=PN.index("mcross"),PN.index("Ne")
h2=np.where(y==2)[0]; h3=np.where(y==3)[0]
Nm=PP[h3,iNe]*PP[h3,imc]
tgt=np.log10(np.maximum(Nm,1e-9))

print("="*78); print("MC3  CALIBRATION OF THE MIGRATION-RATE REGRESSOR"); print("="*78)
rg=RandomForestRegressor(n_estimators=1000,n_jobs=-1,random_state=3,oob_score=True).fit(X[h3],tgt)
print(f"  training OOB R2 (on H3' where true Nm is known): {rg.oob_score_:.3f}")
est_obs=10**float(rg.predict(obs.reshape(1,-1))[0])
print(f"  estimate at the OBSERVED data                  : {est_obs:.4f} migrants/gen")
# apply to H2' simulations, whose true migration is exactly zero
pred_h2=10**rg.predict(X[h2])
q=np.percentile(pred_h2,[2.5,25,50,75,97.5])
print(f"\n  applied to {len(h2)} H2' simulations with TRUE Nm = 0:")
print(f"     median {q[2]:.4f}   IQR [{q[1]:.4f}, {q[3]:.4f}]   95% [{q[0]:.4f}, {q[4]:.4f}]")
frac=100*np.mean(pred_h2>=est_obs)
print(f"     fraction of zero-migration sims estimated at or ABOVE the observed {est_obs:.4f}: {frac:.1f}%")
print(f"\n  VERDICT: " + ("the observed estimate is INDISTINGUISHABLE from the estimator's zero-migration floor"
      if frac>10 else "the observed estimate is ABOVE what zero-migration data produce"))
# and the converse: what do LOW but non-zero H3' sims give?
lo=h3[(Nm>0.05)&(Nm<0.5)]
if len(lo)>30:
    pl=10**rg.predict(X[lo])
    print(f"  for reference, H3' sims with true Nm in 0.05-0.5 (n={len(lo)}): median estimate {np.median(pl):.4f}")

print("\n"+"="*78); print("MC2  POSTERIOR-PREDICTIVE CHECK"); print("="*78)
# retain simulations closest to the observed point under the SELECTED class (two refugia),
# then ask whether they reproduce each statistic
two=np.where(np.isin(y,[2,3]))[0]
sd=np.nanstd(X[two],axis=0); sd[sd==0]=1
Z=(X[two]-obs)/sd
dist=np.sqrt(np.nansum(Z**2,axis=1))
keep=two[np.argsort(dist)[:int(0.02*len(two))]]        # closest 2% = posterior sample
print(f"  posterior sample: {len(keep)} simulations closest to the observed vector\n")
print(f"  {'statistic':<28}{'observed':>10}{'post 2.5%':>11}{'post 50%':>10}{'post 97.5%':>12}{'pctile':>8}  flag")
bad=[]
for j,s in enumerate(STAT):
    v=X[keep,j]; v=v[np.isfinite(v)]
    if len(v)<10: continue
    qq=np.percentile(v,[2.5,50,97.5]); pc=100*np.mean(v<=obs[j])
    fl="" if 2.5<=pc<=97.5 else "  <-- OUTSIDE"
    if fl: bad.append(s)
    print(f"  {s:<28}{obs[j]:>10.4f}{qq[0]:>11.4f}{qq[1]:>10.4f}{qq[2]:>12.4f}{pc:>7.1f}%{fl}")
print(f"\n  {len(bad)}/{len(STAT)} statistics outside the 95% POSTERIOR-predictive range"
      + (f": {bad}" if bad else " -- the selected model reproduces the data"))
