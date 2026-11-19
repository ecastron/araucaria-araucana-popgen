#!/usr/bin/env python
"""Detection limit for cross-gap gene flow, and posterior for the migration rate itself.

H2' is H3' with mcross -> 0: the models are NESTED, so a model-choice contrast is ill-posed.
Every H3' draw with negligible migration is physically identical to H2' and MUST be misclassified.
Rather than ask 'H2 or H3', ask how much gene flow is detectable, and what the data imply about it.

Gene flow is summarised as migrants per generation through the corridor, Nm = Ne * mcross,
and as the total exchanged since contact, Ne * mcross * Tcontact.
"""
import sys, glob, pickle, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
P="/home/ecastron/araucaria"
X2,X3,PP=[],[],[]
STAT=None
for f in sorted(glob.glob(f"{P}/results/Q1/abc/sims_*.npz")):
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    STAT=list(d["stat_order"]); names=list(d["param_names"])
    y=d["y"]; Xf=d["X"]; Pf=d["P"]
    X2.append(Xf[y==1])
    X3.append(Xf[y==2]); PP.append(Pf[y==2])
X2=np.vstack(X2); X3=np.vstack(X3); PP=np.vstack(PP)
iNe=names.index("Ne"); imc=names.index("mcross"); itc=names.index("Tcontact")
Ne=PP[:,iNe]; mc=PP[:,imc]; tc=PP[:,itc]
Nm = Ne*mc                       # migrants per generation through the corridor
tot= Ne*mc*tc                    # total migrants exchanged since contact
print(f"H2 sims {X2.shape[0]}   H3 sims {X3.shape[0]}")
print(f"H3 prior on corridor Nm : {np.percentile(Nm,[2.5,50,97.5])}")
print(f"H3 prior on total migrants since contact: {np.percentile(tot,[2.5,50,97.5])}\n")

X=np.vstack([X2,X3]); y=np.r_[np.zeros(len(X2),int), np.ones(len(X3),int)]
rf=RandomForestClassifier(n_estimators=1000,oob_score=True,n_jobs=-1,random_state=1).fit(X,y)
po=rf.oob_decision_function_.argmax(1)
p3=po[len(X2):]                  # OOB calls for the H3 simulations only
print("DETECTION LIMIT: how often is a true H3' simulation recognised as H3',")
print("as a function of how much gene flow it actually had?\n")
print(f"{'corridor Nm (migrants/gen)':>34} {'n':>6} {'% called H3':>12}")
edges=[0,1e-3,1e-2,0.05,0.2,1,5,1e9]
lab=["< 0.001","0.001-0.01","0.01-0.05","0.05-0.2","0.2-1","1-5","> 5"]
for i in range(len(edges)-1):
    m=(Nm>=edges[i])&(Nm<edges[i+1])
    if m.sum()<10: continue
    print(f"{lab[i]:>34} {m.sum():>6} {100*p3[m].mean():>11.1f}%")
print()
print(f"{'total migrants since contact':>34} {'n':>6} {'% called H3':>12}")
edges2=[0,1,10,100,1e3,1e4,1e9]; lab2=["< 1","1-10","10-100","100-1e3","1e3-1e4","> 1e4"]
for i in range(len(edges2)-1):
    m=(tot>=edges2[i])&(tot<edges2[i+1])
    if m.sum()<10: continue
    print(f"{lab2[i]:>34} {m.sum():>6} {100*p3[m].mean():>11.1f}%")

# where is the observed point?
obs=pickle.load(open(f"{P}/results/Q1/abc/observed.pkl","rb"))["obs"]
print("\n" + "="*72)
print("POSTERIOR FOR THE CORRIDOR MIGRATION RATE (RF regression, Raynal et al. 2019)")
print("="*72)
# train on H3 sims only: statistics -> log10 Nm, then predict at the observed point
tgt=np.log10(np.maximum(Nm,1e-9))
rg=RandomForestRegressor(n_estimators=1000,n_jobs=-1,random_state=3,oob_score=True).fit(X3,tgt)
print(f"  regression OOB R^2 on log10(Nm): {rg.oob_score_:.3f}   "
      f"({'informative' if rg.oob_score_>0.2 else 'NOT informative -- statistics carry little information about the rate'})")
pred=float(rg.predict(obs.reshape(1,-1))[0])
# quantile posterior from the forest's leaf-weighted empirical distribution
leaves=rg.apply(X3); oleaf=rg.apply(obs.reshape(1,-1))[0]
w=np.zeros(len(X3))
for t in range(leaves.shape[1]):
    same=leaves[:,t]==oleaf[t]
    if same.sum(): w[same]+=1.0/same.sum()
w/=w.sum()
o=np.argsort(tgt); cw=np.cumsum(w[o])
q=lambda p: 10**tgt[o][np.searchsorted(cw,p)]
print(f"  point estimate  corridor Nm = {10**pred:.4f} migrants/generation")
print(f"  95% credible interval       = [{q(0.025):.4f}, {q(0.975):.4f}]")
print(f"  median                      = {q(0.5):.4f}")
