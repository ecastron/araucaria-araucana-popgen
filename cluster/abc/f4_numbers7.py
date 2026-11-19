#!/usr/bin/env python
"""Every number F4 prints, recomputed with the SAME procedure as final_abc.py.

The posterior is the ABC-RF posterior of Pudlo et al. (2016) -- a regression
forest fitted to the out-of-bag misclassification indicator and evaluated at the
observed point -- NOT the raw vote proportion, which is a different and biased
quantity. Vote shares are reported separately because they are what panel A
displays per model.
"""
import glob, pickle, json, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
P="/home/ecastron/araucaria"
LAB=["H1' Andean source","H1b' Coastal source","H2' isolation","H3' contact"]
Xs,ys,Ps=[],[],[]; STAT=None; PN=None
for f in sorted(glob.glob(f"{P}/results/Q1/abc_floor10/s7_*.npz")):
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    Xs.append(d["X"]); ys.append(d["y"]); Ps.append(d["P"])
    STAT=list(d["stat_order"]); PN=list(d["param_names"])
X=np.vstack(Xs); y=np.concatenate(ys); PP=np.vstack(Ps)
obs=pickle.load(open(f"{P}/results/Q1/abc/observed2.pkl","rb"))["obs"]
rng=np.random.default_rng(0)
n=min((y==i).sum() for i in range(4))
sel=np.concatenate([rng.choice(np.where(y==i)[0],n,replace=False) for i in range(4)]); rng.shuffle(sel)
Xb,yb,Pb=X[sel],y[sel],PP[sel]
bad=[s for j,s in enumerate(STAT) if not 2.5<=100*np.mean(Xb[:,j]<=obs[j])<=97.5]
USE=[j for j,s in enumerate(STAT) if s not in bad]
out={"n_per_model":int(n),"n_sims":int(Xb.shape[0]),"dropped":bad,"n_stats":len(USE)}

def contrast(Xa,ya,labels,obsv):
    rf=RandomForestClassifier(n_estimators=1000,oob_score=True,n_jobs=-1,random_state=1).fit(Xa,ya)
    po=rf.oob_decision_function_.argmax(1); acc=float((po==ya).mean())
    v=rf.predict_proba(obsv.reshape(1,-1))[0]
    rg=RandomForestRegressor(n_estimators=1000,n_jobs=-1,random_state=2).fit(Xa,(po!=ya).astype(float))
    post=1.0-float(rg.predict(obsv.reshape(1,-1))[0])
    cm=np.zeros((len(labels),len(labels)),int)
    for i in range(len(labels)):
        for j in range(len(labels)): cm[i,j]=int(((ya==i)&(po==j)).sum())
    return dict(votes=[float(x) for x in v],posterior=post,oob=acc,
                selected=labels[int(v.argmax())],cm=cm.tolist())

Xu=Xb[:,USE]; obsu=obs[USE]
out["four"]=contrast(Xu,yb,LAB,obsu)
out["stage1"]=contrast(Xu,(~np.isin(yb,[0,1])).astype(int),["one refugium","two refugia"],obsu)
m=np.isin(yb,[2,3])
out["stage2"]=contrast(Xu[m],(yb[m]==3).astype(int),["isolation","contact"],obsu)

# Divergence is estimated ONLY under the selected model class. T1 is "T_anc" in the
# one-refugium models and "T_div" in the two-refugium models -- different quantities,
# so pooling all four models estimates nothing.
iT=PN.index("T1"); iNe=PN.index("Ne")
sel2=np.isin(yb,[2,3]); X2=Xu[sel2]; P2=Pb[sel2]
tgt=np.log10(P2[:,iT]/(2*P2[:,iNe]))
rg=RandomForestRegressor(n_estimators=800,oob_score=True,n_jobs=-1,random_state=3).fit(X2,tgt)
out["tau_med"]=float(10**rg.predict(obsu.reshape(1,-1))[0]); out["tau_r2"]=float(rg.oob_score_)
lvt=np.array([t.predict(X2) for t in rg.estimators_]).T
olt=np.array([t.predict(obsu.reshape(1,-1))[0] for t in rg.estimators_])
wt=np.zeros(len(tgt))
for t in range(lvt.shape[1]):
    sm=lvt[:,t]==olt[t]
    if sm.sum(): wt[sm]+=1.0/sm.sum()
wt/=wt.sum(); ot=np.argsort(tgt); cwt=np.cumsum(wt[ot])
qt=lambda pp: float(10**tgt[ot][np.searchsorted(cwt,pp)])
out["tau_lo"],out["tau_hi"]=qt(0.025),qt(0.975)
tgtT=np.log10(P2[:,iT])
rgT=RandomForestRegressor(n_estimators=800,oob_score=True,n_jobs=-1,random_state=3).fit(X2,tgtT)
out["T_med"]=float(10**rgT.predict(obsu.reshape(1,-1))[0]); out["T_r2"]=float(rgT.oob_score_)
lv=np.array([t.predict(X2) for t in rgT.estimators_]).T
ol=np.array([t.predict(obsu.reshape(1,-1))[0] for t in rgT.estimators_])
w=np.zeros(len(tgtT))
for t in range(lv.shape[1]):
    sm=lv[:,t]==ol[t]
    if sm.sum(): w[sm]+=1.0/sm.sum()
w/=w.sum(); o=np.argsort(tgtT); cw=np.cumsum(w[o])
q=lambda pp: float(10**tgtT[o][np.searchsorted(cw,pp)])
out["T_lo"],out["T_hi"]=q(0.025),q(0.975)
out["prior_T"]=[float(10**1.0),float(10**4.8)]

# posterior-predictive percentile of the observed value under the two-refugia models
pp={}
for j,s in enumerate(STAT):
    v=Xb[np.isin(yb,[2,3]),j]; v=v[np.isfinite(v)]
    pp[s]=float(100*np.mean(v<=obs[j])) if len(v) else float("nan")
out["pp_pct"]=pp; out["stat"]=STAT
json.dump(out,open(f"{P}/results/Q1/abc_floor10/f4_numbers_floor10.json","w"),indent=1)
print(json.dumps({k:v for k,v in out.items() if k not in ("pp_pct","four","stage1","stage2")},indent=1))
print("\nfour-way votes:", ["%s %.1f%%"%(LAB[i],100*v) for i,v in enumerate(out["four"]["votes"])])
print("stage1 posterior %.3f oob %.3f"%(out["stage1"]["posterior"],out["stage1"]["oob"]))
print("stage2 posterior %.3f oob %.3f (selected %s)"%(out["stage2"]["posterior"],out["stage2"]["oob"],out["stage2"]["selected"]))
print("\nposterior-predictive percentiles:")
for s in STAT: print("  %-26s %6.1f"%(s,pp[s]))

# ---- added 17 Nov 2026: posterior mode, HDI, and mass below the OLD floor (1,000 generations) ----
def mode_hdi(target, wts, bw=0.06):
    grid=np.linspace(target.min(),target.max(),2000)
    dens=(wts[None,:]*np.exp(-0.5*((grid[:,None]-target[None,:])/bw)**2)).sum(1)
    od=np.argsort(dens)[::-1]; cs=np.cumsum(dens[od])/dens.sum()
    h50=grid[od[:np.searchsorted(cs,0.5)+1]]; h95=grid[od[:np.searchsorted(cs,0.95)+1]]
    return 10**grid[dens.argmax()], (10**h50.min(),10**h50.max()), (10**h95.min(),10**h95.max())
mT,h50T,h95T=mode_hdi(tgtT,w); mt,h50t,h95t=mode_hdi(tgt,wt)
extra={"T_mode":mT,"T_hdi50":h50T,"T_hdi95":h95T,"T_mass_below_1000":float(w[P2[:,iT]<1000].sum()),
       "tau_mode":mt,"tau_hdi50":h50t,"tau_hdi95":h95t}
out.update(extra); json.dump(out,open(f"{P}/results/Q1/abc_floor10/f4_numbers_floor10.json","w"),indent=1)
print(json.dumps(extra,indent=1))
