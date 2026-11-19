#!/usr/bin/env python
"""Final ABC analysis: four models, twelve statistics (nine scale-free + three cross-gap)."""
import sys, glob, pickle, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
P="/home/ecastron/araucaria"
LAB=["H1 (Andean source)","H1b (Coastal source)","H2 (isolation)","H3 (contact)"]
Xs,ys,Ps=[],[],[]; STAT=None; PN=None
for f in sorted(glob.glob(f"{P}/results/Q1/abc/s5_*.npz")):
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
print(f"reference table: {X.shape[0]} sims x {X.shape[1]} statistics; balanced to {n}/model\n")

print("=== PRIOR-PREDICTIVE CHECK ===")
print(f"{'statistic':<28}{'observed':>10}{'2.5%':>10}{'50%':>10}{'97.5%':>10}{'pctile':>9}  flag")
bad=[]
for j,s in enumerate(STAT):
    q=np.nanpercentile(Xb[:,j],[2.5,50,97.5]); pc=100*np.mean(Xb[:,j]<=obs[j])
    fl="" if 2.5<=pc<=97.5 else "  <-- OUTSIDE"
    if fl: bad.append(s)
    print(f"{s:<28}{obs[j]:>10.4f}{q[0]:>10.4f}{q[1]:>10.4f}{q[2]:>10.4f}{pc:>8.1f}%{fl}")
print(f"\n  {len(bad)}/{len(STAT)} outside the 95% range" + (f": {bad}" if bad else ""))
USE=[j for j,s in enumerate(STAT) if s not in bad]
if bad:
    print(f"  -> statistics the model cannot reproduce are DROPPED from model choice: {bad}")
    print(f"  -> {len(USE)} statistics retained")

def contrast(Xa,ya,labels,title,obsv):
    print(f"\n{'='*76}\n{title}\n{'='*76}")
    rf=RandomForestClassifier(n_estimators=1000,oob_score=True,n_jobs=-1,random_state=1).fit(Xa,ya)
    po=rf.oob_decision_function_.argmax(1); acc=(po==ya).mean()
    print(f"  OOB accuracy {acc:.3f} (chance {1/len(labels):.3f})")
    print(f"  {'':>26}"+"".join(f"{l[:11]:>13}" for l in labels))
    for i,l in enumerate(labels):
        row=[(po[ya==i]==j).sum() for j in range(len(labels))]
        print(f"  {l:>26}"+"".join(f"{r:>13}" for r in row)+f"  ({100*row[i]/max(sum(row),1):.1f}%)")
    v=rf.predict_proba(obsv.reshape(1,-1))[0]; s=int(v.argmax())
    print("  votes: "+"  ".join(f"{l}={v[i]:.3f}" for i,l in enumerate(labels)))
    rg=RandomForestRegressor(n_estimators=1000,n_jobs=-1,random_state=2).fit(Xa,(po!=ya).astype(float))
    post=1.0-float(rg.predict(obsv.reshape(1,-1))[0])
    print(f"  SELECTED: {labels[s]}   posterior {post:.3f}")
    imp=np.argsort(rf.feature_importances_)[::-1][:5]
    print("  top statistics: "+", ".join(f"{[STAT[k] for k in USE][j]} ({rf.feature_importances_[j]:.3f})" for j in imp))
    return post,acc
Xu=Xb[:,USE]; obsu=obs[USE]
contrast(Xu,yb,LAB,"FOUR-WAY",obsu)
contrast(Xu,(~np.isin(yb,[0,1])).astype(int),["one refugium","two refugia"],
         "STAGE 1  single refugium (either direction) vs two refugia",obsu)
m=np.isin(yb,[2,3])
contrast(Xu[m],(yb[m]==3).astype(int),["H2 isolation","H3 contact"],
         "STAGE 2  isolation vs secondary contact",obsu)

print(f"\n{'='*76}\nDETECTION LIMIT AND MIGRATION RATE (12 statistics)\n{'='*76}")
imc=PN.index("mcross"); iNe=PN.index("Ne"); itc=PN.index("Tcontact")
h3=np.where(yb==3)[0]; h2=np.where(yb==2)[0]
Nm=Pb[h3,iNe]*Pb[h3,imc]
Xc=np.vstack([Xu[h2],Xu[h3]]); yc=np.r_[np.zeros(len(h2),int),np.ones(len(h3),int)]
rf=RandomForestClassifier(n_estimators=1000,oob_score=True,n_jobs=-1,random_state=1).fit(Xc,yc)
p3=rf.oob_decision_function_.argmax(1)[len(h2):]
print(f"{'corridor Nm':>26}{'n':>7}{'% called H3':>13}")
for lo,hi,lab in ((0,.2,"< 0.2"),(.2,1,"0.2-1"),(1,5,"1-5"),(5,1e9,"> 5")):
    k=(Nm>=lo)&(Nm<hi)
    if k.sum()>=10: print(f"{lab:>26}{k.sum():>7}{100*p3[k].mean():>12.1f}%")
tgt=np.log10(np.maximum(Nm,1e-9))
rg=RandomForestRegressor(n_estimators=1000,n_jobs=-1,random_state=3,oob_score=True).fit(Xu[h3],tgt)
lv=rg.apply(Xu[h3]); ol=rg.apply(obsu.reshape(1,-1))[0]
w=np.zeros(len(h3))
for t in range(lv.shape[1]):
    sm=lv[:,t]==ol[t]
    if sm.sum(): w[sm]+=1.0/sm.sum()
w/=w.sum(); o=np.argsort(tgt); cw=np.cumsum(w[o])
q=lambda pp: 10**tgt[o][np.searchsorted(cw,pp)]
print(f"\n  regression OOB R2 on log10(Nm): {rg.oob_score_:.3f}")
print(f"  corridor Nm at the observed point: {10**float(rg.predict(obsu.reshape(1,-1))[0]):.4f} migrants/gen")
print(f"  95% CI [{q(0.025):.4f}, {q(0.975):.4f}]   median {q(0.5):.4f}")
