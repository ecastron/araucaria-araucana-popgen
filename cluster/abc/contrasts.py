#!/usr/bin/env python
"""Two-stage ABC-RF decomposition.

The three-way comparison identifies H1 well (94.8% OOB) but confuses H2 with H3, so a single
'selected model' misrepresents the evidence. Split it into the two contrasts the data can and
cannot answer:
  stage 1  one refugium (H1) vs two refugia (H2+H3)
  stage 2  persistent isolation (H2) vs secondary contact (H3)
Posterior probabilities follow Pudlo et al. 2016: regress the OOB misclassification indicator on
the statistics and evaluate at the observed point. Vote shares are not posterior probabilities.
"""
import sys, glob, pickle, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
P="/home/ecastron/araucaria"
MODELS=["H1","H2","H3"]; Xs,ys=[],[]; STAT=None
for f in sorted(glob.glob(f"{P}/results/Q1/abc/sims_*.npz")):
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    k=d["y"]!=0
    if k.sum(): Xs.append(d["X"][k]); ys.append(d["y"][k])
    STAT=list(d["stat_order"])
for f in sorted(glob.glob(f"{P}/results/Q1/abc/h1_*.npz")):
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    Xs.append(d["X"]); ys.append(np.zeros(d["X"].shape[0],int)); STAT=list(d["stat_order"])
X=np.vstack(Xs); y=np.concatenate(ys)
obs=pickle.load(open(f"{P}/results/Q1/abc/observed.pkl","rb"))["obs"]
rng=np.random.default_rng(0)
nmin=min((y==i).sum() for i in range(3))
sel=np.concatenate([rng.choice(np.where(y==i)[0],nmin,replace=False) for i in range(3)]); rng.shuffle(sel)
X,y=X[sel],y[sel]

def contrast(Xa,ya,labels,title):
    print(f"\n{'='*74}\n{title}\n{'='*74}")
    print(f"  training: " + ", ".join(f"{l}={int((ya==i).sum())}" for i,l in enumerate(labels)))
    rf=RandomForestClassifier(n_estimators=1000,oob_score=True,n_jobs=-1,random_state=1).fit(Xa,ya)
    po=rf.oob_decision_function_.argmax(1)
    acc=(po==ya).mean()
    print(f"  OOB accuracy: {acc:.3f}  (chance = {1/len(labels):.3f})")
    print(f"  {'':>10}" + "".join(f"{l:>10}" for l in labels))
    for i,l in enumerate(labels):
        row=[(po[ya==i]==j).sum() for j in range(len(labels))]
        print(f"  {l:>10}" + "".join(f"{r:>10}" for r in row) + f"   ({100*row[i]/max(sum(row),1):.1f}% correct)")
    v=rf.predict_proba(obs.reshape(1,-1))[0]
    s=int(v.argmax())
    print("  votes: " + "  ".join(f"{l}={v[i]:.3f}" for i,l in enumerate(labels)))
    rg=RandomForestRegressor(n_estimators=1000,n_jobs=-1,random_state=2).fit(Xa,(po!=ya).astype(float))
    post=1.0-float(rg.predict(obs.reshape(1,-1))[0])
    print(f"  SELECTED: {labels[s]}   posterior probability = {post:.3f}  (Pudlo estimate)")
    imp=np.argsort(rf.feature_importances_)[::-1][:4]
    print("  top statistics: " + ", ".join(f"{STAT[j]} ({rf.feature_importances_[j]:.3f})" for j in imp))
    return labels[s], post, acc

r1=contrast(X,(y!=0).astype(int),["H1 (one refugium)","H2+H3 (two refugia)"],
            "STAGE 1  single refugium vs two refugia")
m=y!=0
r2=contrast(X[m],(y[m]==2).astype(int),["H2 (isolation)","H3 (secondary contact)"],
            "STAGE 2  persistent isolation vs secondary contact")
print(f"\n{'='*74}\nSUMMARY\n{'='*74}")
print(f"  stage 1: {r1[0]}  posterior {r1[1]:.3f}  (classifier OOB {r1[2]:.3f})")
print(f"  stage 2: {r2[0]}  posterior {r2[1]:.3f}  (classifier OOB {r2[2]:.3f})")
pickle.dump(dict(stage1=r1,stage2=r2,stat=STAT,obs=obs),open(f"{P}/results/Q1/abc/abcrf_contrasts.pkl","wb"))
print(f"\nsaved -> results/Q1/abc/abcrf_contrasts.pkl")
