#!/usr/bin/env python
"""Four-model ABC-RF including BOTH colonisation directions.

H1'  single refugium, ANDEAN source -> Coastal derived
H1b' single refugium, COASTAL source -> Andean derived   (the 'or vice versa' case)
H2'  two refugia, persistent isolation
H3'  two refugia, secondary contact

Reported as contrasts, because a single 'selected model' hides which comparisons the data can
actually resolve. Posteriors follow Pudlo et al. 2016 (RF regression on the OOB error indicator).
"""
import sys, glob, pickle, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
P="/home/ecastron/araucaria"
LAB=["H1 (Andean source)","H1b (Coastal source)","H2 (isolation)","H3 (contact)"]
Xs,ys=[],[]; STAT=None
for f in sorted(glob.glob(f"{P}/results/Q1/abc/sims_*.npz")):     # H2/H3 only (H1 there was buggy)
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    k=d["y"]!=0
    if k.sum(): Xs.append(d["X"][k]); ys.append(d["y"][k]+1)      # 1,2 -> 2,3
    STAT=list(d["stat_order"])
for f in sorted(glob.glob(f"{P}/results/Q1/abc/h1_*.npz")):
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    Xs.append(d["X"]); ys.append(np.zeros(d["X"].shape[0],int)); STAT=list(d["stat_order"])
for f in sorted(glob.glob(f"{P}/results/Q1/abc/h1b_*.npz")):
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    Xs.append(d["X"]); ys.append(np.ones(d["X"].shape[0],int)); STAT=list(d["stat_order"])
X=np.vstack(Xs); y=np.concatenate(ys)
obs=pickle.load(open(f"{P}/results/Q1/abc/observed.pkl","rb"))["obs"]
rng=np.random.default_rng(0)
n=min((y==i).sum() for i in range(4))
sel=np.concatenate([rng.choice(np.where(y==i)[0],n,replace=False) for i in range(4)]); rng.shuffle(sel)
X,y=X[sel],y[sel]
print(f"balanced reference table: {n} per model, {X.shape[0]} total\n")

pr=STAT.index("private_ratio")
print("private_ratio by model (the statistic that separates colonisation direction):")
for i,l in enumerate(LAB):
    v=X[y==i,pr]
    print(f"  {l:24s} median {np.median(v):7.3f}   2.5-97.5% [{np.percentile(v,2.5):6.3f}, {np.percentile(v,97.5):8.3f}]")
print(f"  {'OBSERVED':24s}        {obs[pr]:7.3f}")

def contrast(Xa,ya,labels,title):
    print(f"\n{'='*76}\n{title}\n{'='*76}")
    rf=RandomForestClassifier(n_estimators=1000,oob_score=True,n_jobs=-1,random_state=1).fit(Xa,ya)
    po=rf.oob_decision_function_.argmax(1); acc=(po==ya).mean()
    print(f"  OOB accuracy {acc:.3f} (chance {1/len(labels):.3f})")
    print(f"  {'':>24}"+"".join(f"{l[:10]:>12}" for l in labels))
    for i,l in enumerate(labels):
        row=[(po[ya==i]==j).sum() for j in range(len(labels))]
        print(f"  {l:>24}"+"".join(f"{r:>12}" for r in row)+f"  ({100*row[i]/max(sum(row),1):.1f}%)")
    v=rf.predict_proba(obs.reshape(1,-1))[0]; s=int(v.argmax())
    print("  votes: "+"  ".join(f"{l}={v[i]:.3f}" for i,l in enumerate(labels)))
    rg=RandomForestRegressor(n_estimators=1000,n_jobs=-1,random_state=2).fit(Xa,(po!=ya).astype(float))
    post=1.0-float(rg.predict(obs.reshape(1,-1))[0])
    print(f"  SELECTED: {labels[s]}   posterior {post:.3f}")
    return labels[s],post,acc

contrast(X,y,LAB,"FOUR-WAY")
one=np.isin(y,[0,1])
contrast(X,(~one).astype(int),["one refugium (H1 or H1b)","two refugia (H2 or H3)"],
         "STAGE 1  single refugium, EITHER direction, vs two refugia")
contrast(X[one],y[one],["H1 (Andean source)","H1b (Coastal source)"],
         "DIRECTION  which range would be the source, if there were only one?")
