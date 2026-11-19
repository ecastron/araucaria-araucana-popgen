#!/usr/bin/env python
"""ABC random forest for model choice (Pudlo et al. 2016, Bioinformatics 32:859).

Step 1  RF classifier on (summary statistics -> model), OOB error as the honest error rate.
Step 2  Predict the model at the observed summary vector.
Step 3  Posterior probability of the selected model, obtained by REGRESSING the out-of-bag
        misclassification indicator on the statistics with a second RF and evaluating that
        regression at the observed point. The vote share is NOT a posterior probability.

Before any of that: a prior-predictive check. If the observed vector lies outside the simulated
cloud, the classification is extrapolation and the result is not interpretable -- this is reported
first and loudly, because it is the failure mode that actually bites.
"""
import sys, glob, pickle, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

P = "/home/ecastron/araucaria"
MODELS = ["H1", "H2", "H3"]
STAT = None
Xs, ys = [], []
# H1 in the first production run was corrupted by the inactive-population bug: only 20.5% of draws
# completed, and they survived only when their parameters happened not to trigger it, so that
# subset is biased on the parameters. Drop it and use the cleanly regenerated H1 files instead.
for f in sorted(glob.glob(f"{P}/results/Q1/abc/sims_*.npz")):
    d = np.load(f, allow_pickle=True)
    if d["X"].size == 0: continue
    keep = d["y"] != 0
    if keep.sum(): Xs.append(d["X"][keep]); ys.append(d["y"][keep])
    STAT = list(d["stat_order"])
for f in sorted(glob.glob(f"{P}/results/Q1/abc/h1_*.npz")):
    d = np.load(f, allow_pickle=True)
    if d["X"].size == 0: continue
    Xs.append(d["X"]); ys.append(np.zeros(d["X"].shape[0], dtype=int))
    STAT = list(d["stat_order"])
if not Xs: sys.exit("no simulation files")
X = np.vstack(Xs); y = np.concatenate(ys)
# balance the classes so the model prior is uniform, as ABC model choice assumes
nmin = min((y == i).sum() for i in range(len(MODELS)))
rng0 = np.random.default_rng(0)
sel = np.concatenate([rng0.choice(np.where(y == i)[0], nmin, replace=False) for i in range(len(MODELS))])
rng0.shuffle(sel)
print(f"balancing classes to {nmin} each (from {[int((y==i).sum()) for i in range(len(MODELS))]})")
X = X[sel]; y = y[sel]
D = pickle.load(open(f"{P}/results/Q1/abc/observed.pkl","rb"))
obs = D["obs"]
print(f"training set: {X.shape[0]} simulations x {X.shape[1]} statistics")
for i, m in enumerate(MODELS): print(f"   {m}: {(y==i).sum()}")

print("\n=== PRIOR-PREDICTIVE CHECK: is the observed point inside the simulated cloud? ===")
print(f"{'statistic':<28} {'observed':>10} {'sim 2.5%':>10} {'sim 50%':>10} {'sim 97.5%':>10} {'pctile':>8}  flag")
bad = 0
for j, s in enumerate(STAT):
    q = np.nanpercentile(X[:, j], [2.5, 50, 97.5])
    pc = 100.0*np.mean(X[:, j] <= obs[j])
    flag = "" if 2.5 <= pc <= 97.5 else "  <-- OUTSIDE"
    if flag: bad += 1
    print(f"{s:<28} {obs[j]:>10.4f} {q[0]:>10.4f} {q[1]:>10.4f} {q[2]:>10.4f} {pc:>7.1f}%{flag}")
print(f"\n  {bad}/{len(STAT)} statistics fall outside the 95% prior-predictive range.")
if bad: print("  WARNING: model choice at this point is EXTRAPOLATION. Treat with suspicion.")

print("\n=== RF MODEL CHOICE ===")
rf = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=1)
rf.fit(X, y)
oob = rf.oob_decision_function_
pred_oob = oob.argmax(1)
print(f"OOB accuracy: {(pred_oob==y).mean():.3f}   (chance = {1/len(MODELS):.3f})")
print("\nOOB confusion matrix (row = truth, col = predicted):")
print(f"{'':>6}" + "".join(f"{m:>8}" for m in MODELS))
for i, m in enumerate(MODELS):
    row = [(pred_oob[y==i]==j).sum() for j in range(len(MODELS))]
    tot = max(sum(row),1)
    print(f"{m:>6}" + "".join(f"{r:>8}" for r in row) + f"   ({100*row[i]/tot:.1f}% correct)")

print("\nvariable importance:")
for j in np.argsort(rf.feature_importances_)[::-1]:
    print(f"   {STAT[j]:<28} {rf.feature_importances_[j]:.4f}")

votes = rf.predict_proba(obs.reshape(1,-1))[0]
sel = int(votes.argmax())
print(f"\nselected model: {MODELS[sel]}")
print("vote shares (NOT posterior probabilities):")
for i, m in enumerate(MODELS): print(f"   {m}: {votes[i]:.3f}")

# Pudlo posterior probability
err = (pred_oob != y).astype(float)
rg = RandomForestRegressor(n_estimators=1000, n_jobs=-1, random_state=2)
rg.fit(X, err)
post = 1.0 - float(rg.predict(obs.reshape(1,-1))[0])
print(f"\nposterior probability of {MODELS[sel]} (Pudlo regression estimate): {post:.3f}")
pickle.dump(dict(models=MODELS, stat=STAT, obs=obs, votes=votes, selected=MODELS[sel],
                 posterior=post, oob_acc=float((pred_oob==y).mean()), outside=bad,
                 importance=rf.feature_importances_),
            open(f"{P}/results/Q1/abc/abcrf_result.pkl","wb"))
print(f"\nsaved -> results/Q1/abc/abcrf_result.pkl")
