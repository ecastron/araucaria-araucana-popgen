#!/usr/bin/env python
"""Clock ABC analysis (17 Nov 2026): model choice and divergence in YEARS, with the validation agreed
before the run. Same model-choice procedure as f4_numbers6 (Pudlo et al. 2016 posterior; Raynal et al.
2019 quantile-weight posterior for parameters), on 12 shape statistics + 3 log clock statistics.

Validation, reported before any number is quoted:
  1. identifiability: OOB R^2 for log T_years, and 95%-interval coverage on 150 held-out simulations
  2. generation time: refit on g < 180 and g >= 180 halves; T_years at the observation should agree
  3. prior floor: posterior mass near the 1,000-yr floor, and refit with training T_years >= 10^4
  4. observed clock values: H_obs <= 0.20 loci (primary) against all loci (sensitivity)
"""
import glob, pickle, json, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
P = "/home/ecastron/araucaria"; OUTD = f"{P}/results/Q1/abc_clock"
LAB = ["H1' Andean source", "H1b' Coastal source", "H2' isolation", "H3' contact"]
Xs, ys, Ps = [], [], []
for f in sorted(glob.glob(f"{OUTD}/s8_*.npz")):
    d = np.load(f, allow_pickle=True)
    if d["X"].size == 0: continue
    Xs.append(d["X"]); ys.append(d["y"]); Ps.append(d["P"])
    STAT = list(d["stat_order"]); PN = list(d["param_names"])
X = np.vstack(Xs); y = np.concatenate(ys); PP = np.vstack(Ps)
shape = pickle.load(open(f"{P}/results/Q1/abc/observed2.pkl", "rb"))["obs"]
CLOCK_PRIMARY = np.log10([0.00269, 0.00275, 0.00296])   # loci with H_obs <= 0.20 (clock_da.py)
CLOCK_ALL     = np.log10([0.00284, 0.00288, 0.00310])   # all loci
obs = np.concatenate([shape, CLOCK_PRIMARY]); obs_all = np.concatenate([shape, CLOCK_ALL])

rng = np.random.default_rng(0)
n = min((y == i).sum() for i in range(4))
sel = np.concatenate([rng.choice(np.where(y == i)[0], n, replace=False) for i in range(4)]); rng.shuffle(sel)
Xb, yb, Pb = X[sel], y[sel], PP[sel]
out = {"n_per_model": int(n), "n_sims": int(Xb.shape[0]), "stat": STAT}

# prior-predictive check on all 15 statistics; statistics the prior cannot reach are dropped
prior_pct = {s: float(100*np.mean(Xb[:, j] <= obs[j])) for j, s in enumerate(STAT)}
bad = [s for s in STAT if not 2.5 <= prior_pct[s] <= 97.5]
USE = [j for j, s in enumerate(STAT) if s not in bad]
out.update(prior_pct=prior_pct, dropped=bad, n_stats=len(USE))
Xu, obsu, obsu_all = Xb[:, USE], obs[USE], obs_all[USE]

def contrast(Xa, ya, labels, o):
    rf = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=1).fit(Xa, ya)
    po = rf.oob_decision_function_.argmax(1)
    v = rf.predict_proba(o.reshape(1, -1))[0]
    rg = RandomForestRegressor(n_estimators=1000, n_jobs=-1, random_state=2).fit(Xa, (po != ya).astype(float))
    return dict(votes=[float(x) for x in v], posterior=1.0-float(rg.predict(o.reshape(1, -1))[0]),
                oob=float((po == ya).mean()), selected=labels[int(v.argmax())])
out["four"]   = contrast(Xu, yb, LAB, obsu)
out["stage1"] = contrast(Xu, (~np.isin(yb, [0, 1])).astype(int), ["one refugium", "two refugia"], obsu)
m = np.isin(yb, [2, 3])
out["stage2"] = contrast(Xu[m], (yb[m] == 3).astype(int), ["isolation", "contact"], obsu)

def qrf(Xa, t, o, seed=3):
    """Raynal et al. 2019: median, 95% interval and weights from leaf co-membership with the observation."""
    rg = RandomForestRegressor(n_estimators=800, oob_score=True, min_samples_leaf=5, n_jobs=-1,
                               random_state=seed).fit(Xa, t)
    lv = rg.apply(Xa); ol = rg.apply(o.reshape(1, -1))[0]
    w = np.zeros(len(t))
    for k in range(lv.shape[1]):
        sm = lv[:, k] == ol[k]; w[sm] += 1.0/sm.sum()
    w /= w.sum(); od = np.argsort(t); cw = np.cumsum(w[od])
    q = lambda pp: float(t[od][min(np.searchsorted(cw, pp), len(t)-1)])
    return dict(med=q(0.5), lo=q(0.025), hi=q(0.975), r2=float(rg.oob_score_)), w, rg

# parameters under the two-refugia class (both H2 and H3 define T1 as the range split)
X2, P2 = Xu[m], Pb[m]
col = lambda k: P2[:, PN.index(k)]
est = {}
for k in ("T_years", "T1", "tau", "Ne", "Nanc", "mu_yr", "g"):
    r, w, _ = qrf(X2, np.log10(col(k)), obsu)
    est[k] = {kk: (10**vv if kk != "r2" else vv) for kk, vv in r.items()}
    if k == "T_years":
        wT = w
        est[k]["mass_below_10k_yr"] = float(w[col(k) < 1e4].sum())
        est[k]["mass_below_3k_yr"]  = float(w[col(k) < 10**3.5].sum())
out["two_refugia_params"] = est
out["prior_T_years"] = [1e3, 1e7]

# ---- validation ----
V = {}
lT = np.log10(col("T_years"))
# 1. coverage on held-out simulations
hold = rng.choice(len(lT), 150, replace=False); keep = np.setdiff1d(np.arange(len(lT)), hold)
rg = RandomForestRegressor(n_estimators=500, min_samples_leaf=5, n_jobs=-1, random_state=4).fit(X2[keep], lT[keep])
lvK = rg.apply(X2[keep]); lvH = rg.apply(X2[hold]); od = np.argsort(lT[keep]); cov = []; err = []
for i in range(len(hold)):
    w = np.zeros(len(keep))
    for k in range(lvK.shape[1]):
        sm = lvK[:, k] == lvH[i, k]; w[sm] += 1.0/sm.sum()
    cw = np.cumsum(w[od]/w.sum()); tt = lT[keep][od]
    lo, me, hi = (tt[min(np.searchsorted(cw, p), len(tt)-1)] for p in (0.025, 0.5, 0.975))
    cov.append(lo <= lT[hold][i] <= hi); err.append(me - lT[hold][i])
V["coverage95_T_years"] = float(np.mean(cov)); V["median_abs_log10_error"] = float(np.median(np.abs(err)))
# 2. generation-time halves
for nm, mk in (("g_lt_180", col("g") < 180), ("g_ge_180", col("g") >= 180)):
    r, _, _ = qrf(X2[mk], lT[mk], obsu); V[nm] = {kk: (10**vv if kk != "r2" else vv) for kk, vv in r.items()}
# 3. prior floor raised to 10^4 years
mk = lT >= 4.0
r, _, _ = qrf(X2[mk], lT[mk], obsu); V["floor_1e4"] = {kk: (10**vv if kk != "r2" else vv) for kk, vv in r.items()}
# 4. observed clock values from all loci
r, _, _ = qrf(X2, lT, obsu_all); V["clock_all_loci"] = {kk: (10**vv if kk != "r2" else vv) for kk, vv in r.items()}
out["validation"] = V

# posterior-predictive percentiles under the two-refugia class, weighted by the T_years posterior
out["pp_pct"] = {s: float(100*np.sum(wT[X2[:, jj] <= obsu[jj]])) for jj, s in enumerate([STAT[j] for j in USE])}
json.dump(out, open(f"{OUTD}/f4_numbers_clock.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in out.items() if k not in ("stat",)}, indent=1, default=float))
