"""Reporting items standard in ABC-RF papers (Fraimout et al. 2017; Chapuis et al. 2020; Collin et al. 2021):
votes and posterior over 10 replicate forests (mean and s.d.), prior error rate, and global and local NMAE for
the split age (years) and, under contact, migrants per generation. sim9 reference table, as f4_numbers9.py."""
import json, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
g = {}; exec(open("/home/ecastron/araucaria/scripts/abc/f4_numbers9.py").read().split("# ---- validation ----")[0], g)
Xu, yb, obsu, m, X2, P2, PN, qrf, OUTD = (g[k] for k in ("Xu","yb","obsu","m","X2","P2","PN","qrf","OUTD"))
o = obsu.reshape(1, -1); out = {}
def rep(Xa, ya, tag):
    V, Pp, E = [], [], []
    for s in range(10):
        rf = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=100+s).fit(Xa, ya)
        po = rf.oob_decision_function_.argmax(1); v = rf.predict_proba(o)[0]; k = int(v.argmax())
        rg = RandomForestRegressor(n_estimators=1000, n_jobs=-1, random_state=200+s).fit(Xa, (po != ya).astype(float))
        V.append(v); Pp.append(1-float(rg.predict(o)[0])); E.append(1-float((po == ya).mean()))
    V = np.array(V)
    out[tag] = dict(votes_mean=V.mean(0).tolist(), votes_sd=V.std(0).tolist(), selected=[int(x) for x in V.argmax(1)],
                    posterior_mean=float(np.mean(Pp)), posterior_sd=float(np.std(Pp)),
                    prior_error_rate=float(np.mean(E)), prior_error_sd=float(np.std(E)))
rep(Xu, yb, "four"); rep(Xu, (~np.isin(yb, [0, 1])).astype(int), "stage1"); rep(Xu[m], (yb[m] == 3).astype(int), "stage2")
def nmae(Xa, t, tag):
    rg = RandomForestRegressor(n_estimators=800, oob_score=True, min_samples_leaf=5, n_jobs=-1, random_state=3).fit(Xa, np.log10(t))
    pred = 10**rg.oob_prediction_; err = np.abs(pred - t)/t
    _, w, _ = qrf(Xa, np.log10(t), obsu)
    out[f"nmae_{tag}"] = dict(global_=float(err.mean()), global_median=float(np.median(err)), local=float((w*err).sum()))
nmae(X2, P2[:, PN.index("T_years")], "T_years")
h3 = yb[m] == 3
nmae(X2[h3], (P2[:, PN.index("Ne")]*P2[:, PN.index("mcross")])[h3], "Nm_contact")
json.dump(out, open(f"{OUTD}/f4_conventions_clock_lit.json", "w"), indent=1); print(json.dumps(out, indent=1))
