"""Figure data for F4 and S1 from the literature-mu clock ABC (sim9), the primary demographic analysis
from 18 Nov 2026. Same balancing, statistics and forests as f4_numbers9.py. Also exports the
model-averaged split posterior from the broad-mu run (sim8) for the sensitivity curve in F4C."""
import json, numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import confusion_matrix
S = "/home/ecastron/araucaria/scripts/abc"
def load(script):
    g = {}; exec(open(f"{S}/{script}").read().split("# ---- validation ----")[0], g); return g
g9 = load("f4_numbers9.py")
Xu, yb, obsu, m, qrf, Pb, PN = (g9[k] for k in ("Xu", "yb", "obsu", "m", "qrf", "Pb", "PN"))
X2, P2 = g9["X2"], g9["P2"]
STAT = [g9["STAT"][j] for j in g9["USE"]]
out, arr = {}, {}

def cls(Xa, ya, seed=1):
    rf = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=seed).fit(Xa, ya)
    po = rf.oob_decision_function_.argmax(1)
    return rf, po, confusion_matrix(ya, po), float((po == ya).mean())
rf4, po4, cm4, acc4 = cls(Xu, yb)
arr.update(cm4=cm4, imp4=rf4.feature_importances_)
out.update(acc4=acc4, recall4=[float(cm4[i, i]/cm4[i].sum()) for i in range(4)])
y1 = (~np.isin(yb, [0, 1])).astype(int); _, _, cm1, acc1 = cls(Xu, y1); arr["cm1"] = cm1; out["acc1"] = acc1
y2 = (yb[m] == 3).astype(int); _, po2, cm2, acc2 = cls(Xu[m], y2); arr["cm2"] = cm2; out["acc2"] = acc2

# detection limit: share of contact simulations called "contact", by migrants per generation across the gap
Nm = P2[:, PN.index("Ne")]*P2[:, PN.index("mcross")]
h3 = y2 == 1; p3 = po2[h3]; Nm3 = Nm[h3]
out["power"] = [[lab, int(k.sum()), float(p3[k].mean())] for lo, hi, lab in
                ((0, .2, "< 0.2"), (.2, 1, "0.2-1"), (1, 5, "1-5"), (5, 1e9, "> 5"))
                for k in [(Nm3 >= lo) & (Nm3 < hi)]]
# migration under contact, with its zero-migration calibration (estimator applied to isolation sims)
tgt = np.log10(Nm3)
r, w, rg = qrf(X2[h3], tgt, obsu)
point = float(rg.predict(obsu.reshape(1, -1))[0])
pred_iso = rg.predict(X2[~h3])
out["migration"] = dict(Nm_point=10**point, Nm_median=10**r["med"], ci=[10**r["lo"], 10**r["hi"]], r2=r["r2"],
                        P_gt_1=float(w[Nm3 > 1].sum()), P_gt_1_prior=float(np.mean(Nm3 > 1)),
                        iso_sims_at_or_above_obs=float(np.mean(pred_iso >= point)))

# split posteriors (log10 years) for the density curves
lT = np.log10(P2[:, PN.index("T_years")])
for tag, mk in (("avg", np.ones(len(lT), bool)), ("iso", ~h3), ("con", h3)):
    _, w, _ = qrf(X2[mk], lT[mk], obsu); arr[f"lT_{tag}"] = lT[mk]; arr[f"w_{tag}"] = w
g8 = load("f4_numbers8.py")
lT8 = np.log10(g8["P2"][:, g8["PN"].index("T_years")])
_, w8, _ = g8["qrf"](g8["X2"], lT8, g8["obsu"]); arr["lT_broad"] = lT8; arr["w_broad"] = w8

arr.update(X=Xu, obs=obsu, stat=np.array(STAT))
out.update(stat=STAT, n_sims=int(len(yb)))
OUTD = g9["OUTD"]
np.savez(f"{OUTD}/fig9.npz", **arr)
json.dump(out, open(f"{OUTD}/fig9.json", "w"), indent=1)
print(json.dumps(out, indent=1))
