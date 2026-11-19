"""Strength and timing of secondary contact under H3 (literature-mu clock ABC, sim9),
plus mode and 95% HDI of the model-averaged split date. Same data, balancing and
statistics as f4_numbers9.py; quantile-weight posteriors (Raynal et al. 2019).
Prior mass is reported next to posterior mass so informativeness is visible."""
import sys, json, numpy as np
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
src = open("/home/ecastron/araucaria/scripts/abc/f4_numbers9.py").read()
exec(src.split("# ---- validation ----")[0])       # rebuilds X2, P2, col, obsu, qrf, yb, m
res = {}

# --- mode and HDI of the model-averaged split (years) ---
lT = np.log10(col("T_years"))
r, w, _ = qrf(X2, lT, obsu)
grid = np.linspace(3, 7, 2000)
dens = (w[None, :]*np.exp(-0.5*((grid[:, None]-lT[None, :])/0.06)**2)).sum(1)
od = np.argsort(dens)[::-1]; cs = np.cumsum(dens[od])/dens.sum(); h = grid[od[:np.searchsorted(cs, 0.95)+1]]
res["split_model_averaged"] = dict(median=10**r["med"], lo=10**r["lo"], hi=10**r["hi"],
                                   mode=10**grid[dens.argmax()], hdi95=[10**h.min(), 10**h.max()],
                                   P_25_to_75pct=[float(10**np.interp(q, np.cumsum(w[np.argsort(lT)]), np.sort(lT))) for q in (0.25, 0.75)])

# --- contact parameters, H3 only ---
h3 = (yb[m] == 3)
Xc, Pc = X2[h3], P2[h3]
c = lambda k: Pc[:, PN.index(k)]
targets = {
    "mcross": c("mcross"),                                  # per-generation crossing rate between adjacent Coastal/Andean demes
    "Nm_cross": c("mcross")*c("Ne"),                        # migrants per generation across the contact
    "Tcontact_years": c("Tcontact")*c("g"),                 # how long ago contact began
    "contact_fraction_of_split": c("Tcontact")/c("T1"),     # share of the post-split time spent in contact
}
thr = {"mcross": (1e-4, 1e-3, 1e-2), "Nm_cross": (0.1, 1.0, 10.0),
       "Tcontact_years": (1e3, 1e4, 2.65e4), "contact_fraction_of_split": (0.05, 0.2, 0.5)}
for k, v in targets.items():
    t = np.log10(v)
    r, w, _ = qrf(Xc, t, obsu)
    res[k] = dict(prior_median=float(np.median(v)), prior_lo=float(np.quantile(v, .025)), prior_hi=float(np.quantile(v, .975)),
                  post_median=10**r["med"], post_lo=10**r["lo"], post_hi=10**r["hi"], r2=r["r2"],
                  **{f"P_above_{b:g}": dict(prior=float(np.mean(v > b)), post=float(w[v > b].sum())) for b in thr[k]})
json.dump(res, open(f"{OUTD}/f4_contact_clock_lit.json", "w"), indent=1)
print(json.dumps(res, indent=1))
