"""Posterior mass of the clock-ABC divergence (years) relative to glacial benchmarks; mode and HDI.
Same data, balancing and statistics as f4_numbers8.py; quantile-weight posterior (Raynal et al. 2019)."""
import sys, json, numpy as np
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
src = open("/home/ecastron/araucaria/scripts/abc/f4_numbers8.py").read()
exec(src.split("# ---- validation ----")[0])       # rebuild X2, obsu, qrf exactly as in the main analysis
lT = np.log10(col("T_years")); res = {}
def summ(Xa, t, o, tag):
    r, w, _ = qrf(Xa, t, o)
    grid = np.linspace(3, 7, 2000)
    dens = (w[None, :]*np.exp(-0.5*((grid[:, None]-t[None, :])/0.06)**2)).sum(1)
    od = np.argsort(dens)[::-1]; cs = np.cumsum(dens[od])/dens.sum(); h = grid[od[:np.searchsorted(cs, 0.95)+1]]
    res[tag] = dict(median=10**r["med"], lo=10**r["lo"], hi=10**r["hi"], mode=10**grid[dens.argmax()],
                    hdi95=[10**h.min(), 10**h.max()],
                    **{f"P_below_{int(b/1000)}kyr": float(w[t < np.log10(b)].sum()) for b in (10e3, 19e3, 26.5e3, 115e3)})
summ(X2, lT, obsu, "primary")
mk = lT >= 4.0; summ(X2[mk], lT[mk], obsu, "floor_1e4")
summ(X2, lT, obsu_all, "clock_all_loci")
for s in (3, 11, 29): r, _, _ = qrf(X2, lT, obsu, seed=s); res[f"forest_seed_{s}_median"] = 10**r["med"]
json.dump(res, open(f"{OUTD}/f4_lgm_clock.json", "w"), indent=1)
print(json.dumps(res, indent=1))
