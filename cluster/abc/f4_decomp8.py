"""Where does the width of the clock-ABC divergence interval come from? (18 Nov 2026)
Refit the T_years posterior on subsets of the existing simulations: H2 only, H3 only, and narrow mu windows.
Subsetting a log-uniform prior gives a log-uniform prior on the narrower range, so no new simulations are needed."""
import sys, json, numpy as np
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
exec(open("/home/ecastron/araucaria/scripts/abc/f4_numbers8.py").read().split("# parameters under the two-refugia class")[0])
X2, P2 = Xu[m], Pb[m]; y2 = yb[m]
col = lambda k: P2[:, PN.index(k)]
lT = np.log10(col("T_years")); lmu = np.log10(col("mu_yr")); res = {}
def run(mk, tag):
    r, w, _ = qrf(X2[mk], lT[mk], obsu)
    res[tag] = dict(n=int(mk.sum()), median=10**r["med"], lo=10**r["lo"], hi=10**r["hi"],
                    width_log10=r["hi"]-r["lo"], P_before_26kyr=float(w[lT[mk] >= np.log10(26.5e3)].sum()))
    print(tag, json.dumps(res[tag]), flush=True)
run(np.ones(len(lT), bool), "both_full_mu")
run(y2 == 2, "H2_isolation_only"); run(y2 == 3, "H3_contact_only")
for lo, hi in ((1e-10, 3e-10), (3e-10, 7e-10), (4.1e-10, 1.4e-9), (7e-10, 2e-9)):
    run((lmu >= np.log10(lo)) & (lmu <= np.log10(hi)), f"mu_{lo:.0e}_{hi:.0e}")
json.dump(res, open(f"{OUTD}/f4_decomp_clock.json", "w"), indent=1)
