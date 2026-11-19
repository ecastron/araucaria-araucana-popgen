"""Export the balanced reference tables (sim9 primary, sim8 broad-mu) as CSV for abcrf in R.
Same balancing (seed 0) and statistic set (prior-predictive filter) as f4_numbers9/8.py."""
import numpy as np, pandas as pd
for script, tag in (("f4_numbers9.py", "sim9"), ("f4_numbers8.py", "sim8")):
    g = {}; exec(open(f"/home/ecastron/araucaria/scripts/abc/{script}").read().split("# ---- validation ----")[0], g)
    STAT = [g["STAT"][j] for j in g["USE"]]
    df = pd.DataFrame(g["Xu"], columns=STAT); df.insert(0, "model", ["H1", "H1b", "H2", "H3"][0:0] or [ ["H1","H1b","H2","H3"][i] for i in g["yb"] ])
    for k in g["PN"]: df[f"p_{k}"] = g["Pb"][:, g["PN"].index(k)]
    df.to_csv(f"{g['OUTD']}/reftable_{tag}.csv", index=False)
    pd.DataFrame([g["obsu"], g["obsu_all"]], columns=STAT, index=["primary", "all_loci"]).to_csv(f"{g['OUTD']}/observed_{tag}.csv")
    print(tag, df.shape)
