#!/usr/bin/env python3
"""Summarize the sNMF K = 2 ancestry null (revision, 26 Nov 2026).
Input : figdata/revision/ancnull/snmf_task*.tsv (cluster/revision/ancestry_null.R), selected.tsv
Output: figdata/revision/ancestry_null_summary.json and printed table. For each metric: observed (163 trees,
all SNPs and 20,000-SNP subset), H2 (zero gene flow) and H3 (contact) medians and 90% ranges, and
the percentile of the observed value within H2.
"""
import json, glob
from pathlib import Path
import pandas as pd, numpy as np
D = Path(__file__).resolve().parents[1] / "figdata/revision"
t = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(D / "ancnull/snmf_task*.tsv")))])
t["set"] = t.dataset.str.extract(r"^(H2|H3|obs_all|obs20k)")[0]
sel = pd.read_csv(D / "ancnull/selected.tsv", sep="\t")
met = [c for c in t.columns if c not in ("dataset", "set")]
out = {"n": t.set.value_counts().to_dict(), "metrics": {}}
print(f"{'metric':20s} {'obs_all':>8s} {'obs20k':>8s} {'H2 med (5-95%)':>22s} {'H3 med (5-95%)':>22s} {'obs20k pct in H2':>16s}")
for m in met:
    o1 = float(t.loc[t.set == "obs_all", m].iloc[0]); o2 = float(t.loc[t.set == "obs20k", m].iloc[0])
    h2 = t.loc[t.set == "H2", m].to_numpy(float); h3 = t.loc[t.set == "H3", m].to_numpy(float)
    q = lambda v: [float(np.median(v)), float(np.percentile(v, 5)), float(np.percentile(v, 95))]
    pct = float(100 * np.mean(h2 <= o2))
    out["metrics"][m] = dict(obs_all=o1, obs20k=o2, H2=q(h2), H3=q(h3), obs20k_pct_in_H2=pct,
                             obs_all_pct_in_H2=float(100 * np.mean(h2 <= o1)))
    a, b = q(h2), q(h3)
    print(f"{m:20s} {o1:8.3f} {o2:8.3f} {a[0]:8.3f} ({a[1]:.3f}-{a[2]:.3f}) {b[0]:8.3f} ({b[1]:.3f}-{b[2]:.3f}) {pct:16.0f}")
h3 = sel[sel.model == "H3"]
out["H3_selected_mcross_median"] = float(h3.p_mcross.median())
out["H2_T_years_median"] = float(sel[sel.model == "H2"].p_T_years.median())
print("H3 selected mcross median", out["H3_selected_mcross_median"], " H2 T_years median", out["H2_T_years_median"])
json.dump(out, open(D / "ancestry_null_summary.json", "w"), indent=1)
