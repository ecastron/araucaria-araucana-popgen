#!/usr/bin/env python3
"""Sensitivity of the Coastal-Andean environmental contrast (Figure 1B) to single stands.

Revision (25 Nov 2026). With three Coastal stands, one
of them the isolated VA remnant, is any d driven by a single stand?

Inputs (repository environment/ folder, same files as env_contrast_all.py):
  chelsa_per_tree.tsv  per-tree CHELSA v2.1 values, averaged within stand
  soil_per_stand.tsv   SoilGrids 2.0 point values per stand
  terrain_per_tree.tsv per-tree SRTM elevation, averaged within stand (reported, not tested)

Method (same d as the paper): d = (Coastal mean - Andean mean) / population SD of the 12
stand means. For each of the seven retained variables:
  * full d (must reproduce Figure 1B)
  * leave-one-stand-out d for each of the 12 stands; reported: d without each Coastal stand,
    and the min/max over all 12 deletions
  * stand bootstrap: resample 3 Coastal and 9 Andean stands with replacement within range,
    10,000 draws, seed 1; 95% percentile interval of d (draws with zero total SD dropped)

Outputs (figdata/revision/):
  env_loo.tsv           one row per variable
  env_stand_values.tsv  raw stand values of the seven variables plus elevation (Table S6)
"""
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]   # repository root
ENV = REPO / "environment"
OUT = REPO / "figdata/revision"
OUT.mkdir(parents=True, exist_ok=True)
COASTAL = ["VA", "TR", "PN"]

chelsa = pd.read_csv(ENV / "chelsa_per_tree.tsv", sep="\t")
soil = pd.read_csv(ENV / "soil_per_stand.tsv", sep="\t").set_index("site_code")
terr = pd.read_csv(ENV / "terrain_per_tree.tsv", sep="\t")

# the seven variables of Figure 1B, with the labels used in the paper
VARS = {
    "bio5": "max temperature of warmest month (C)",
    "bio6": "min temperature of coldest month (C)",
    "bio18": "precipitation of warmest quarter (mm)",
    "clay_0-5cm": "clay (%)",
    "cec_0-5cm": "cation exchange capacity (mmol(c)/kg)",
    "phh2o_0-5cm": "pH",
    "nitrogen_0-5cm": "total nitrogen (g/kg)",
}
clim = chelsa.groupby("site_code")[["bio5", "bio6", "bio18"]].mean()
tab = clim.join(soil[["clay_0-5cm", "cec_0-5cm", "phh2o_0-5cm", "nitrogen_0-5cm"]])
tab["elevation_m"] = terr.groupby("site_code")["dem_elev_m"].mean()
tab["n_trees_env"] = chelsa.groupby("site_code").size()
tab["range"] = np.where(tab.index.isin(COASTAL), "Coastal", "Andean")
tab[["lat", "lon"]] = soil[["lat", "lon"]]
assert len(tab) == 12


def d_of(x, is_c):
    sd = x.std()  # population SD (ddof 0), as in the paper
    return (x[is_c].mean() - x[~is_c].mean()) / sd if sd > 0 else np.nan


rng = np.random.default_rng(1)
rows = []
for v, lab in VARS.items():
    x = tab[v].to_numpy(float)
    is_c = tab["range"].eq("Coastal").to_numpy()
    full = d_of(x, is_c)
    loo = {}
    for i, s in enumerate(tab.index):
        keep = np.arange(12) != i
        loo[s] = d_of(x[keep], is_c[keep])
    xc, xa = x[is_c], x[~is_c]
    boots = []
    for _ in range(10_000):
        bc = rng.choice(xc, 3, replace=True)
        ba = rng.choice(xa, 9, replace=True)
        allv = np.r_[bc, ba]
        if allv.std() > 0:
            boots.append((bc.mean() - ba.mean()) / allv.std())
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rows.append({
        "variable": v, "label": lab, "d": round(full, 3),
        **{f"d_without_{s}": round(loo[s], 3) for s in COASTAL},
        "d_loo_min": round(min(loo.values()), 3), "d_loo_max": round(max(loo.values()), 3),
        "boot_lo95": round(lo, 3), "boot_hi95": round(hi, 3),
        "coastal_values": ", ".join(f"{s} {tab.loc[s, v]:.3g}" for s in COASTAL),
        "andean_range": f"{xa.min():.3g} to {xa.max():.3g}",
    })

res = pd.DataFrame(rows)
res.to_csv(OUT / "env_loo.tsv", sep="\t", index=False)
cols = ["range", "lat", "lon", "n_trees_env", "elevation_m"] + list(VARS)
tab[cols].sort_values(["range", "lat"]).round(3).to_csv(OUT / "env_stand_values.tsv", sep="\t")
pd.set_option("display.width", 250)
print(res.drop(columns=["label"]).to_string(index=False))
