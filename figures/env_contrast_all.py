"""Coastal versus Andean environmental contrast for all 39 extracted variables.

Produces: figdata/env_contrast_all.tsv (output path given on the command line,
default ./env_contrast_all.tsv).
Columns: layer, variable, coastal_mean, andean_mean, cohen_d, perm_p, separates.

Inputs (environment/ in this repository):
  chelsa_per_tree.tsv   19 CHELSA v2.1 bioclim variables per tree (scaled columns bio1..bio19)
  terrain_per_tree.tsv  6 SRTM terrain variables per tree
  soil_per_stand.tsv    14 SoilGrids 2.0 properties per stand

Method:
  * Unit is the STAND (12 stands). Per-tree climate and terrain values are averaged
    within stand over every tree in the file (199 trees, not only the 163 genotyped).
  * Coastal = VA, TR, PN (3 stands); Andean = the other 9.
  * cohen_d = (Coastal mean - Andean mean) / SD of the 12 stand means, with the
    population SD (ddof = 0, numpy's default). Positive = higher on the Coastal range.
    NOTE: this is a bespoke standardized difference, not the textbook pooled-SD
    Cohen's d (which gives larger values, e.g. clay 0-5 cm 2.90 instead of 1.87).
  * perm_p: exact permutation over all C(12,3) = 220 assignments of 3 stands to the
    Coastal range; statistic |mean difference|; p = share of assignments with a
    statistic >= the observed one (observed included, ties counted), so the minimum is 1/220.
  * separates = 1 when the Coastal and Andean stand-mean ranges do not overlap.
  * Rows sorted by perm_p, then by |cohen_d| descending.

"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ENV = Path(__file__).resolve().parents[1] / "environment"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("env_contrast_all.tsv")
COASTAL = {"VA", "TR", "PN"}

chelsa = pd.read_csv(ENV / "chelsa_per_tree.tsv", sep="\t")
terrain = pd.read_csv(ENV / "terrain_per_tree.tsv", sep="\t")
soil = pd.read_csv(ENV / "soil_per_stand.tsv", sep="\t")

bio = [f"bio{i}" for i in range(1, 20)]
ter = ["dem_elev_m", "slope_deg", "northness", "eastness", "tpi_m", "heat_load"]
layers = {
    "climate": chelsa.groupby("site_code")[bio].mean(),
    "terrain": terrain.groupby("site_code")[ter].mean(),
    "soil": soil.set_index("site_code").drop(columns=["lat", "lon"]),
}

stands = sorted(layers["soil"].index)
assert len(stands) == 12 and COASTAL <= set(stands)
combos = list(itertools.combinations(range(12), 3))
assert len(combos) == 220

rows = []
for layer, df in layers.items():
    df = df.loc[stands]
    is_c = np.array([s in COASTAL for s in stands])
    for v in df.columns:
        x = df[v].to_numpy(float)
        c, a = x[is_c], x[~is_c]
        diff = c.mean() - a.mean()
        d = diff / x.std()  # population SD of the 12 stand means
        # exact permutation of the range labels over stands
        tot = x.sum()
        null = np.array([abs(x[list(k)].mean() - (tot - x[list(k)].sum()) / 9) for k in combos])
        # ties counted as extreme (tolerance absorbs float noise in equal sums)
        p = np.mean(null >= abs(diff) - 1e-12)
        sep = int(c.max() < a.min() or c.min() > a.max())
        rows.append((layer, v, c.mean(), a.mean(), d, p, sep))

out = pd.DataFrame(rows, columns=["layer", "variable", "coastal_mean", "andean_mean",
                                  "cohen_d", "perm_p", "separates"])
out["absd"] = out["cohen_d"].abs()
out = out.sort_values(["perm_p", "absd"], ascending=[True, False], kind="mergesort").drop(columns="absd")
with open(OUT, "w") as fh:
    fh.write("\t".join(out.columns) + "\n")
    for r in out.itertuples(index=False):
        fh.write(f"{r.layer}\t{r.variable}\t{r.coastal_mean:.4f}\t{r.andean_mean:.4f}\t"
                 f"{r.cohen_d:.3f}\t{r.perm_p:.4f}\t{r.separates}\n")
print(f"wrote {OUT} ({len(out)} variables)")
