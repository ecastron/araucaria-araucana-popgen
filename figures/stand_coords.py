#!/usr/bin/env python3
"""Stand coordinates for Figures 1 and S1: mean longitude and latitude of the 163 kept trees per stand
(metadata/sra_runs_radseq.csv; figdata/keep163.txt, a copy of results/Q1/feems/keep.txt written by
cluster/make_keep163.py). Coastal stands: VA, TR, PN.
Reconstructed 19 Nov 2026 to replace lost inline code; matches figdata/stand_coords.npz to 3e-14.
Usage: python figures/stand_coords.py [OUT]   (default: a new file, never the original)"""
import sys, numpy as np, pandas as pd
OUT = sys.argv[1] if len(sys.argv) > 1 else "stand_coords_rebuilt.npz"
m = pd.read_csv("metadata/sra_runs_radseq.csv").set_index("run_accession")
ids = [l.split()[0] for l in open("figdata/keep163.txt") if l.strip()]
g = m.loc[ids].groupby("site_code")[["longitude", "latitude"]].mean().sort_index()
np.savez(OUT, coords=g.values, stands=np.array(g.index), coastal=np.isin(g.index, ["VA", "TR", "PN"]))
ref = np.load("figdata/stand_coords.npz", allow_pickle=True)
print("max |difference| from figdata/stand_coords.npz:", np.abs(ref["coords"] - g.values).max(),
      "| stands and range labels identical:", list(ref["stands"]) == list(g.index) and (ref["coastal"] == np.isin(g.index, ["VA", "TR", "PN"])).all())
