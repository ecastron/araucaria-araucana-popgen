#!/usr/bin/env python3
"""Stand coordinates for Figures 1 and S1: mean longitude and latitude of the 163 kept trees per stand
(metadata/sra_runs_radseq.csv; figdata/keep163.txt, a copy of results/Q1/feems/keep.txt written by
cluster/make_keep163.py). Coastal stands: VA, TR, PN.
Usage: python figures/stand_coords.py [OUT]   (default: figdata/stand_coords.npz)"""
import sys, numpy as np, pandas as pd
OUT = sys.argv[1] if len(sys.argv) > 1 else "figdata/stand_coords.npz"
m = pd.read_csv("metadata/sra_runs_radseq.csv").set_index("run_accession")
ids = [l.split()[0] for l in open("figdata/keep163.txt") if l.strip()]
g = m.loc[ids].groupby("site_code")[["longitude", "latitude"]].mean().sort_index()
np.savez(OUT, coords=g.values, stands=np.array(g.index), coastal=np.isin(g.index, ["VA", "TR", "PN"]))
