#!/usr/bin/env python
"""feems effective migration surface on a buffered-hull triangular grid.
Fits at several lambda across three grid densities; only features stable across
BOTH lambda and grid resolution should be believed."""
import os, sys, glob, numpy as np
import feems
from feems import SpatialGraph
from feems.utils import prepare_graph_inputs
from pandas_plink import read_plink
from sklearn.impute import SimpleImputer

D = "/home/ecastron/araucaria/results/Q1/feems"
os.chdir(D)
bim, fam, G = read_plink(os.path.join(D, "geno"))
gen = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(np.array(G).T)
coord = np.loadtxt(os.path.join(D, "coords.txt"))
print(f"genotypes {gen.shape} | coords {coord.shape}\n", flush=True)

for km in (20, 12, 8):
    g = f"grids/tri_{km}km.shp"
    try:
        outer, edges, grid, ipmap = prepare_graph_inputs(coord=coord, ggrid=g,
                                                        translated=False, buffer=0.0)
        sg = SpatialGraph(gen, coord, grid, edges, scale_snps=True)
        occ = len(np.unique(ipmap))
        print(f"=== tri_{km}km: {grid.shape[0]} demes, {edges.shape[0]} edges, "
              f"{occ} occupied by samples ===", flush=True)
        for lamb in (20.0, 2.0, 0.5):
            try:
                sg.fit(lamb=lamb, lamb_q=lamb)
                w = np.array(sg.w); lw = np.log10(w)
                np.save(f"w_tri{km}_lamb{lamb}.npy", w)
                print(f"   lamb={lamb:<5} log10(w): min {lw.min():.2f} max {lw.max():.2f} "
                      f"sd {lw.std():.3f} range {lw.max()-lw.min():.2f}", flush=True)
            except Exception as e:
                print(f"   lamb={lamb} FAILED: {e}", flush=True)
        np.save(f"grid_tri{km}.npy", grid); np.save(f"edges_tri{km}.npy", edges)
        np.save(f"ipmap_tri{km}.npy", ipmap)
    except Exception as e:
        print(f"=== tri_{km}km FAILED: {e} ===", flush=True)
print("\nDONE")
