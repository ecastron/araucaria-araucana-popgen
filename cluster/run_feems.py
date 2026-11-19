#!/usr/bin/env python
"""feems: effective migration surface -- where is gene flow reduced between stands?
Fit at several smoothing values; trust only features stable across them."""
import sys, os, numpy as np
try:
    import feems
    from feems import SpatialGraph
    from feems.utils import prepare_graph_inputs
except Exception as e:
    print("feems NOT IMPORTABLE:", e); sys.exit(0)
from pandas_plink import read_plink
from sklearn.impute import SimpleImputer

D="/home/ecastron/araucaria/results/Q1/feems"
os.chdir(D)
bim, fam, G = read_plink(os.path.join(D,"geno"))
gen = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(np.array(G).T)
coord = np.loadtxt(os.path.join(D,"coords.txt"))
print(f"genotypes {gen.shape} | coords {coord.shape}")

# locate a bundled discrete global grid
import glob
cand = glob.glob(os.path.join(os.path.dirname(feems.__file__), "data", "*.shp"))
print("bundled grids:", [os.path.basename(c) for c in cand])
if not cand:
    print("no bundled grid shapefile found - cannot build the spatial graph"); sys.exit(0)
grids = sorted(cand, key=lambda p: os.path.basename(p))
for gpath in grids[:3]:
    tag = os.path.basename(gpath).replace(".shp","")
    try:
        outer, edges, grid, _ = prepare_graph_inputs(coord=coord, ggrid=gpath,
                                                     translated=False, buffer=1.0)
        sg = SpatialGraph(gen, coord, grid, edges, scale_snps=True)
        print(f"\n=== grid {tag}: {grid.shape[0]} demes, {edges.shape[0]} edges ===")
        for lamb in (20.0, 2.0, 0.2):
            sg.fit(lamb=lamb)
            w = sg.w
            print(f"  lamb={lamb:<5}  log10(w) range {np.log10(w).min():.2f}..{np.log10(w).max():.2f} "
                  f"| sd {np.log10(w).std():.3f}")
            np.save(os.path.join(D, f"w_{tag}_lamb{lamb}.npy"), w)
        np.save(os.path.join(D, f"grid_{tag}.npy"), grid)
        np.save(os.path.join(D, f"edges_{tag}.npy"), edges)
    except Exception as e:
        print(f"  grid {tag} FAILED: {e}")
print("\nDONE")
