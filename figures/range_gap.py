#!/usr/bin/env python3
"""Shortest distance between the Coastal and Andean Araucaria araucana forests (CONAF Catastro 2024).

Reads figdata/distribution_cache.pkl (written by cache_distribution.py from the CONAF shapefiles) and
measures the minimum great-circle-equivalent distance between any vertex of a Coastal polygon and any
vertex of an Andean polygon, for two definitions:
  dominant       Araucaria is the forest type (tipo_fores), the definition behind 201,264 ha
  any component  dominant polygons plus those where Araucaria is listed as a component species
Polygons are assigned to a range by mean ring longitude: west of 72.5 W is Coastal (Nahuelbuta),
east is Andean. The longitude histogram has no polygons between 72.95 and 72.15 W except one
component-species ring near 72.6 W, which this rule assigns to the Coastal side; the script also
reports the gap without it. Distances use UTM zone 18S (EPSG:32718) and a KD-tree (scipy).
Also reports the closest sampled cross-range stand pair (PN-TH) for comparison.
Output: figdata/range_gap.json
"""
import json
import pickle
import numpy as np
import pyproj
from scipy.spatial import cKDTree

D = pickle.load(open("figdata/distribution_cache.pkl", "rb"))
tf = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32718", always_xy=True)
SPLIT = -72.5


def split_rings(rings):
    coast, andes = [], []
    for r in rings:
        (coast if r[:, 0].mean() < SPLIT else andes).append(r)
    return coast, andes


def to_xy(rings):
    pts = np.vstack(rings).astype(float)
    x, y = tf.transform(pts[:, 0], pts[:, 1])
    return np.column_stack([x, y]), pts


def min_gap(coast, andes):
    xa, la = to_xy(andes)
    xc, lc = to_xy(coast)
    d, j = cKDTree(xa).query(xc, k=1)
    i = int(np.argmin(d))
    return dict(km=float(d[i] / 1000), coastal_point=lc[i].tolist(), andean_point=la[j[i]].tolist(),
                n_coastal_rings=len(coast), n_andean_rings=len(andes))


out = {}
c_dom, a_dom = split_rings(D["rings_dom"])
out["dominant"] = min_gap(c_dom, a_dom)
c_x, a_x = split_rings(D["rings_extra"])
out["any_component"] = min_gap(c_dom + c_x, a_dom + a_x)
lone = [r for r in c_x if r[:, 0].mean() > -73.0]
if lone:
    out["any_component_without_valley_ring"] = min_gap(
        c_dom + [r for r in c_x if r[:, 0].mean() <= -73.0], a_dom + a_x)
    out["valley_rings_lon"] = [float(r[:, 0].mean()) for r in lone]

S = np.load("figdata/stand_coords.npz", allow_pickle=True)
st = list(S["stands"]) if "stands" in S.files else None
if st is not None:
    sc = S["coords"]
    p, t = sc[st.index("PN")], sc[st.index("TH")]
    xy = np.column_stack(tf.transform([p[0], t[0]], [p[1], t[1]]))
    out["closest_sampled_pair_PN_TH_km"] = float(np.linalg.norm(xy[0] - xy[1]) / 1000)

json.dump(out, open("figdata/range_gap.json", "w"), indent=1)
for k, v in out.items():
    print(k, v if not isinstance(v, dict) else f"{v['km']:.1f} km  coast {v['coastal_point']}  andes {v['andean_point']}")
