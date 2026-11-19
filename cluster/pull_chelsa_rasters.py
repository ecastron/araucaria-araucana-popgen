#!/usr/bin/env python3
"""Windowed CHELSA v2.1 surfaces for the study region, for F1's climate panels.

Same remote rasters pull_chelsa.py samples per tree, but read as a WINDOW over
/vsicurl instead of at points, so we get a continuous surface without
downloading the global GeoTIFF. Region is the F1 map extent.

Three variables, chosen on evidence rather than convention:
  bio4  temperature seasonality   -- largest Coastal/Andean separation (d = -2.05)
  bio7  temperature annual range  -- the axis McIntosh et al. (2025) link to the
                                     trait divergence; separates our stands with
                                     no overlap between stand means (d = -1.66)
  bio12 annual precipitation      -- the contrast a reader expects; included with
                                     the caveat that stand ranges DO overlap here
"""
import numpy as np, rasterio
from rasterio.windows import from_bounds

BASE = ("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/"
        "1981-2010/bio/CHELSA_bio{n}_1981-2010_V.2.1.tif")
# The FINAL variable set is bio5 (summer heat), bio18 (growing-season water) and
# bio6 (winter cold) — bio4, bio7 and bio12 were dropped as redundant (|r| = 0.80,
# 0.80 and 0.96 with bio6/bio18). Temperature variables are K*10, so offset -273.15.
SCALE = {5: (0.1, -273.15), 6: (0.1, -273.15), 18: (0.1, 0.0)}
W, E, S, N = -74.4, -70.6, -40.1, -36.9

out = {}
for n, (sc, off) in SCALE.items():
    url = BASE.format(n=n)
    with rasterio.open(f"/vsicurl/{url}") as src:
        win = from_bounds(W, S, E, N, src.transform)
        a = src.read(1, window=win).astype("float32")
        nod = src.nodata
        if nod is not None:
            a[a == nod] = np.nan
        a = a * sc + off
        t = src.window_transform(win)
    out[f"bio{n}"] = a
    out[f"bio{n}_extent"] = np.array([t.c, t.c + a.shape[1] * t.a,
                                      t.f + a.shape[0] * t.e, t.f])
    print("bio%-3d %s  range %.1f to %.1f" % (n, a.shape, np.nanmin(a), np.nanmax(a)))
np.savez_compressed("figdata/chelsa_surfaces.npz", **out)
print("written figdata/chelsa_surfaces.npz")
