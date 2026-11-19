#!/usr/bin/env python3
"""SoilGrids 2.0 surfaces for the study region, for F1 panel B.

EXTENDED 10 Sep 2026: all four soil variables of the final literature set are now pulled
(clay, CEC, pH, nitrogen), so panel B can map all seven model variables rather than four.
Conversion factors follow the SoilGrids 2.0 FAQ: clay g/kg->%%, CEC mmol(c)/kg->cmol(c)/kg,
phh2o pH*10->pH, nitrogen cg/kg->g/kg.

pull_soil.py used the ISRIC point API at 12 stand centroids, which cannot be
mapped. The COGs are served over webdav but in Homolosine (ESRI:54052), so the
window has to be reprojected to WGS84 -- done here with a WarpedVRT rather than
downloading the global raster.
"""
import numpy as np, rasterio
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds
from rasterio.enums import Resampling

B = "/vsicurl/https://files.isric.org/soilgrids/latest/data"
SRC_LIST = [("clay",  f"{B}/clay/clay_0-5cm_mean.vrt",     10.0),
            ("cec",   f"{B}/cec/cec_0-5cm_mean.vrt",       10.0),
            ("ph",    f"{B}/phh2o/phh2o_0-5cm_mean.vrt",   10.0),
            ("nitro", f"{B}/nitrogen/nitrogen_0-5cm_mean.vrt", 100.0)]
W, E, S, N = -74.4, -70.6, -40.1, -36.9

out = {}
for NAME, SRC, DIV in SRC_LIST:
  with rasterio.open(SRC) as src:
    with WarpedVRT(src, crs="EPSG:4326", resampling=Resampling.bilinear) as vrt:
        win = from_bounds(W, S, E, N, vrt.transform)
        a = vrt.read(1, window=win).astype("float32")
        if vrt.nodata is not None:
            a[a == vrt.nodata] = np.nan
        t = vrt.window_transform(win)
  a = a / DIV
  out[NAME] = a
  out[NAME + "_extent"] = np.array([t.c, t.c + a.shape[1] * t.a, t.f + a.shape[0] * t.e, t.f])
  print("%-5s %s range %.1f to %.1f | nan %.1f%%"
        % (NAME, a.shape, np.nanmin(a), np.nanmax(a), 100 * np.isnan(a).mean()))
np.savez_compressed("figdata/soilgrids_soil.npz", **out)
print("written figdata/soilgrids_soil.npz")
