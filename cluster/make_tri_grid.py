#!/usr/bin/env python
"""Triangular lattice shapefile for feems, over a BUFFERED CONVEX HULL of the sample coords.

feems' create_tile_dict hardcodes `for i in range(3)`, so tiles MUST be triangles.
Rows are offset by half a step (hexagonal packing) so triangles are near-equilateral in
locally-projected km, not in raw degrees.
Usage: make_tri_grid.py <coords.txt> <out.shp> <spacing_km> <buffer_deg>
"""
import sys, math, numpy as np, fiona
from fiona.crs import from_epsg
from shapely.geometry import MultiPoint, Polygon, mapping

coords_f, out_shp, km, buf = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
C = np.loadtxt(coords_f)
hull = MultiPoint([(x, y) for x, y in C]).convex_hull.buffer(buf)
minx, miny, maxx, maxy = hull.bounds
latm = float(np.mean(C[:, 1]))
dlon = km / (111.320 * math.cos(math.radians(latm)))   # horizontal step
dlat = km * math.sqrt(3) / 2 / 110.574                 # row step -> equilateral
nx = int((maxx - minx) / dlon) + 3
ny = int((maxy - miny) / dlat) + 3
def P(i, j):
    return (minx - dlon + i * dlon + (j % 2) * dlon / 2.0, miny - dlat + j * dlat)

tris = []
for j in range(ny):
    for i in range(nx):
        if j % 2 == 0:
            t1 = [P(i, j), P(i + 1, j), P(i, j + 1)]
            t2 = [P(i + 1, j), P(i + 1, j + 1), P(i, j + 1)]
        else:
            t1 = [P(i, j), P(i + 1, j), P(i + 1, j + 1)]
            t2 = [P(i, j), P(i + 1, j + 1), P(i, j + 1)]
        for t in (t1, t2):
            poly = Polygon(t)
            if poly.is_valid and hull.intersects(poly):
                tris.append(poly)

schema = {"geometry": "Polygon", "properties": {"id": "int"}}
with fiona.open(out_shp, "w", driver="ESRI Shapefile", crs=from_epsg(4326), schema=schema) as dst:
    for k, t in enumerate(tris):
        dst.write({"geometry": mapping(t), "properties": {"id": k}})
verts = {(round(x,3), round(y,3)) for t in tris for x, y in t.exterior.coords}
print(f"{out_shp}: spacing {km} km | {len(tris)} triangles | {len(verts)} unique vertices "
      f"| bbox lon {minx:.2f}..{maxx:.2f} lat {miny:.2f}..{maxy:.2f}")
