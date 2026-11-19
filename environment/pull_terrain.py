#!/usr/bin/env python3
"""Per-tree terrain covariates from SRTM 30 m, via the OpenTopoData point API.

Within-stand microsite cannot come from a climate database: WorldClim/CHELSA at 30 arc-sec
gives a median of 2 distinct cells per stand for 8-10 trees. A 30 m DEM gives ~8, roughly one
per tree. Finer DEMs buy nothing because coordinate error here is tens of metres, so cell
assignment - not DEM resolution - is the limiting uncertainty.

Slope and aspect are computed by the Horn method from a 3x3 neighbourhood sampled at 30 m
spacing. Aspect is decomposed into northness and eastness because aspect is circular and
cannot enter a linear model directly. Topographic position index (centre minus neighbourhood
mean) separates ridges from hollows. Heat load index follows McCune & Keon (2002), which
folds slope, aspect and latitude into a single insolation proxy.
"""
import csv, json, math, time, urllib.request, urllib.parse

D = "/Users/ecastron/Dropbox/01_papers/ARAUCARIA"
OUT = "environment/terrain_per_tree.tsv"
API = "https://api.opentopodata.org/v1/srtm30m"
CELL = 30.0          # metres; matches SRTM resolution and our positional accuracy
BATCH = 100          # API maximum locations per request


def tree_coords():
    """Prefer microbiome coordinates (used for within-stand spatial terms), else RAD."""
    out = []
    for r in csv.DictReader(open(f"{D}/design trees.csv")):
        if r["mic_latitude"]:
            lat, lon, src = float(r["mic_latitude"]), float(r["mic_longitude"]), "mic"
        elif r["rad_latitude"]:
            lat, lon, src = float(r["rad_latitude"]), float(r["rad_longitude"]), "rad"
        else:
            continue
        out.append({"tree_id": r["tree_id"], "site": r["site_code"], "lat": lat, "lon": lon,
                    "src": src, "recorded_alt": r["altitude_m"]})
    return out


def neighbourhood(lat, lon):
    """3x3 grid, rows north->south, cols west->east."""
    dlat = CELL / 111320.0
    dlon = CELL / (111320.0 * math.cos(math.radians(lat)))
    pts = []
    for i in (1, 0, -1):            # north row first
        for j in (-1, 0, 1):        # west to east
            pts.append((lat + i * dlat, lon + j * dlon))
    return pts


def fetch(points):
    """Batch elevation lookups, respecting the public API's 1 call/second limit."""
    vals = []
    for i in range(0, len(points), BATCH):
        chunk = points[i:i + BATCH]
        loc = "|".join(f"{a:.6f},{b:.6f}" for a, b in chunk)
        url = f"{API}?locations={urllib.parse.quote(loc, safe='|,')}"
        for attempt in range(4):
            try:
                d = json.loads(urllib.request.urlopen(url, timeout=120).read().decode())
                if d.get("status") == "OK":
                    vals.extend(r["elevation"] for r in d["results"])
                    break
                raise RuntimeError(d.get("error", "non-OK status"))
            except Exception as e:
                if attempt == 3:
                    raise
                time.sleep(3 * (attempt + 1))
        time.sleep(1.1)
        print(f"  fetched {min(i+BATCH, len(points))}/{len(points)}", end="\r")
    print()
    return vals


def horn(z):
    """z is the 9-cell neighbourhood, row-major from NW. Returns slope deg, aspect deg."""
    a, b, c, d, e, f, g, h, i = z
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * CELL)
    dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * CELL)
    slope = math.degrees(math.atan(math.hypot(dzdx, dzdy)))
    asp = math.degrees(math.atan2(dzdy, -dzdx))
    asp = (90.0 - asp) % 360.0
    return slope, asp


def heat_load(lat, slope_deg, aspect_deg):
    """McCune & Keon (2002) heat load index; folded aspect about the SW/NE axis."""
    L = math.radians(abs(lat))
    S = math.radians(slope_deg)
    A = math.radians(abs(180.0 - abs(aspect_deg - 225.0)))
    return (-1.467 + 1.582 * math.cos(L) * math.cos(S)
            - 1.500 * math.cos(A) * math.sin(S) * math.sin(L)
            - 0.262 * math.sin(L) * math.sin(S)
            + 0.607 * math.sin(A) * math.sin(S))


def main():
    trees = tree_coords()
    print(f"trees with coordinates: {len(trees)}")
    pts = [p for t in trees for p in neighbourhood(t["lat"], t["lon"])]
    print(f"elevation points to fetch: {len(pts)} ({len(pts)//BATCH + 1} requests)")
    elev = fetch(pts)
    assert len(elev) == len(pts), f"expected {len(pts)} elevations, got {len(elev)}"

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["tree_id", "site_code", "coord_src", "recorded_alt_m", "dem_elev_m",
                    "slope_deg", "aspect_deg", "northness", "eastness", "tpi_m", "heat_load"])
        for k, t in enumerate(trees):
            z = elev[k * 9:(k + 1) * 9]
            if any(v is None for v in z):
                continue
            slope, asp = horn(z)
            centre = z[4]
            tpi = centre - (sum(z) - centre) / 8.0
            w.writerow([t["tree_id"], t["site"], t["src"], t["recorded_alt"], f"{centre:.0f}",
                        f"{slope:.2f}", f"{asp:.1f}",
                        f"{math.cos(math.radians(asp)):.4f}",
                        f"{math.sin(math.radians(asp)):.4f}",
                        f"{tpi:.2f}", f"{heat_load(t['lat'], slope, asp):.4f}"])
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
