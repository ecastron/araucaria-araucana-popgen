#!/usr/bin/env python3
"""CHELSA V2.1 bioclimatic variables per tree, read remotely without downloading.

CHELSA is 30 arc-sec (~0.93 km), which gives a median of 2 distinct cells per stand for 8-10
trees. That is ample for the between-stand environmental comparison and useless within stands;
within-stand microsite is handled separately by terrain derived from a 30 m DEM.

Each raster is opened once over /vsicurl and all points sampled from it, so the whole pull is
19 remote opens rather than 19 x 199 point queries. Raw integer values are written alongside
converted ones, with the scale and offset recorded per variable, so the conversion is auditable
rather than buried.
"""
import csv, os, sys
import rasterio

BASE = ("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/"
        "1981-2010/bio/CHELSA_bio{n}_1981-2010_V.2.1.tif")

# CHELSA V2.1 technical specification. Temperature variables are stored in K*10.
SCALE_OFFSET = {
    1: (0.1, -273.15, "degC"),   2: (0.1, 0.0, "degC"),    3: (0.1, 0.0, "percent"),
    4: (0.1, 0.0, "degC"),       5: (0.1, -273.15, "degC"), 6: (0.1, -273.15, "degC"),
    7: (0.1, 0.0, "degC"),       8: (0.1, -273.15, "degC"), 9: (0.1, -273.15, "degC"),
    10: (0.1, -273.15, "degC"), 11: (0.1, -273.15, "degC"), 12: (0.1, 0.0, "mm"),
    13: (0.1, 0.0, "mm"),       14: (0.1, 0.0, "mm"),      15: (0.1, 0.0, "unitless"),
    16: (0.1, 0.0, "mm"),       17: (0.1, 0.0, "mm"),      18: (0.1, 0.0, "mm"),
    19: (0.1, 0.0, "mm"),
}

DESIGN = "/home/ecastron/araucaria/metadata/design trees.csv"
OUT = "/home/ecastron/araucaria/metadata/chelsa_per_tree.tsv"


def trees():
    out = []
    for r in csv.DictReader(open(DESIGN)):
        lat = r["mic_latitude"] or r["rad_latitude"]
        lon = r["mic_longitude"] or r["rad_longitude"]
        if lat:
            out.append((r["tree_id"], r["site_code"], float(lon), float(lat)))
    return out


def main():
    t = trees()
    pts = [(lon, lat) for _, _, lon, lat in t]
    print(f"trees: {len(t)}", flush=True)
    cols = {}
    for n in range(1, 20):
        url = BASE.format(n=n)
        for attempt in range(3):
            try:
                with rasterio.open(f"/vsicurl/{url}") as src:
                    vals = [v[0] for v in src.sample(pts)]
                break
            except Exception as e:
                if attempt == 2:
                    print(f"  bio{n}: FAILED {e}", file=sys.stderr)
                    vals = [None] * len(pts)
        cols[n] = vals
        print(f"  bio{n} done", flush=True)

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        hdr = ["tree_id", "site_code", "lon", "lat"]
        for n in range(1, 20):
            hdr += [f"bio{n}_raw", f"bio{n}"]
        w.writerow(hdr)
        for i, (tid, site, lon, lat) in enumerate(t):
            row = [tid, site, f"{lon:.5f}", f"{lat:.5f}"]
            for n in range(1, 20):
                raw = cols[n][i]
                sc, off, _ = SCALE_OFFSET[n]
                row += ["" if raw is None else raw,
                        "" if raw is None else f"{raw * sc + off:.3f}"]
            w.writerow(row)
    with open(OUT.replace(".tsv", "_units.tsv"), "w") as fh:
        fh.write("variable\tscale\toffset\tunit\n")
        for n, (sc, off, u) in SCALE_OFFSET.items():
            fh.write(f"bio{n}\t{sc}\t{off}\t{u}\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
