#!/usr/bin/env python3
"""Soil covariates per stand from SoilGrids 2.0, via the ISRIC REST point API.

SoilGrids is 250 m, giving a median of 4 distinct cells per stand for 8-10 trees, so it is
usable between stands but only marginal within them. Queried at stand centroids, which is the
scale it supports and the scale the environmental positive control needs.

Same variable family Varas-Myrik et al. (2022) used, so the covariate set is already accepted
by reviewers of the prior paper on this system.
"""
import csv, json, time, urllib.request, statistics as st

D = "/Users/ecastron/Dropbox/01_papers/ARAUCARIA"
OUT = "environment/soil_per_stand.tsv"
API = "https://rest.isric.org/soilgrids/v2.0/properties/query"
PROPS = ["phh2o", "soc", "clay", "sand", "cec", "nitrogen", "bdod"]
DEPTHS = ["0-5cm", "5-15cm"]


def centroids():
    rows = list(csv.DictReader(open(f"{D}/design trees.csv")))
    out = {}
    for r in rows:
        lat = r["mic_latitude"] or r["rad_latitude"]
        lon = r["mic_longitude"] or r["rad_longitude"]
        if not lat:
            continue
        out.setdefault(r["site_code"], []).append((float(lat), float(lon)))
    return {s: (st.mean([p[0] for p in v]), st.mean([p[1] for p in v])) for s, v in out.items()}


def query(lat, lon):
    q = f"{API}?lon={lon:.5f}&lat={lat:.5f}"
    for p in PROPS:
        q += f"&property={p}"
    for d in DEPTHS:
        q += f"&depth={d}"
    q += "&value=mean"
    for attempt in range(4):
        try:
            return json.loads(urllib.request.urlopen(q, timeout=180).read().decode())
        except Exception:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))


def main():
    cen = centroids()
    print(f"stands: {len(cen)}")
    rows = []
    for s, (lat, lon) in sorted(cen.items()):
        d = query(lat, lon)
        rec = {"site_code": s, "lat": round(lat, 5), "lon": round(lon, 5)}
        for layer in d.get("properties", {}).get("layers", []):
            name = layer["name"]
            # d_factor converts SoilGrids mapped units to conventional units
            df = layer["unit_measure"]["d_factor"]
            for dep in layer["depths"]:
                v = dep["values"].get("mean")
                rec[f"{name}_{dep['label']}"] = None if v is None else round(v / df, 3)
        rows.append(rec)
        print(f"  {s}: {len(rec)-3} values")
        time.sleep(1.5)
    cols = ["site_code", "lat", "lon"] + [k for k in rows[0] if k not in ("site_code", "lat", "lon")]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
