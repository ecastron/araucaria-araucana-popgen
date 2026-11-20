#!/usr/bin/env python3
"""Add the BioSample collection_date to the run-level metadata tables.

Reads the BioSample accessions in sra_runs_radseq.csv and sra_runs_amplicon.csv,
fetches each record from NCBI E-utilities (efetch, db=biosample), and writes a
collection_date column next to latitude/longitude. Dates are copied exactly as
submitted to BioSample. The RAD-Seq BioSamples (PRJNA634877) carry no collection_date; those
trees were sampled in 2018 (E. Castro-Nallar, field records), so missing dates in the RAD-Seq table
are written as 2018. Missing dates in any other table stay NA.

Usage: python3 metadata/add_collection_dates.py   (run from the repo root)
"""
import csv
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id="
TABLES = ["sra_runs_radseq.csv", "sra_runs_amplicon.csv"]
FALLBACK = {"sra_runs_radseq.csv": "2018"}  # year of sampling when BioSample has no date


def fetch_dates(accessions):
    """Return {BioSample accession: collection_date} for the given accessions."""
    dates, acc = {}, sorted(accessions)
    for i in range(0, len(acc), 100):  # E-utilities accepts batches of ids
        url = EFETCH + ",".join(acc[i:i + 100])
        req = urllib.request.Request(url, headers={"User-Agent": "metadata-script/1.0"})
        root = ET.fromstring(urllib.request.urlopen(req, timeout=60).read())
        for bs in root.iter("BioSample"):
            date = next((a.text for a in bs.iter("Attribute")
                         if a.get("attribute_name") == "collection_date"), None)
            dates[bs.get("accession")] = date if date else "NA"
        time.sleep(0.4)  # stay under the 3 requests per second limit
    missing = set(acc) - set(dates)
    assert not missing, f"BioSamples not returned by NCBI: {sorted(missing)[:5]}"
    return dates


for name in TABLES:
    path = HERE / name
    rows = list(csv.DictReader(path.open()))
    dates = fetch_dates({r["biosample"] for r in rows})
    fields = [f for f in rows[0].keys() if f != "collection_date"]
    pos = fields.index("longitude") + 1
    fields.insert(pos, "collection_date")
    for r in rows:
        r["collection_date"] = dates[r["biosample"]]
        if r["collection_date"] == "NA" and name in FALLBACK:
            r["collection_date"] = FALLBACK[name]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    tally = {}
    for r in rows:
        tally[r["collection_date"]] = tally.get(r["collection_date"], 0) + 1
    print(name, dict(sorted(tally.items())))
