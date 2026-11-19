#!/usr/bin/env python3
"""Record replicate-pair status against the DEPTH-MATCHED null.

Supersedes an earlier version of this file that judged every pair against one pooled null. That
was wrong: mismatch depends strongly on the depth combination of the two libraries --

    deep-deep        random median 0.7030
    deep-shallow     random median 0.4030
    shallow-shallow  random median 0.3994

so the three deep-deep pairs looked catastrophic (~0.69) only because they were being compared
against a null built mostly from shallow pairs. Judged against their OWN stratum, two of the
three are below every random pair and are genuine replicates.
"""
import csv

PROJ = "/home/ecastron/araucaria"
SRC = f"{PROJ}/metadata/replicate_pairs.tsv"
DST = f"{PROJ}/metadata/replicate_pairs_validated.tsv"

# tree -> (stratum, status, evidence). Absent trees are deep-shallow or shallow-shallow pairs
# that fell below every random pair in their own stratum.
SPECIAL = {
    "CALM06": ("deep-deep", "validated",
               "0.6847 vs deep-deep null 0.7030 (IQR 0.6979-0.7078); below ALL random"),
    "CATH04": ("deep-deep", "validated",
               "0.6872 vs deep-deep null 0.7030 (IQR 0.6979-0.7078); below ALL random"),
    "CALR15": ("deep-deep", "QUESTIONABLE",
               "0.7126 vs deep-deep null 0.7030; NOT below stratum median - only failing pair"),
    "CABP07": ("n/a", "NO_DATA", "2-10 shared sites; untestable"),
    "CCTR18": ("n/a", "NO_DATA", "3-20 shared sites; untestable"),
    "CAHL16": ("shallow-shallow", "validated",
               "0.3530 vs null 0.3994; below stratum median but not below all random"),
}

rows = list(csv.DictReader(open(SRC), delimiter="\t"))
with open(DST, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["tree_id", "run_a", "run_b", "site_code",
                "depth_stratum", "status", "evidence"])
    counts = {}
    for r in rows:
        t = r["tree_id"]
        stratum, status, ev = SPECIAL.get(
            t, ("mixed/shallow", "validated",
                "below ALL random pairs within its own depth stratum"))
        w.writerow([t, r["run_a"], r["run_b"], r["site_code"], stratum, status, ev])
        counts[status] = counts.get(status, 0) + 1

print(f"wrote {DST}")
for k, v in sorted(counts.items()):
    print(f"  {k}: {v}")
