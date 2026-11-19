#!/usr/bin/env python3
"""Table S6: sensitivity of the Coastal-Andean split age to analysis choices.

Reads the abcrf outputs (abcrf/abcrf_analysis.R -> figdata/abcrf_results.json; abcrf/abcrf_logscale.R ->
figdata/abcrf_logscale.json) and writes figdata/split_sensitivity.tsv. No new estimation is done here.
Columns: posterior median and 95% credible interval (kyr), posterior probability that the split predates the
onset of the LGM (26.5 kyr), and out-of-bag R2 of the regression forest (on the log10 scale for that row).
"""
import csv
import json

A = json.load(open("figdata/abcrf_results.json"))
L = json.load(open("figdata/abcrf_logscale.json"))["T_years"]

ROWS = [
    ("Primary: isolation and contact simulations pooled", "split_model_averaged"),
    ("Isolation simulations only", "split_isolation"),
    ("Secondary-contact simulations only", "split_contact"),
    ("Estimation on a log10 scale", None),
    ("Broad substitution-rate prior (1 × 10^−10^ to 2 × 10^−9^)", "split_broad_prior"),
    ("Split-age prior floor raised to 10 kyr", "split_floor_1e4"),
    ("Observed time-scale statistics from all loci", "split_all_loci"),
    ("Generation time below 180 years", "split_g_lt_180"),
    ("Generation time 180 years or more", "split_g_ge_180"),
]


def med(x):
    while isinstance(x, list):
        x = x[0]
    return x


out = []
for label, key in ROWS:
    if key is None:
        r = dict(median=med(L["median"]), lo=L["lo"], hi=L["hi"], p=L["P_gt_26500"], r2=L["r2"])
    else:
        e = A[key]
        r = dict(median=e["median"], lo=e["lo"], hi=e["hi"], p=e["P_older_than_26.5kyr"], r2=e["r2"])
    out.append([label, f"{r['median'] / 1e3:.0f}", f"{r['lo'] / 1e3:.0f} to {r['hi'] / 1e3:,.0f}",
                f"{r['p']:.2f}", f"{r['r2']:.2f}"])

with open("figdata/split_sensitivity.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Analysis", "Median (kyr)", "95% credible interval (kyr)", "P(before LGM onset)", "OOB R2"])
    w.writerows(out)
for row in out:
    print("\t".join(row))
