#!/usr/bin/env python3
"""19 Nov 2026: correct figdata/het_diagnosis.json with the reproduced S2E data (figdata/s2e_mechanisms.tsv,
from cluster scripts s2e_*.py/.R, job 45389). Replaces three mechanism strings whose numbers had no
producer or used mismatched baselines, and the spread_after_intervention block (superseded by the table, which
Figure S2E now reads directly). The previous file is kept as figdata/het_diagnosis.pre20261119.json."""
import json, shutil, pandas as pd
F="figdata/het_diagnosis.json"; shutil.copy(F,"figdata/het_diagnosis.pre20261119.json")
H=json.load(open(F)); T=pd.read_csv("figdata/s2e_mechanisms.tsv",sep="\t")
row=lambda s: T[T.mechanism.str.startswith(s)].iloc[0]
pb=row("beta-binomial"); fb=row("independent pipeline")
fix={"De novo assembly artifact": f"reproduces, spread {fb.filtered_spread:.3f} vs {fb.baseline_spread:.3f} de novo (same estimator)",
     "Genotype-model misfit": f"spread {pb.filtered_spread:.3f} vs {pb.baseline_spread:.3f} for binomial calls on the same loci; the model widens it",
     "Copy number, measured": "F identical across 1, 2, 3-5 and 6+ copy classes (spread 1.405-1.448; unmapped loci 1.404)"}
for m in H["mechanisms"]:
    if m[0] in fix: m[2]=fix[m[0]]
assert all(any(m[0]==k for m in H["mechanisms"]) for k in fix)
H["spread_after_intervention"]="superseded 19 Nov 2026 by figdata/s2e_mechanisms.tsv (per-route baselines)"
H["_source"]=str(H.get("_source",""))+" | corrected 19 Nov 2026 by scripts_local/het_diagnosis_update.py from s2e_mechanisms.tsv"
json.dump(H,open(F,"w"),indent=1); print("updated", F)
