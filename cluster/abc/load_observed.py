#!/usr/bin/env python
"""Compute the observed summary-statistic vector, and cache the missingness mask + design
so simulations can be given the SAME missingness pattern (otherwise the observed point sits
outside the simulated cloud for trivial reasons)."""
import sys, csv, re, numpy as np, pickle
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
from sumstats import summarise, vec, STAT_ORDER
P="/home/ecastron/araucaria"
V=f"{P}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
COAST={"VA","TR","PN"}
site={r[0]:r[1] for r in csv.reader(open(f"{P}/metadata/popmap_177.tsv"),delimiter="\t") if len(r)>=2}
keep={l.split()[0] for l in open(f"{P}/results/Q1/feems/keep.txt")}   # relatedness-pruned 163
lat={};lon={}
for r in csv.DictReader(open(f"{P}/metadata/manifest rad runs.csv")):
    try: lat[r["run_accession"]]=float(r["latitude"]); lon[r["run_accession"]]=float(r["longitude"])
    except ValueError: pass
samples=[]; cols=None; rows=[]
for line in open(V):
    if line.startswith("##"): continue
    f=line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples=[re.sub(r"\.[12]$","",s) for s in f[9:]]
        cols=[i for i,s in enumerate(samples) if s in keep and s in site]
        continue
    if len(f[3])!=1 or len(f[4])!=1 or "," in f[4]: continue
    g=[]
    for i in cols:
        t=f[9+i].split(":",1)[0].replace("|","/")
        if t in ("./.","."): g.append(-1)
        else:
            a,b=t.split("/"); g.append(int(a!="0")+int(b!="0"))
    rows.append(g)
G=np.array(rows,dtype=np.int8).T.astype(float)
G[G<0]=np.nan
sel=[samples[i] for i in cols]
stands=sorted({site[s] for s in sel})
deme=np.array([stands.index(site[s]) for s in sel])
coords=np.array([[np.mean([lon[s] for s in sel if site[s]==st]),
                  np.mean([lat[s] for s in sel if site[s]==st])] for st in stands])
coastal=np.array([st in COAST for st in stands])

def load():
    return G, deme, coords, coastal, stands
