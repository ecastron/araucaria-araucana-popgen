#!/usr/bin/env python3
"""Faithful replication of Varas-Myrik et al. (2022) filter chain on our assemblies.

Paper: ipyrad ct 85%, then VCFtools excluding sites with heterozygosity >50%, MAF <5%,
missing >50%; individuals with >50% missing also discarded. Result 134 individuals,
49,122 biallelic SNPs. Table 1 per-stand N is reproduced exactly in ARM B.
"""
import sys, csv, re
from collections import defaultdict
CT=sys.argv[1]
V=f"rad/ipyrad_full/full_ct{CT}_outfiles/full_ct{CT}.vcf"
TABLE1={"BP":12,"HL":7,"LM":11,"LR":14,"MC":8,"MM":10,"PC":9,"PN":14,"RA":14,"TH":9,"TR":11,"VA":15}
site_of={r[0]:r[1] for r in csv.reader(open("metadata/popmap_177.tsv"),delimiter="\t") if len(r)>=2}

def sitestats(gts, cols):
    """returns (n_called, n_het, maf) over the given column indices"""
    called=het=alt=0
    for i in cols:
        g=gts[i]
        if g==3: continue
        called+=1; alt+=g
        if g==1: het+=1
    if called==0: return 0,0,0.0
    p=alt/(2*called)
    return called,het,min(p,1-p)

def code(c):
    g=c.split(":",1)[0].replace("|","/")
    if g in ("./.","."): return 3
    a,_,b=g.partition("/")
    return (a!="0")+(b!="0")

# ---- pass 1: per-site stats over ALL libraries, and per-library called counts ----
samples=[]; per_all=None; nsite=0
keepA=[]                      # site indices passing filters computed on ALL libraries
rows=[]                       # store coded genotypes (compact) for pass 2
for line in open(V):
    if line.startswith("##"): continue
    f=line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples=[re.sub(r"\.[12]$","",s) for s in f[9:]]
        per_all=[0]*len(samples); N=len(samples); ALL=list(range(N)); continue
    if len(f[4].split(","))>1 or len(f[3])!=1 or len(f[4])!=1: continue   # biallelic SNPs only
    nsite+=1
    g=[code(c) for c in f[9:]]
    for i in ALL:
        if g[i]!=3: per_all[i]+=1
    rows.append(bytes(g))

print(f"ct 0.{CT[1:]} | {nsite} biallelic SNP sites x {len(samples)} libraries")

def apply_sites(cols):
    out=[]
    for idx,g in enumerate(rows):
        called,het,maf = sitestats(g, cols)
        if called < 0.5*len(cols): continue          # missing > 50% excluded
        if maf < 0.05: continue                       # MAF < 5% excluded
        if het/called > 0.5: continue                 # heterozygosity > 50% excluded
        out.append(idx)
    return out

# ---- ARM A: paper's stated order — sites first (all libraries), then drop individuals ----
sA=apply_sites(list(range(len(samples))))
callA=[sum(1 for i in sA if rows[i][k]!=3) for k in range(len(samples))]
keptA=[k for k in range(len(samples)) if len(sA) and callA[k]/len(sA) >= 0.5]
print(f"  ARM A (sites first, then individuals): {len(keptA)} individuals x {len(sA)} SNPs")

# ---- ARM B: Table 1 per-stand N (best by overall call rate within stand), then sites ----
bystand=defaultdict(list)
for k,s in enumerate(samples):
    st=site_of.get(s)
    if st: bystand[st].append(k)
sel=[]
for st,n in TABLE1.items():
    c=sorted(bystand.get(st,[]), key=lambda k:-per_all[k])[:n]
    sel.extend(c)
sB=apply_sites(sel)
callB=[sum(1 for i in sB if rows[i][k]!=3) for k in sel]
keptB=sum(1 for v in callB if len(sB) and v/len(sB)>=0.5)
print(f"  ARM B (Table 1 per-stand N first, then sites): {len(sel)} individuals x {len(sB)} SNPs"
      f"   [{keptB}/{len(sel)} would survive the >50%-missing individual filter]")
print(f"  per-stand selected: " + " ".join(f"{st}={len([k for k in sel if site_of.get(samples[k])==st])}"
                                          for st in sorted(TABLE1)))
