#!/usr/bin/env python
"""Build a MIGRATE-N infile: two populations (Coastal vs Andean), PSEUDO-HAPLOID.

One allele is drawn at random from each called genotype, which avoids depending on the
heterozygote calls (the class shown to be unreliable in this dataset). Missing genotypes are
written as '?'. Loci are chosen for high call rate in BOTH populations so MIGRATE sees dense data.
Usage: make_migrate_infile.py <out.infile> <n_loci> <datatype a|u> [seed]
"""
import sys, csv, re, random
OUT, NLOCI, DT = sys.argv[1], int(sys.argv[2]), sys.argv[3]
random.seed(int(sys.argv[4]) if len(sys.argv)>4 else 2026)
V="rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
COAST={"VA","TR","PN"}
site={r[0]:r[1] for r in csv.reader(open("metadata/popmap_177.tsv"),delimiter="\t") if len(r)>=2}
keep={l.split()[0] for l in open("results/Q1/feems/keep.txt")}      # relatedness-pruned 163

samples=[]; cols=None; loci=[]
for line in open(V):
    if line.startswith("##"): continue
    f=line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples=[re.sub(r"\.[12]$","",s) for s in f[9:]]
        cols=[i for i,s in enumerate(samples) if s in keep and s in site]
        pops={"Coastal":[i for i in cols if site[samples[i]] in COAST],
              "Andean":[i for i in cols if site[samples[i]] not in COAST]}
        continue
    if len(f[3])!=1 or len(f[4])!=1 or "," in f[4]: continue
    g=[]
    for i in cols:
        t=f[9+i].split(":",1)[0].replace("|","/")
        if t in ("./.","."): g.append(None)
        else:
            a,b=t.split("/"); g.append(random.choice([a,b]))
    loci.append((f[0],f[1],g))

idx={c:j for j,c in enumerate(cols)}
def callrate(l, members):
    v=[l[2][idx[m]] for m in members]
    return sum(1 for x in v if x is not None)/len(v)
MINMAF=0.05
scored=[]
for l in loci:
    cC=callrate(l,pops["Coastal"]); cA=callrate(l,pops["Andean"])
    if cC<0.5 or cA<0.5: continue
    obs=[x for x in l[2] if x is not None]
    if not obs: continue
    p1=sum(1 for x in obs if x=="1")/len(obs)
    maf=min(p1,1-p1)
    if maf<MINMAF: continue                      # must be segregating
    scored.append((min(cC,cA)*maf, l))           # rank on call rate x informativeness
scored.sort(key=lambda x:-x[0])
sel=[l for _,l in scored[:NLOCI]]
print(f"{len(loci)} biallelic SNPs -> {len(scored)} pass (call rate >=0.5 in both pops, MAF >=0.05) "
      f"-> {len(sel)} selected")
print(f"populations: Coastal n={len(pops['Coastal'])}, Andean n={len(pops['Andean'])}")

with open(OUT,"w") as o:
    o.write(f"   2 {len(sel)} Araucaria two-cordillera pseudo-haploid ({DT})\n")
    for pname in ("Coastal","Andean"):
        mem=pops[pname]
        o.write(f"{len(mem):>5} {pname}\n")
        for m in mem:
            nm=("%s%03d"%(pname[0],mem.index(m)+1)).ljust(10)
            alle=[]
            for l in sel:
                a=l[2][idx[m]]
                alle.append("?" if a is None else str(int(a)+1))   # alleles as 1/2
            o.write(nm + " " + " ".join(alle) + "\n")
print(f"wrote {OUT}")
