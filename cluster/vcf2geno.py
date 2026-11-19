#!/usr/bin/env python3
"""VCF -> LEA .geno: one line per SNP, one char per individual (0/1/2, 9=missing).
Optional third arg: file of sample IDs to keep (one per line)."""
import sys, re
V,OUT = sys.argv[1], sys.argv[2]
keep=None
if len(sys.argv)>3:
    keep={l.strip() for l in open(sys.argv[3]) if l.strip()}
cols=None; n=0
with open(V) as fh, open(OUT,"w") as o:
    for line in fh:
        if line.startswith("##"): continue
        f=line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            s=[re.sub(r"\.[12]$","",x) for x in f[9:]]
            cols=[i for i,x in enumerate(s) if (keep is None or x in keep)]
            open(OUT+".samples","w").write("\n".join(s[i] for i in cols)+"\n")
            print(f"{len(cols)} individuals kept"); continue
        if len(f[3])!=1 or len(f[4])!=1 or "," in f[4]: continue
        row=[]
        for i in cols:
            g=f[9+i].split(":",1)[0].replace("|","/")
            if g in ("./.","."): row.append("9")
            else:
                a,_,b=g.partition("/"); row.append(str((a!="0")+(b!="0")))
        o.write("".join(row)+"\n"); n+=1
print(f"{n} SNPs written to {OUT}")
