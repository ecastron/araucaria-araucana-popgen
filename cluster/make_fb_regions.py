"""Regions for the reference-based FreeBayes run (Figure S3, S2D): rad_regions.bed,
rad_regions_clamped.bed and the 200 fb_chunks/chunk_NNN files.
Input : results/sv/loci_vs_genome.paf (minimap2, RAD loci against AraucariaDB v1; map_loci_to_genome.sh)
Rule  : primary alignments with query coverage >= 0.8 and identity >= 0.9 (the filters of
        map_loci_to_genome.sh), target interval padded 50 bp each side, overlapping intervals merged.
Clamp : interval end limited to the contig length.
Chunks: the clamped BED split with GNU `split -n l/200 -d -a 3` (200 line-aligned pieces of near-equal size).
Usage : python make_fb_regions.py OUTDIR"""
import sys, os, collections
P="/home/ecastron/araucaria/results/sv/"; OUT=sys.argv[1] if len(sys.argv)>1 else "."; os.makedirs(OUT,exist_ok=True)
iv=collections.defaultdict(list); clen={}
for ln in open(P+"loci_vs_genome.paf"):
    f=ln.rstrip("\n").split("\t"); clen[f[5]]=int(f[6])
    tp=[x for x in f[12:] if x.startswith("tp:A:")][0][5]
    if tp!="P" or (int(f[3])-int(f[2]))/int(f[1])<0.8 or int(f[9])/int(f[10])<0.9: continue
    iv[f[5]].append((max(0,int(f[7])-50),int(f[8])+50))
reg=[]
for t in sorted(iv):
    cur=None
    for s,e in sorted(iv[t]):
        if cur and s<=cur[1]: cur[1]=max(cur[1],e)
        else:
            if cur: reg.append((t,cur[0],cur[1]))
            cur=[s,e]
    reg.append((t,cur[0],cur[1]))
txt=lambda R:"".join(f"{t}\t{s}\t{e}\n" for t,s,e in R)
clamped=[(t,s,min(e,clen[t])) for t,s,e in reg]
open(os.path.join(OUT,"rad_regions.bed"),"w").write(txt(reg)); open(os.path.join(OUT,"rad_regions_clamped.bed"),"w").write(txt(clamped))
# chunking: GNU split into 200 line-aligned pieces of near-equal bytes
import subprocess
cd=os.path.join(OUT,"fb_chunks"); os.makedirs(cd,exist_ok=True)
subprocess.run(["split","-n","l/200","-d","-a","3",os.path.join(OUT,"rad_regions_clamped.bed"),os.path.join(cd,"chunk_")],check=True)
