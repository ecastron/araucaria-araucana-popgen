#!/usr/bin/env python3
"""Write a dp20 VCF: genotypes with DP < MINDP become ./., then sites called in fewer than
MINCALL libraries are dropped. Uses the VCF DP field (p[1]) -- the same field and semantics
q1_popgen.py already uses for its per-stand diversity table, so the comparison is like-for-like.
Usage: mask_depth.py <in.vcf> <out.vcf> <mindp> <mincall>
"""
import sys
IN, OUT, MINDP, MINCALL = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
kept = dropped = masked = total = 0
per_ind = None
with open(IN) as fh, open(OUT, "w") as out:
    for line in fh:
        if line.startswith("##"):
            out.write(line); continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            per_ind = [0]*(len(f)-9); out.write(line); continue
        total += 1
        called = 0
        for k, c in enumerate(f[9:]):
            p = c.split(":")
            if p[0] in ("./.", "."): continue
            try: dp = int(p[1])
            except (IndexError, ValueError): dp = 0
            if dp < MINDP:
                p[0] = "./."; f[9+k] = ":".join(p); masked += 1
            else:
                called += 1; per_ind[k] += 1
        if called >= MINCALL:
            out.write("\t".join(f)+"\n"); kept += 1
        else:
            dropped += 1
import statistics as st
print(f"sites: {total} in, {kept} kept, {dropped} dropped (<{MINCALL} called at {MINDP}x)")
print(f"genotypes masked to missing: {masked}")
print(f"per-library called sites: median {int(st.median(per_ind))} min {min(per_ind)} max {max(per_ind)}")
print(f"libraries with >=2000: {sum(1 for x in per_ind if x>=2000)}/{len(per_ind)}")
