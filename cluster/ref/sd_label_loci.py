#!/usr/bin/env python3
"""Label each RAD locus by whether it sits inside a segmental duplication.

SDs come from the genome self-alignment (40 chunks, minimap2 -x asm20, blocks
>=1 kb at <=10% divergence, self-overlaps removed). Both sides of every pair are
duplicated sequence, so both are recorded.

This is INDEPENDENT of the earlier copy-number test. That one mapped RAD locus
consensuses back to the genome and counted placements, so it could only see
duplications a 150 bp locus still matches. A locus can map uniquely and still lie
inside a diverged or partially collapsed duplicated block -- which is exactly the
blind spot left open on 28 August.
"""
import glob, gzip, collections, sys

SD = collections.defaultdict(list)
n = 0
for f in sorted(glob.glob("/home/ecastron/araucaria/results/sv/selfaln/self_*.paf.gz")):
    with gzip.open(f, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            q, qs, qe, t, ts, te = p[0], int(p[1]), int(p[2]), p[4], int(p[5]), int(p[6])
            SD[q].append((qs, qe)); SD[t].append((ts, te)); n += 1
print("SD pairs read: %d" % n, flush=True)

def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out
tot = 0
for c in SD:
    SD[c] = merge(SD[c]); tot += sum(e - s for s, e in SD[c])
print("contigs with SD: %d | SD bases: %.2f Gb (%.1f%% of 19.32 Gb)"
      % (len(SD), tot / 1e9, 100 * tot / 19.32e9), flush=True)

import bisect
starts = {c: [s for s, _ in v] for c, v in SD.items()}
def inside(c, pos):
    v = SD.get(c)
    if not v: return False
    i = bisect.bisect_right(starts[c], pos) - 1
    return i >= 0 and v[i][0] <= pos < v[i][1]

out = open("/home/ecastron/araucaria/results/sv/loci_sd.tsv", "w")
out.write("locus\tin_sd\n")
nin = ntot = 0
for line in open("/home/ecastron/araucaria/results/sv/dosage/sc_loci.bed"):
    c, s, e, name = line.rstrip("\n").split("\t")
    hit = inside(c, (int(s) + int(e)) // 2)
    out.write("%s\t%d\n" % (name, int(hit)))
    ntot += 1; nin += hit
out.close()
print("RAD loci labelled: %d | inside an SD: %d (%.1f%%)" % (ntot, nin, 100 * nin / ntot))
print("NOTE: these are loci that map to exactly ONE place in the genome, so any")
print("      SD overlap here is invisible to the copy-number test.")
