#!/usr/bin/env python3
"""BED of RAD loci with exactly ONE confident genomic placement.

Restricting to single-copy loci is what makes read depth interpretable as
haplotype dosage: at a multi-copy locus, depth confounds copy number with
dosage and the test would be meaningless.
"""
import sys
paf, copies, out = sys.argv[1], sys.argv[2], sys.argv[3]
single = {l.split('\t')[0] for l in open(copies) if l.rstrip('\n').split('\t')[1] == '1'}
best = {}
for l in open(paf):
    f = l.rstrip('\n').split('\t')
    q = f[0]
    if q not in single:
        continue
    qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
    nm, bl = int(f[9]), int(f[10])
    if (qe - qs) / qlen < 0.8 or nm / bl < 0.90:
        continue
    if q not in best or nm > best[q][3]:
        best[q] = (f[5], int(f[7]), int(f[8]), nm)
with open(out, 'w') as fh:
    for q, (c, s, e, _) in sorted(best.items(), key=lambda kv: (kv[1][0], kv[1][1])):
        fh.write(f"{c}\t{s}\t{e}\t{q}\n")
print("single-copy loci with a placement: %d" % len(best))
