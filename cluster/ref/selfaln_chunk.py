#!/usr/bin/env python3
"""Partition contigs into N chunks of roughly equal total bp (greedy, largest first)."""
import sys
fai, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
ctg = sorted(((l.split('\t')[0], int(l.split('\t')[1])) for l in open(fai)),
             key=lambda x: -x[1])
bins = [[] for _ in range(n)]; load = [0]*n
for name, ln in ctg:                      # greedy bin packing
    i = load.index(min(load)); bins[i].append(name); load[i] += ln
for i, b in enumerate(bins, 1):
    with open(f"{out}/chunk_{i:03d}.txt", "w") as fh:
        fh.write("\n".join(b) + "\n")
print("chunks: %d | bp per chunk min %.2f Gb max %.2f Gb"
      % (n, min(load)/1e9, max(load)/1e9))
