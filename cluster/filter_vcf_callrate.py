#!/usr/bin/env python3
"""Apply the per-analysis call-rate filter the spec defers to this stage.

The assembly deliberately leaves min_samples_locus permissive (4 of 177 = 2.3%) so each
downstream analysis can take the missingness threshold it needs from one assembly. The structure
test was mistakenly run on the raw output, which retains loci present in as few as 4 libraries --
exactly where non-homologous clusters live. This writes filtered VCFs at several call-rate
floors so the test can be repeated as the design intended.

STREAMS the input. An earlier version buffered every line in a list, which needs >7.7 GB for
this VCF against 3 GB on the login node, so it was killed without producing output. All output
files are written in a single pass.
"""
import sys

V = sys.argv[1]
cuts = sorted(float(x) for x in sys.argv[2:]) or [0.3, 0.5, 0.7]

outs, kept, total, n = {}, {c: 0 for c in cuts}, 0, 0
for c in cuts:
    outs[c] = open(V.replace(".vcf", f".cr{int(c * 100)}.vcf"), "w")

with open(V) as fh:
    for line in fh:
        if line.startswith("#"):
            for c in cuts:
                outs[c].write(line)
            if line.startswith("#CHROM"):
                n = len(line.rstrip("\n").split("\t")) - 9
                print(f"{V}: {n} libraries", flush=True)
            continue
        total += 1
        # Count called genotypes without splitting the whole line into a list of fields:
        # the GT is the first token of each sample column, and missing is './.' or '.'.
        parts = line.split("\t")
        called = 0
        for x in parts[9:]:
            g = x.partition(":")[0]
            if g != "./." and g != ".":
                called += 1
        for c in cuts:
            if called >= c * n:
                outs[c].write(line)
                kept[c] += 1
        if total % 500000 == 0:
            print(f"  ...{total} sites scanned", flush=True)

for c in cuts:
    outs[c].close()
print(f"\n{total} sites scanned")
for c in cuts:
    print(f"  call rate >= {c:.0%}: {kept[c]} sites kept ({100 * kept[c] / total:.2f}%)"
          f" -> {V.replace('.vcf', f'.cr{int(c * 100)}.vcf')}")
