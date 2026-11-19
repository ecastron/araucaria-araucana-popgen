#!/usr/bin/env python3
"""VCF -> TreeMix input: per-population counts of the reference and alternate allele.

TreeMix models the covariance of allele FREQUENCIES among populations, so it needs allele counts
per population per SNP, one column per population, formatted "refCount,altCount".

Sites are kept only where every population has at least one called genotype. TreeMix cannot use a
site with a missing population, and dropping them here is explicit rather than leaving TreeMix to
discard them silently.

Note the VCF names samples with a .1 mate suffix (ipyrad takes the name from the symlink target)
while the popmap uses bare accessions; the mapping strips the suffix.
"""
import csv, re, sys

VCF, POPMAP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

pop_of = {}
for r in csv.reader(open(POPMAP), delimiter="\t"):
    if len(r) >= 2:
        pop_of[r[0]] = r[1]

samples, cols, pops = None, None, None
kept = dropped = 0
with open(VCF) as fh, open(OUT, "w") as out:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            pops = sorted({pop_of[s] for s in samples if s in pop_of})
            cols = {p: [i for i, s in enumerate(samples) if pop_of.get(s) == p] for p in pops}
            out.write(" ".join(pops) + "\n")
            print(f"populations: {len(pops)} ({', '.join(pops)})")
            print("libraries per population: " +
                  ", ".join(f"{p}={len(cols[p])}" for p in pops))
            continue
        cells = f[9:]
        counts, ok = [], True
        for p in pops:
            ref = alt = 0
            for i in cols[p]:
                g = cells[i].partition(":")[0].replace("|", "/")
                if g in ("./.", "."):
                    continue
                a, _, b = g.partition("/")
                for x in (a, b):
                    if x == "0":
                        ref += 1
                    else:
                        alt += 1
            if ref + alt == 0:
                ok = False
                break
            counts.append(f"{ref},{alt}")
        if ok:
            out.write(" ".join(counts) + "\n")
            kept += 1
        else:
            dropped += 1
print(f"sites written: {kept}; dropped for a population with no calls: {dropped}")
