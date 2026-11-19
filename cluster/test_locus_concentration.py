#!/usr/bin/env python3
"""Are the opposite-homozygote calls concentrated in a few loci, or spread across all of them?

Established by the previous test: ~80% of disagreements between two libraries of the SAME tree
are opposite homozygotes (0/0 vs 1/1), and 38% remain so at >=30x depth. Dropout cannot produce
those, so the sites being compared are not homologous.

This asks whether the damage is localised:

  CONCENTRATED -> a minority of loci are collapsed repeats/paralogs. Excluding them should leave
                  a subset with low replicate mismatch, and Q6 is recoverable by filtering.
  UNIFORM      -> non-homology is pervasive. No filter recovers a usable subset, and de novo
                  genotyping of this species fails without a reference.

The decisive output is not the concentration statistic itself but the LAST table: residual
mismatch and how many loci survive after removing the worst offenders. A filter that leaves
2,000 clean loci is a usable dataset; one that leaves 40 is not.
"""
import csv, glob, os, re, collections

PROJ = "/home/ecastron/araucaria"
VCF = sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf"))[0]

pairs = [(r["tree_id"], r["run_a"], r["run_b"])
         for r in csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"),
                                 delimiter="\t") if r["status"] == "validated"]

samples, recs = [], []
with open(VCF) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            continue
        recs.append((f[0], f[9:]))          # (locus, genotype cells)
idx = {s: i for i, s in enumerate(samples)}
print(f"{os.path.basename(VCF)}: {len(samples)} libraries, {len(recs)} sites, "
      f"{len(set(l for l, _ in recs))} loci, {len(pairs)} validated pairs\n")

def parse(c):
    p = c.split(":")
    try:
        return p[0].replace("|", "/"), int(p[1])
    except (IndexError, ValueError):
        return p[0].replace("|", "/"), 0

def het(g):
    a, _, b = g.partition("/")
    return a != b

# Per locus: comparable sites, mismatches, and opposite-homozygote count across all pairs.
loc = collections.defaultdict(lambda: [0, 0, 0])      # locus -> [compared, mismatch, opp_hom]
for locus, cells in recs:
    for _, a, b in pairs:
        if a not in idx or b not in idx:
            continue
        ga, da = parse(cells[idx[a]])
        gb, db = parse(cells[idx[b]])
        if ga in ("./.", ".") or gb in ("./.", "."):
            continue
        s = loc[locus]
        s[0] += 1
        if ga != gb:
            s[1] += 1
            if not het(ga) and not het(gb):
                s[2] += 1

tot_cmp = sum(v[0] for v in loc.values())
tot_mm = sum(v[1] for v in loc.values())
tot_oh = sum(v[2] for v in loc.values())
print(f"overall: {tot_cmp} comparisons, mismatch {tot_mm/tot_cmp:.4f}, "
      f"opposite-hom share of mismatches {100*tot_oh/tot_mm:.1f}%\n")

# Concentration: how few loci carry how much of the damage?
ranked = sorted(loc.items(), key=lambda kv: -kv[1][2])
n_loci = len(ranked)
print("CONCENTRATION of opposite-homozygote calls")
print(f"  {'top % of loci':>15}{'loci':>9}{'% of all opp-hom':>19}")
cum = 0
targets = [0.01, 0.05, 0.10, 0.25, 0.50]
ti = 0
for i, (_, v) in enumerate(ranked, 1):
    cum += v[2]
    while ti < len(targets) and i >= targets[ti] * n_loci:
        print(f"  {100*targets[ti]:>14.0f}%{i:>9}{100*cum/tot_oh:>18.1f}%")
        ti += 1
loci_with_oh = sum(1 for _, v in ranked if v[2] > 0)
print(f"\n  loci with at least one opposite homozygote: {loci_with_oh} / {n_loci} "
      f"({100*loci_with_oh/n_loci:.1f}%)")

# THE decisive table: what survives if offending loci are removed?
print("\nRESIDUAL after excluding loci with >= k opposite homozygotes")
print(f"  {'k':>4}{'loci kept':>11}{'sites kept':>12}{'comparisons':>13}{'mismatch':>10}{'opp-hom %':>11}")
for k in (1, 2, 3, 5, 10):
    keep = {l for l, v in loc.items() if v[2] < k}
    c = sum(loc[l][0] for l in keep)
    m = sum(loc[l][1] for l in keep)
    o = sum(loc[l][2] for l in keep)
    sites = sum(1 for l, _ in recs if l in keep)
    if c:
        print(f"  {k:>4}{len(keep):>11}{sites:>12}{c:>13}{m/c:>10.4f}"
              f"{100*o/m if m else 0:>10.1f}%")
