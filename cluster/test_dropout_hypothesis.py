#!/usr/bin/env python3
"""Does allelic dropout explain the replicate mismatch, or does something else?

HYPOTHESIS: at low depth a true heterozygote often shows only one allele and is called
homozygous. That is a systematic bias, not noise, and it is what genotype likelihoods correct.

The ipyrad VCF carries GT:DP:CATG per sample, so depth is available per genotype and the test
can be run WITHIN a replicate pair -- controlling for individual, library prep and locus set --
instead of correlating across libraries where everything varies at once.

THREE PREDICTIONS, all stated before running:

  P1  Mismatch falls as the minimum depth required of BOTH members rises.
      If it plateaus high, dropout is not the explanation.

  P2  Mismatch at low depth is dominated by het-vs-hom disagreements. Dropout can turn a het
      into a hom, so it PRODUCES het/hom conflicts. It CANNOT produce opposite homozygotes
      (0/0 vs 1/1) -- that needs a different mechanism (paralogy, contamination, wrong sample).
      So the opposite-homozygote share should be small and should shrink with depth.

  P3  Within a library, heterozygosity rises with depth, because dropout suppresses het calls.

FAILURE MEANS: if mismatch stays high at high depth, or opposite homozygotes stay common, then
the problem is not dropout, genotype likelihoods will not fix it, and Q6 needs rescoping rather
than a better estimator.
"""
import csv, glob, os, re, collections

PROJ = "/home/ecastron/araucaria"
VCF = sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf"))[0]

VALID = f"{PROJ}/metadata/replicate_pairs_validated.tsv"
pairs = [(r["tree_id"], r["run_a"], r["run_b"])
         for r in csv.DictReader(open(VALID), delimiter="\t")
         if r["status"] == "validated"]

samples, recs = [], []
with open(VCF) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            fmt_i = 8
            continue
        recs.append(f[9:])
idx = {s: i for i, s in enumerate(samples)}
print(f"{os.path.basename(VCF)}: {len(samples)} libraries, {len(recs)} sites, "
      f"{len(pairs)} validated pairs\n")

def parse(cell):
    p = cell.split(":")
    gt = p[0].replace("|", "/")
    try:
        dp = int(p[1])
    except (IndexError, ValueError):
        dp = 0
    return gt, dp

def is_het(g):
    a, _, b = g.partition("/")
    return a != b

# ---- P1 + P2: mismatch and its composition, as a function of minimum depth ---------------
print("P1/P2  mismatch vs minimum depth required of BOTH members (validated pairs pooled)")
print(f"  {'min DP':>7}{'sites':>10}{'mismatch':>10}{'het-vs-hom':>12}{'opposite hom':>14}")
for cut in (1, 2, 3, 5, 8, 10, 15, 20, 30):
    n = m = hh = oh = 0
    for _, a, b in pairs:
        if a not in idx or b not in idx:
            continue
        ia, ib = idx[a], idx[b]
        for r in recs:
            ga, da = parse(r[ia])
            gb, db = parse(r[ib])
            if ga in ("./.", ".") or gb in ("./.", ".") or da < cut or db < cut:
                continue
            n += 1
            if ga == gb:
                continue
            m += 1
            if is_het(ga) != is_het(gb):
                hh += 1
            elif not is_het(ga) and not is_het(gb):
                oh += 1
    if n:
        print(f"  {cut:>7}{n:>10}{m/n:>10.4f}{100*hh/m if m else 0:>11.1f}%"
              f"{100*oh/m if m else 0:>13.1f}%")

# ---- P3: heterozygosity vs depth, within library ----------------------------------------
print("\nP3  heterozygosity by depth bin, within library (all libraries pooled)")
bins = [(1, 2), (3, 4), (5, 7), (8, 11), (12, 19), (20, 10**9)]
tally = collections.defaultdict(lambda: [0, 0])
for r in recs:
    for s in samples:
        g, d = parse(r[idx[s]])
        if g in ("./.", ".") or d < 1:
            continue
        for lo, hi in bins:
            if lo <= d <= hi:
                tally[(lo, hi)][0] += 1
                tally[(lo, hi)][1] += is_het(g)
                break
print(f"  {'depth':>10}{'genotypes':>12}{'het rate':>10}")
for lo, hi in bins:
    tot, het = tally[(lo, hi)]
    if tot:
        lab = f"{lo}-{hi}" if hi < 10**9 else f"{lo}+"
        print(f"  {lab:>10}{tot:>12}{het/tot:>10.4f}")
