#!/usr/bin/env python3
"""Is the replicate-mismatch statistic measuring identity, or just depth pairing?

The libraries fall into two clean groups: 12 "deep" (28-45% of sites called, ~19% heterozygous)
and 18 "shallow" (1-11% called, ~11% heterozygous). Crucially:

  every VALIDATED replicate pair is deep + shallow
  all three "NON-REPLICATE" pairs are deep + deep

So depth combination is perfectly confounded with the verdict, and the 0.69 mismatch attributed
to those three pairs may be nothing more than what a deep-vs-deep comparison yields for ANY two
libraries. That is testable directly: compute mismatch for random pairs within each depth
combination and see whether declared replicates stand out from their own stratum.

  deep-deep random ~= 0.69  -> the three pairs are unremarkable; the "not replicates" verdict
                              is an artefact and must be withdrawn
  deep-deep random ~= 0.39  -> the three pairs really are anomalous, independent of depth

The same stratification also re-frames the headline number: if within-pair mismatch only looks
low because those pairs are deep+shallow, then ~0.27 is a property of the depth combination
rather than a genotyping error rate, and every parameter selected against it is suspect.
"""
import csv, glob, os, re, random, statistics as st, collections

PROJ = "/home/ecastron/araucaria"
random.seed(1)
VCF = sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf"))[0]

samples, rows = [], []
with open(VCF) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            continue
        rows.append([c.split(":")[0].replace("|", "/") for c in f[9:]])

idx = {s: i for i, s in enumerate(samples)}
called = {s: sum(1 for r in rows if r[idx[s]] not in ("./.", ".")) for s in samples}
DEEP_CUT = 0.20 * len(rows)
deep = {s for s in samples if called[s] >= DEEP_CUT}
shallow = {s for s in samples if 0.005 * len(rows) <= called[s] < DEEP_CUT}
print(f"{os.path.basename(VCF)}: {len(rows)} sites | deep={len(deep)} shallow={len(shallow)}\n")

def mm(a, b):
    ia, ib = idx[a], idx[b]
    n = m = 0
    for r in rows:
        x, y = r[ia], r[ib]
        if x in ("./.", ".") or y in ("./.", "."):
            continue
        n += 1
        if x != y:
            m += 1
    return (m / n if n >= 500 else None), n

pairs = [(r["tree_id"], r["run_a"], r["run_b"])
         for r in csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs.tsv"), delimiter="\t")]
known = {frozenset((a, b)) for _, a, b in pairs}

def stratum(x, y):
    d = (x in deep) + (y in deep)
    return {2: "deep-deep", 1: "deep-shallow", 0: "shallow-shallow"}[d]

# Null distribution per depth stratum, excluding declared pairs.
null = collections.defaultdict(list)
for _ in range(4000):
    x, y = random.sample(samples, 2)
    if frozenset((x, y)) in known or x not in deep | shallow or y not in deep | shallow:
        continue
    s = stratum(x, y)
    if len(null[s]) >= 120:
        continue
    e, n = mm(x, y)
    if e is not None:
        null[s].append(e)

print("NULL: random NON-replicate pairs, by depth combination")
for s in ("deep-deep", "deep-shallow", "shallow-shallow"):
    v = null.get(s, [])
    if v:
        print(f"  {s:<16} n={len(v):<4} median {st.median(v):.4f}   "
              f"IQR {st.quantiles(v, n=4)[0]:.4f}-{st.quantiles(v, n=4)[2]:.4f}")

print("\nDECLARED replicate pairs, each against its OWN depth stratum")
print(f"  {'tree':<9}{'stratum':<17}{'mismatch':>9}{'null median':>13}{'verdict':>26}")
for t, a, b in pairs:
    if a not in idx or b not in idx:
        continue
    e, n = mm(a, b)
    if e is None:
        print(f"  {t:<9}{'(insufficient data)':<17}{'':>9}{'':>13}{'untestable':>26}")
        continue
    s = stratum(a, b)
    v = null.get(s, [])
    nm = st.median(v) if v else float("nan")
    if v and e < min(v):
        verdict = "below ALL random - real"
    elif v and e < nm:
        verdict = "below stratum median"
    else:
        verdict = "NOT below stratum median"
    print(f"  {t:<9}{s:<17}{e:>9.4f}{nm:>13.4f}{verdict:>26}")
