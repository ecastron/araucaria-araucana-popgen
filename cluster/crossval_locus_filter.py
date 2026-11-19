#!/usr/bin/env python3
"""Cross-validate the opposite-homozygote locus filter.

Filtering loci on replicate disagreement and then measuring replicate disagreement on the SAME
pairs is selection on the dependent variable: mismatch falls from 0.484 to 0.149 by
construction, whether or not the filter identifies anything real.

The honest test splits the pairs. Loci are blacklisted using TRAIN pairs only, then mismatch is
measured on HELD-OUT pairs the filter never saw.

  held-out mismatch stays low  -> the bad loci are a stable property of the assembly (genuinely
                                  non-homologous clusters). Filtering gives a usable dataset.
  held-out mismatch rebounds   -> the filter memorised which sites happened to disagree in the
                                  training pairs. It generalises to nothing and the apparent
                                  improvement is an artefact.

Reported against two baselines on the same held-out pairs: no filter at all, and a filter that
removes an equal NUMBER of randomly chosen loci. The random-filter baseline matters because
discarding half the loci changes the site composition on its own, and any honest gain has to
beat that, not just beat zero.
"""
import csv, glob, os, re, random, statistics as st, collections

PROJ = "/home/ecastron/araucaria"
random.seed(7)
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
        recs.append((f[0], f[9:]))
idx = {s: i for i, s in enumerate(samples)}
pairs = [p for p in pairs if p[1] in idx and p[2] in idx]
all_loci = sorted({l for l, _ in recs})
print(f"{os.path.basename(VCF)}: {len(recs)} sites, {len(all_loci)} loci, {len(pairs)} pairs\n")

def gt(cell):
    return cell.split(":")[0].replace("|", "/")

def het(g):
    a, _, b = g.partition("/")
    return a != b

def blacklist(train):
    bad = set()
    for locus, cells in recs:
        for _, a, b in train:
            ga, gb = gt(cells[idx[a]]), gt(cells[idx[b]])
            if ga in ("./.", ".") or gb in ("./.", "."):
                continue
            if ga != gb and not het(ga) and not het(gb):
                bad.add(locus)
                break
    return bad

def mismatch(test, drop):
    n = m = 0
    for locus, cells in recs:
        if locus in drop:
            continue
        for _, a, b in test:
            ga, gb = gt(cells[idx[a]]), gt(cells[idx[b]])
            if ga in ("./.", ".") or gb in ("./.", "."):
                continue
            n += 1
            if ga != gb:
                m += 1
    return (m / n if n else float("nan")), n

print(f"  {'fold':>5}{'train':>7}{'test':>6}{'blacklisted':>13}{'no filter':>11}"
      f"{'CV filter':>11}{'random filter':>15}")
res = []
order = pairs[:]
random.shuffle(order)
half = len(order) // 2
for fold, (tr, te) in enumerate([(order[:half], order[half:]),
                                 (order[half:], order[:half])], 1):
    bad = blacklist(tr)
    base, n0 = mismatch(te, set())
    filt, n1 = mismatch(te, bad)
    rnd = set(random.sample(all_loci, len(bad)))
    rndm, _ = mismatch(te, rnd)
    res.append((base, filt, rndm))
    print(f"  {fold:>5}{len(tr):>7}{len(te):>6}{len(bad):>13}{base:>11.4f}"
          f"{filt:>11.4f}{rndm:>15.4f}")

b = st.mean(r[0] for r in res); f = st.mean(r[1] for r in res); r_ = st.mean(r[2] for r in res)
print(f"\n  mean over folds: no filter {b:.4f} | CV filter {f:.4f} | random filter {r_:.4f}")
print(f"  gain attributable to the filter (vs random of equal size): {r_ - f:+.4f}")
print(f"  in-sample (circular) figure for comparison: 0.1491")
if f > r_ - 0.02:
    print("\n  >>> The filter does not beat dropping the same number of RANDOM loci.")
    print("  >>> It generalises to nothing; the in-sample improvement was an artefact.")
elif f < 0.10:
    print("\n  >>> Held-out mismatch is low: the bad loci are a stable property of the")
    print("  >>> assembly and filtering yields a usable dataset.")
else:
    print("\n  >>> Partial: the filter beats random but held-out mismatch remains too high")
    print("  >>> for reliable individual genotypes. Usable for relative, not absolute, work.")
