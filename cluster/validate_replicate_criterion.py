#!/usr/bin/env python3
"""Are the replicate pairs actually replicates?

Every parameter decision in the RAD half rests on one number: mismatch between the two
libraries of a tree. Both assemblers now put that at ~0.28, and the conclusion drawn from it --
that the data cannot support individual genotype calling -- is only valid if the criterion
itself is sound. It is not sound if the "replicate" pairs are not really the same tree.

The test is a null comparison the sweep never made: compare WITHIN-PAIR mismatch against
mismatch between RANDOM library pairs from the same assembly.

  within << random : the pairs are genuinely the same individual and the residual mismatch is
                     genotyping noise. The verdict stands -- use genotype likelihoods.
  within ~= random : the pairs are no more similar than unrelated trees. The pairs are not
                     replicates, the criterion is meaningless, and every sweep result selected
                     on it must be discarded rather than reinterpreted.

Also reported per pair, because a criterion can be sound in aggregate while individual pairs
are mislabelled -- and a single bad pair dragging the median would be invisible otherwise.
"""
import csv, glob, os, re, random, statistics as st

PROJ = "/home/ecastron/araucaria"
random.seed(1)

def load(vcf):
    samples, cols = [], []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
                continue
            cols.append([c.split(":")[0].replace("|", "/") for c in f[9:]])
    return samples, cols

def mismatch(cols, ia, ib):
    n = m = 0
    for r in cols:
        a, b = r[ia], r[ib]
        if a in ("./.", ".") or b in ("./.", "."):
            continue
        n += 1
        if a != b:
            m += 1
    return (m / n if n else None), n

pairs = [(r["run_a"], r["run_b"], r["tree_id"])
         for r in csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs.tsv"), delimiter="\t")]

for vcf in sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf")):
    ct = os.path.basename(vcf).replace(".vcf", "")
    samples, cols = load(vcf)
    idx = {s: i for i, s in enumerate(samples)}
    print(f"\n{'='*74}\n{ct}  ({len(samples)} libraries, {len(cols)} sites)\n{'='*74}")

    within = []
    print(f"  {'tree':<10}{'lib A':<14}{'lib B':<14}{'shared':>9}{'mismatch':>10}")
    for a, b, t in pairs:
        if a not in idx or b not in idx:
            continue
        e, n = mismatch(cols, idx[a], idx[b])
        if e is None:
            continue
        within.append(e)
        print(f"  {t:<10}{a:<14}{b:<14}{n:>9}{e:>10.4f}")

    # Null: random library pairs, excluding the declared replicate pairs themselves.
    known = {frozenset((a, b)) for a, b, _ in pairs}
    rnd = []
    tries = 0
    while len(rnd) < 200 and tries < 5000:
        tries += 1
        x, y = random.sample(samples, 2)
        if frozenset((x, y)) in known:
            continue
        e, n = mismatch(cols, idx[x], idx[y])
        if e is not None and n >= 500:
            rnd.append(e)

    if within and rnd:
        w, r = st.median(within), st.median(rnd)
        print(f"\n  within-pair  median mismatch : {w:.4f}  (n = {len(within)} pairs)")
        print(f"  RANDOM pair  median mismatch : {r:.4f}  (n = {len(rnd)} draws, "
              f"IQR {st.quantiles(rnd, n=4)[0]:.4f}-{st.quantiles(rnd, n=4)[2]:.4f})")
        print(f"  ratio within/random          : {w / r:.3f}")
        below = sum(1 for e in within if e < min(rnd))
        print(f"  pairs below the MINIMUM random mismatch: {below} / {len(within)}")
        if w / r > 0.9:
            print("\n  >>> The declared replicates are NO MORE SIMILAR than unrelated libraries.")
            print("  >>> The criterion is invalid. Sweep selections made on it cannot be used.")
        elif w / r < 0.5:
            print("\n  >>> Replicates are clearly more similar than unrelated libraries: the")
            print("  >>> pairing is real and the residual mismatch is genotyping error.")
        else:
            print("\n  >>> Intermediate: replicates are somewhat more similar, but far less so")
            print("  >>> than true replicates should be. Suspect a subset of mislabelled pairs.")
