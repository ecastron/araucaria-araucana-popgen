#!/usr/bin/env python3
"""Replicate concordance across the ipyrad clustering-threshold sweep.

Same criterion used for the Stacks sweep, so the two are directly comparable: mismatch between
the two libraries of a tree, over sites called in both. Stacks gave 0.27-0.33 under every
treatment tried. If ipyrad lands near 0.02-0.05 the problem was Stacks' catalog matching; if it
lands near 0.2-0.3 the data cannot support individual genotype calling.
"""
import re
import csv, glob, os, statistics as st, collections

PROJ = "/home/ecastron/araucaria"
WORK = f"{PROJ}/rad/ipyrad"
PAIRS = f"{PROJ}/metadata/replicate_pairs.tsv"

def load(vcf):
    samples, cols = [], []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("##"): continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]; continue
            cols.append([c.split(":")[0].replace("|", "/") for c in f[9:]])
    # ipyrad resolves each symlink to its target and names the sample from THAT filename, so
    # samples arrive as SRR11841624.1 (the mate suffix of the cleaned file) while every other
    # table in the project keys on the bare accession. Normalise, or every replicate pair
    # silently fails to match and the sweep reports zero pairs instead of an error.
    samples = [re.sub(r"\.[12]$", "", s) for s in samples]
    return samples, cols

pairs = [(r["run_a"], r["run_b"], r["tree_id"]) for r in
         csv.DictReader(open(PAIRS), delimiter="\t")]

print(f"{'ct':>6} {'snps':>9} {'libs':>5} {'pairs':>6} {'median_err':>11} {'shared_med':>11}")
rows = []
for vcf in sorted(glob.glob(f"{WORK}/*_outfiles/*.vcf")):
    ct = os.path.basename(vcf).replace(".vcf", "")
    samples, cols = load(vcf)
    idx = {s: i for i, s in enumerate(samples)}
    errs, shared = [], []
    for a, b, t in pairs:
        if a not in idx or b not in idx: continue
        ia, ib = idx[a], idx[b]
        n = m = 0
        for r in cols:
            ga, gb = r[ia], r[ib]
            if ga in (".", "./.") or gb in (".", "./."): continue
            n += 1
            if ga != gb: m += 1
        if n >= 200:
            errs.append(m / n); shared.append(n)
    rows.append((ct, len(cols), len(samples), len(errs),
                 st.median(errs) if errs else None,
                 int(st.median(shared)) if shared else None))
    print(f"{ct:>6} {len(cols):>9} {len(samples):>5} {len(errs):>6} "
          f"{('%.4f' % rows[-1][4]) if rows[-1][4] else 'NA':>11} {str(rows[-1][5]):>11}")

print("\nStacks reference: median replicate mismatch 0.27-0.33 across M=1-8,")
print("unchanged by depth filtering, HDplot paralog filtering, or PCR-duplicate removal.")
best = [r for r in rows if r[4] is not None]
if best:
    b = min(best, key=lambda r: r[4])
    print(f"\nbest ipyrad: {b[0]} at {b[4]:.4f}")
    if b[4] < 0.10:
        print("VERDICT: markedly better than Stacks -> catalog matching was the problem.")
    else:
        print("VERDICT: comparable to Stacks -> the data limit individual genotype calling;")
        print("         Q6 needs genotype likelihoods / depth-aware relatedness.")
