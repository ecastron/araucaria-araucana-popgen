#!/usr/bin/env python3
"""Summarise the M sweep: r80 loci recovered, and genotyping error from replicate pairs.

These are the two criteria the spec fixes for choosing M (and n = M):
  - r80: loci present in at least 80% of samples (Paris et al. 2017)
  - genotyping error measured on the 15 trees sequenced as two independent libraries
    (Mastretta-Yanes et al. 2015), which is the primary criterion

Error is computed only over sites called in BOTH libraries of a pair, so it measures
disagreement rather than missingness, and missingness is reported separately. Both matter:
a setting that minimises error by calling almost nothing is not a good setting.

Output is a table, not a decision. The choice is recorded separately once reviewed.
"""
import csv, gzip, os, sys, statistics as st
from collections import defaultdict

PROJ = "/home/ecastron/araucaria"
SWEEP = f"{PROJ}/rad/sweep"
PAIRS = f"{PROJ}/metadata/replicate_pairs.tsv"
OUT = f"{PROJ}/logs/sweep_summary.tsv"


def read_vcf(path):
    """Return {sample: {site: genotype}} using only unphased diploid calls."""
    gts = defaultdict(dict)
    samples = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            site = (f[0], f[1])
            for s, cell in zip(samples, f[9:]):
                g = cell.split(":")[0]
                if g not in (".", "./.", ".|."):
                    gts[s][site] = g.replace("|", "/")
    return gts, samples


def main():
    pairs = [r for r in csv.DictReader(open(PAIRS), delimiter="\t")]
    rows = []
    for M in range(1, 9):
        vcf = f"{SWEEP}/M{M}/populations.snps.vcf"
        if not os.path.exists(vcf) or os.path.getsize(vcf) == 0:
            print(f"M{M}: no VCF", file=sys.stderr)
            continue
        gts, samples = read_vcf(vcf)
        n_sites = len({s for g in gts.values() for s in g})

        errs, shared, miss = [], [], []
        for p in pairs:
            a, b = p["run_a"], p["run_b"]
            if a not in gts or b not in gts:
                continue
            both = set(gts[a]) & set(gts[b])
            either = set(gts[a]) | set(gts[b])
            if len(both) < 50:
                continue
            mism = sum(1 for s in both if gts[a][s] != gts[b][s])
            errs.append(mism / len(both))
            shared.append(len(both))
            miss.append(1 - len(both) / len(either) if either else float("nan"))

        rows.append({
            "M": M, "n": M,
            "r80_snps": n_sites,
            "samples": len(samples),
            "pairs_used": len(errs),
            "err_mean": round(st.mean(errs), 5) if errs else None,
            "err_sd": round(st.pstdev(errs), 5) if len(errs) > 1 else None,
            "err_se": round(st.pstdev(errs) / (len(errs) ** 0.5), 5) if len(errs) > 1 else None,
            "shared_sites_median": int(st.median(shared)) if shared else None,
            "pair_missing_frac": round(st.mean(miss), 4) if miss else None,
        })
        print(f"M{M}: r80={n_sites} err={rows[-1]['err_mean']} pairs={len(errs)}", flush=True)

    if not rows:
        print("no results", file=sys.stderr)
        return
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {OUT}\n")
    print(f"{'M':>2} {'r80_snps':>9} {'err_mean':>9} {'err_se':>8} {'shared':>8} {'pair_miss':>10}")
    for r in rows:
        print(f"{r['M']:>2} {r['r80_snps']:>9} {str(r['err_mean']):>9} {str(r['err_se']):>8} "
              f"{str(r['shared_sites_median']):>8} {str(r['pair_missing_frac']):>10}")

    # The spec's rule for the filtering grid is "most permissive within 1 SE of the minimum".
    # Report the analogous candidate here for review; do not treat it as the decision.
    good = [r for r in rows if r["err_mean"] is not None]
    if good:
        best = min(good, key=lambda r: r["err_mean"])
        thr = best["err_mean"] + (best["err_se"] or 0)
        within = [r for r in good if r["err_mean"] <= thr]
        pick = max(within, key=lambda r: r["r80_snps"])
        print(f"\nlowest error: M={best['M']} ({best['err_mean']})")
        print(f"within 1 SE and maximising r80: M={pick['M']} "
              f"(err={pick['err_mean']}, r80={pick['r80_snps']})")
        print("\nFor review - not applied.")


if __name__ == "__main__":
    main()
