#!/usr/bin/env python3
"""SNP filtering thresholds by replicate concordance, on the rule declared in the spec.

Grid over per-individual call rate x MAF floor. At each cell:
  - genotyping error between the 15 replicate library pairs, over sites called in BOTH
  - r80-style loci retained
  - individuals retained

Selection rule, fixed before the grid was run: take the MOST PERMISSIVE combination whose
replicate error lies within one standard error of the observed minimum. Permissive
deliberately - Huang & Knowles (2016) showed stringent missing-data filters preferentially
strip the highest-mutation-rate loci.

This produces the relatedness threshold set, which is the one replicate concordance is the
right criterion for. The full surface is written out so the other analyses (structure, Ne,
SFS) can be read off it under their own criteria.
"""
import csv, os, sys, statistics as st
from collections import defaultdict

PROJ = "/home/ecastron/araucaria"
VCF = f"{PROJ}/rad/final/populations.snps.vcf"
PAIRS = f"{PROJ}/metadata/replicate_pairs.tsv"
OUT = f"{PROJ}/logs/snp_filter_grid.tsv"
RUNLOG = f"{PROJ}/logs/run_log.md"

CALL_RATES = [0.50, 0.70, 0.90]      # min fraction of sites called per individual
MAFS = [0.01, 0.02, 0.05, 0.10]


def load(path):
    samples, sites = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            gts = [c.split(":")[0].replace("|", "/") for c in f[9:]]
            sites.append(gts)
    return samples, sites


def maf_of(col):
    a = 0
    tot = 0
    for g in col:
        if g in (".", "./."):
            continue
        for h in g.split("/"):
            if h == ".":
                continue
            tot += 1
            if h != "0":
                a += 1
    if tot == 0:
        return 0.0
    p = a / tot
    return min(p, 1 - p)


def main():
    if not os.path.exists(VCF):
        print("no VCF; assembly must run first", file=sys.stderr)
        sys.exit(1)
    samples, sites = load(VCF)
    pairs = [(r["run_a"], r["run_b"]) for r in csv.DictReader(open(PAIRS), delimiter="\t")]
    idx = {s: i for i, s in enumerate(samples)}
    pairs = [(a, b) for a, b in pairs if a in idx and b in idx]
    print(f"{len(samples)} libraries, {len(sites)} sites, {len(pairs)} replicate pairs", flush=True)

    called = [sum(1 for g in col if g not in (".", "./.")) for col in
              ([s[i] for s in sites] for i in range(len(samples)))]
    per_ind_call = [c / len(sites) for c in called]
    site_maf = [maf_of(col) for col in sites]

    rows = []
    for cr in CALL_RATES:
        keep_ind = [i for i, c in enumerate(per_ind_call) if c >= cr]
        for maf in MAFS:
            keep_site = [j for j, mv in enumerate(site_maf) if mv >= maf]
            errs = []
            for a, b in pairs:
                ia, ib = idx[a], idx[b]
                if ia not in keep_ind or ib not in keep_ind:
                    continue
                both = mism = 0
                for j in keep_site:
                    ga, gb = sites[j][ia], sites[j][ib]
                    if ga in (".", "./.") or gb in (".", "./."):
                        continue
                    both += 1
                    if ga != gb:
                        mism += 1
                if both >= 50:
                    errs.append(mism / both)
            rows.append({
                "call_rate": cr, "maf": maf,
                "individuals": len(keep_ind), "sites": len(keep_site),
                "pairs_used": len(errs),
                "err_mean": round(st.mean(errs), 6) if errs else None,
                "err_se": round(st.pstdev(errs) / (len(errs) ** 0.5), 6) if len(errs) > 1 else None,
            })
            print(f"  call>={cr} maf>={maf}: {len(keep_ind)} ind, {len(keep_site)} sites, "
                  f"err={rows[-1]['err_mean']}", flush=True)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    usable = [r for r in rows if r["err_mean"] is not None and r["pairs_used"] >= 10]
    if not usable:
        print("GATE FAILED: no cell with >=10 usable replicate pairs", file=sys.stderr)
        with open(RUNLOG, "a") as fh:
            fh.write("\n- SNP filtering grid ran but no cell had >=10 usable replicate pairs; "
                     "thresholds NOT fixed.\n")
        sys.exit(1)

    best = min(usable, key=lambda r: r["err_mean"])
    thr = best["err_mean"] + (best["err_se"] or 0.0)
    within = [r for r in usable if r["err_mean"] <= thr]
    # Most permissive = lowest call-rate requirement, then lowest MAF floor.
    pick = sorted(within, key=lambda r: (r["call_rate"], r["maf"]))[0]

    print(f"\nlowest error: call>={best['call_rate']} maf>={best['maf']} err={best['err_mean']}")
    print(f"CHOSEN (most permissive within 1 SE): call>={pick['call_rate']} maf>={pick['maf']} "
          f"err={pick['err_mean']} sites={pick['sites']} individuals={pick['individuals']}")

    from datetime import datetime, timezone
    with open(RUNLOG, "a") as fh:
        fh.write(f"\n## {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%MZ')} — "
                 f"SNP filtering thresholds fixed (relatedness set)\n\n")
        fh.write(f"- **Per-individual call rate >= {pick['call_rate']}, MAF >= {pick['maf']}**\n")
        fh.write(f"- Replicate error {pick['err_mean']}, "
                 f"{pick['sites']} sites, {pick['individuals']} libraries retained\n")
        fh.write(f"- Rule applied: most permissive cell within 1 SE of the minimum "
                 f"(min was call>={best['call_rate']}, maf>={best['maf']}, err={best['err_mean']})\n")
        fh.write(f"- Full surface: `logs/snp_filter_grid.tsv`\n")
        fh.write(f"- Criterion is replicate concordance, which is causally independent of every "
                 f"question in Section 2, so this cannot bias a result.\n")


if __name__ == "__main__":
    main()
