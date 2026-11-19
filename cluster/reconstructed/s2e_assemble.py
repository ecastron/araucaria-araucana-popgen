#!/usr/bin/env python3
"""Assemble the corrected Figure S2E table: every intervention with its own same-route baseline.

PRODUCES  results/reconstructed/s2e_mechanisms.tsv   (copied to figdata/s2e_mechanisms.tsv)

INPUTS    results/reconstructed/s2e/binomial_route.tsv   s2e_binomial_route.py (ipyrad cr30 VCF, binomial)
          results/reconstructed/s2e/polyrad_route.tsv    s2e_polyrad_route.R  (34,674 polyRAD loci)
          results/polyrad_fb/fb_F_by_depth.csv           fb_F_by_depth.R (FreeBayes + polyRAD, reference)

DESIGN    Primary bin set for every row: the 7 genotype-depth bins 6-10 ... 100+.
          Route A, binomial all-SNP (VCF): rarefaction, locus mean depth, depth CV, CNV ratio, H_obs and
            structural filters re-expressed on binomial calls, copy-number classes, SD exclusion.
            Baseline = the unfiltered row of the same script (bins with >= 500 genotypes).
          Route B, polyRAD posterior (34,674 loci, "valid" estimator = the one that gives the drawn
            1.207): H_obs ceiling and structural filter as published. Baseline 1.207 (bins > 2000).
          Model row: polyRAD posterior vs binomial calls on the SAME 34,674 loci and read counts.
          Pipeline row: FreeBayes + polyRAD on the reference is a replication with a different locus
            universe; its comparator is the de novo polyRAD route on the same estimator.
          route_baseline_universe says what the baseline is for each row.
Written 19 Nov 2026 (S2E repair).
"""
import csv

P = "/home/ecastron/araucaria"
S = f"{P}/results/reconstructed/s2e"
OUT = f"{P}/results/reconstructed/s2e_mechanisms.tsv"
NM = ['6-10', '10-15', '15-20', '20-30', '30-50', '50-100', '100+']
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
bn = {r["mechanism"]: r for r in rd(f"{S}/binomial_route.tsv")}
pr = {r["mechanism"]: r for r in rd(f"{S}/polyrad_route.tsv") if r["estimator"] == "valid"}
fb = {r["depth"]: r for r in csv.DictReader(open(f"{P}/results/polyrad_fb/fb_F_by_depth.csv"))}

A = "binomial all-SNP (ipyrad cr30 VCF)"
B = "polyRAD posterior d=5 (34,674 de novo loci)"
rows = []


def add(mech, s2e_label, route, src, n_kept, unit, base, filt, nbins, univ, published):
    rows.append(dict(mechanism=mech, s2e_row_published=s2e_label, published_value=published, route=route,
                     n_kept=n_kept, unit=unit, n_bins=nbins, baseline_spread=round(float(base), 3),
                     filtered_spread=round(float(filt), 3),
                     reduction=round(float(base) - float(filt), 3),
                     pct_reduction=round(100 * (float(base) - float(filt)) / float(base), 1),
                     route_baseline_universe=univ, source=src,
                     **{f"F_{b}": src_row(src, mech).get(f"F_{b}", "") for b in NM}))


def src_row(src, mech):
    return {"A": bn, "B": pr}.get(src, {}).get(mech, {})


def fromA(mech, label, pub):
    r = bn[mech]
    add(mech, label, A, "A", r["n_snps_kept"], "SNPs", r["baseline_spread"], r["filtered_spread"], r["n_bins"],
        "same SNPs, unfiltered, 7 bins", pub)


def fromB(mech, label, pub):
    r = pr[mech]
    add(mech, label, B, "B", r["n_loci_kept"], "loci", r["baseline_spread"], r["filtered_spread"], r["n_bins"],
        "same 34,674 loci, unfiltered polyRAD, 7 bins", pub)


fromA("none - unfiltered baseline", "none - unfiltered baseline", 1.207)
fromA("rarefy every genotype to 30 reads", "rarefy every genotype to 30 reads", 1.093)
fromA("locus mean depth <= 1.25x median", "locus mean depth <= 1.25x median (15% dropped)", 1.138)
fromA("depth-CV filter (CV <= 1.0)", "depth-CV filter (166 of 76,823 loci)", 1.207)
fromA("locus-relative depth (CNV) filter <= 1.5", "locus-relative depth (CNV) filter <= 1.5", 1.337)
fromA("H_obs ceiling 0.20 (binomial calls)", "", "")
fromA("structural filter cv<=0.80 depth<=95 miss<=0.70", "", "")
for m in ("single-copy loci only (genome mapping)", "exclude segmental-duplication loci"):
    fromA(m, "", "")
fromB("none - unfiltered baseline (polyRAD posterior)", "none - unfiltered baseline", 1.207)
fromB("H_obs ceiling 0.20", "H_obs ceiling 0.20 (41% of loci dropped)", 1.146)
fromB("structural filter cv<=0.80 depth<=95 miss<=0.70", "structural locus filter, moderate setting", 1.538)
# genotype model: polyRAD posterior against binomial calls on the same loci and reads
bm = pr["binomial calls, same loci (genotype-model baseline)"]; pb = pr["none - unfiltered baseline (polyRAD posterior)"]
rows.append(dict(mechanism="beta-binomial overdispersion, polyRAD d = 5 (vs binomial calls, same loci)",
                 s2e_row_published="beta-binomial overdispersion, polyRAD d = 5", published_value=1.175,
                 route="binomial calls -> polyRAD posterior, same 34,674 loci", n_kept=pb["n_loci_kept"], unit="loci",
                 n_bins=pb["n_bins"], baseline_spread=round(float(bm["filtered_spread"]), 3),
                 filtered_spread=round(float(pb["filtered_spread"]), 3),
                 reduction=round(float(bm["filtered_spread"]) - float(pb["filtered_spread"]), 3),
                 pct_reduction=round(100 * (float(bm["filtered_spread"]) - float(pb["filtered_spread"]))
                                     / float(bm["filtered_spread"]), 1),
                 route_baseline_universe="binomial calls from the same read counts, same loci, 7 bins",
                 source="B", **{f"F_{b}": pb[f"F_{b}"] for b in NM}))
# pipeline replication: FreeBayes + polyRAD on the reference, same estimator as route B
fF = [float(fb[b]["F_IS"]) for b in NM if b in fb]
fbs = max(fF) - min(fF)
rows.append(dict(mechanism="independent pipeline: bowtie2 + FreeBayes + polyRAD d = 7 (reference)",
                 s2e_row_published="independent pipeline: bowtie2 + FreeBayes", published_value=1.231,
                 route="FreeBayes + polyRAD (reference); comparator = de novo polyRAD, same estimator",
                 n_kept=122051, unit="loci (different universe)", n_bins=len(fF),
                 baseline_spread=round(float(pb["filtered_spread"]), 3), filtered_spread=round(fbs, 3),
                 reduction=round(float(pb["filtered_spread"]) - fbs, 3),
                 pct_reduction=round(100 * (float(pb["filtered_spread"]) - fbs) / float(pb["filtered_spread"]), 1),
                 route_baseline_universe="replication, not a filter: different loci; de novo polyRAD shown as comparator",
                 source="fb_F_by_depth.csv", **{f"F_{b}": fb[b]["F_IS"] if b in fb else "" for b in NM}))

cols = list(rows[0].keys())
with open(OUT, "w") as fo:
    w = csv.DictWriter(fo, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
print("%-74s %-6s %7s %7s %7s %6s" % ("mechanism", "src", "base", "filt", "reduct", "pct"))
for r in rows:
    print("%-74s %-6s %7.3f %7.3f %7.3f %6.1f" % (r["mechanism"][:74], r["source"][:6], r["baseline_spread"],
                                                 r["filtered_spread"], r["reduction"], r["pct_reduction"]))
print(f"written -> {OUT}")
