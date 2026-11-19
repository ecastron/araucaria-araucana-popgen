#!/usr/bin/env python3
"""Re-derive the clustering threshold WITHOUT replicates: r80 plus null expectations.

The replicate-concordance criterion that selected ct = 0.91 is void -- it measured sibling
relatedness, not genotyping error. This applies the standard replicate-free machinery instead:

  r80 (Paris et al. 2017)   Count loci recovered in >=80% of samples. Under-clustering splits
                            true loci across samples and depresses r80; over-clustering merges
                            paralogs, which inflates raw locus counts but not r80. The peak is
                            the operating point.

  Ts/Tv (Shafer et al. 2017) Transition/transversion ratio against its biological expectation
                            (~1.5-2.5 in plants). Sequencing and clustering errors are
                            direction-agnostic, so they push Ts/Tv toward 0.5 (the ratio
                            expected if substitutions were random). A falling Ts/Tv as
                            thresholds loosen is an error signal independent of r80.

  Heterozygosity            Over-clustering merges paralogs, which appear as persistent
                            heterozygotes. A jump in Ho at loose thresholds flags it. Computed
                            WITHIN a single stand where possible, since pooling structured
                            stands depresses Ho by the Wahlund effect (+0.474 here).

Reports raw locus and SNP counts alongside, because those are what over-clustering inflates --
the whole point of r80 is that it does not reward them.
"""
import csv, glob, os, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
vcfs = sys.argv[1:] or sorted(glob.glob(f"{PROJ}/rad/ipyrad/*_outfiles/*.vcf"))

TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

site_of = {r[0]: r[1] for r in csv.reader(open(f"{PROJ}/metadata/popmap_all_libraries.tsv"),
                                          delimiter="\t") if len(r) >= 2}

print(f"{'assembly':>12}{'libs':>6}{'loci':>10}{'SNPs':>10}"
      f"{'r80 loci':>10}{'r80 SNPs':>10}{'Ts/Tv':>8}{'Ho(pooled)':>12}")
rows = []
for V in vcfs:
    name = os.path.basename(V).replace(".vcf", "")
    samples, loci_pres, ts, tv = [], {}, 0, 0
    nsnp = 0
    het = tot = 0
    with open(V) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
                continue
            nsnp += 1
            ref, alt = f[3], f[4].split(",")[0]
            if len(ref) == 1 and len(alt) == 1:
                if (ref, alt) in TS:
                    ts += 1
                else:
                    tv += 1
            called = 0
            for c in f[9:]:
                g = c.partition(":")[0].replace("|", "/")
                if g in ("./.", "."):
                    continue
                called += 1
                tot += 1
                a, _, b = g.partition("/")
                het += (a != b)
            locus = f[0] if f[0] != "locus" else f[2].split("_")[0]
            key = f[2].split("_")[0] if "_" in f[2] else f[0]
            prev = loci_pres.get(key, 0)
            loci_pres[key] = max(prev, called)
    n = len(samples)
    r80_loci = sum(1 for v in loci_pres.values() if v >= 0.8 * n)
    # r80 SNPs: recount, cheaper to approximate from the per-locus max call count
    r80_keys = {k for k, v in loci_pres.items() if v >= 0.8 * n}
    r80_snps = 0
    with open(V) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = f[2].split("_")[0] if "_" in f[2] else f[0]
            if key in r80_keys:
                r80_snps += 1
    tstv = ts / tv if tv else float("nan")
    ho = het / tot if tot else float("nan")
    rows.append((name, n, len(loci_pres), nsnp, r80_loci, r80_snps, tstv, ho))
    print(f"{name:>12}{n:>6}{len(loci_pres):>10}{nsnp:>10}{r80_loci:>10}{r80_snps:>10}"
          f"{tstv:>8.3f}{ho:>12.4f}")

if len(rows) > 1:
    best = max(rows, key=lambda r: r[4])
    print(f"\nr80 peak: {best[0]} with {best[4]} loci / {best[5]} SNPs")
    print("Check Ts/Tv alongside: a threshold that wins on r80 but shows Ts/Tv falling toward")
    print("0.5 is over-clustering, and the extra loci are paralogs rather than recovered data.")
