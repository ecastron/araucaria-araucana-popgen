#!/usr/bin/env python3
"""Why is observed heterozygosity so low? Four candidates, separated by their signatures.

Observed Ho is ~0.016 overall, rising with depth to ~0.095. That is far below expectation for a
wind-pollinated outcrossing conifer. Four mechanisms push Ho the same direction, and they have
distinguishable signatures. Deciding between them determines whether the data are damaged or
simply structured.

  WAHLUND        Pooling differentiated stands inflates apparent homozygosity. We have already
                 demonstrated real structure (isolation by distance, Mantel r = +0.601).
                 SIGNATURE: FIS near zero WITHIN stands, clearly positive when pooled.
                 -> biological, expected, no correction needed for within-stand work.

  INBREEDING     SIGNATURE: FIS positive both within stands AND pooled, fairly uniform.
                 -> biological, plausible in a fragmented relict, must be modelled not removed.

  NULL ALLELES   A mutated restriction site removes an allele from the library entirely.
                 SIGNATURE: DEPTH-INDEPENDENT deficit, concentrated at particular loci, and
                 correlated with missing data (the homozygous-null individuals fail to amplify).
                 -> a RAD artefact; affected loci should be filtered.

  CALLING ERROR  SIGNATURE: strongly DEPTH-DEPENDENT deficit, spread evenly across loci.
                 -> technical; mitigated by depth filtering or genotype likelihoods.

The FIS within-versus-pooled contrast is the decisive comparison, because Wahlund is the only
mechanism that vanishes when computed within stands.
"""
import csv, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1] if len(sys.argv) > 1 else (
    f"{PROJ}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")
MINDP = 20      # high-depth genotypes only, so calling error is largely excluded

samples, G, D = [], [], []
for line in open(V):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
        continue
    g = np.full(len(samples), np.nan)
    d = np.zeros(len(samples), dtype=np.int32)
    for k, c in enumerate(f[9:]):
        p = c.split(":")
        gt = p[0].replace("|", "/")
        try:
            d[k] = int(p[1])
        except (IndexError, ValueError):
            d[k] = 0
        if gt in ("./.", "."):
            continue
        a, _, b = gt.partition("/")
        g[k] = (a != "0") + (b != "0")
    G.append(g); D.append(d)
G = np.vstack(G); D = np.vstack(D)
idx = {s: i for i, s in enumerate(samples)}
print(f"{V.split('/')[-1]}: {G.shape[0]} sites, {G.shape[1]} libraries")

site_of = {r[0]: r[1] for r in csv.reader(open(f"{PROJ}/metadata/popmap_177.tsv"), delimiter="\t")
           if len(r) >= 2}
stands = sorted({site_of[s] for s in samples if s in site_of})

def fis(gsub):
    """Weir & Cockerham style FIS = 1 - Ho/He, per locus, averaged over informative loci."""
    n = np.sum(~np.isnan(gsub), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.nansum(gsub, axis=1) / (2 * n)
        ho = np.nansum(gsub == 1, axis=1) / n
        he = 2 * p * (1 - p)
        ok = (n >= 8) & (he > 0.05)
        f = 1 - ho[ok] / he[ok]
    return float(np.mean(f)), float(np.mean(ho[ok])), float(np.mean(he[ok])), int(ok.sum())

hi = np.where(D >= MINDP, G, np.nan)

print(f"\nusing genotypes at >= {MINDP}x only\n")
print("A. FIS WITHIN each stand (Wahlund cannot act here)")
print(f"  {'stand':<8}{'libs':>6}{'loci':>8}{'Ho':>9}{'He':>9}{'FIS':>9}")
within = []
for st in stands:
    cols = [idx[s] for s in samples if site_of.get(s) == st]
    if len(cols) < 6:
        continue
    f_, ho, he, nl = fis(hi[:, cols])
    within.append(f_)
    print(f"  {st:<8}{len(cols):>6}{nl:>8}{ho:>9.4f}{he:>9.4f}{f_:>9.4f}")

allcols = [idx[s] for s in samples if s in site_of]
f_all, ho_all, he_all, nl_all = fis(hi[:, allcols])
print(f"\nB. FIS POOLED across all stands")
print(f"  {'ALL':<8}{len(allcols):>6}{nl_all:>8}{ho_all:>9.4f}{he_all:>9.4f}{f_all:>9.4f}")

mw = float(np.mean(within)) if within else float("nan")
print(f"\n  mean within-stand FIS : {mw:+.4f}")
print(f"  pooled FIS            : {f_all:+.4f}")
print(f"  Wahlund contribution  : {f_all - mw:+.4f}")

print("\nC. depth dependence of the deficit (calling error vs null alleles)")
print(f"  {'depth bin':>12}{'genotypes':>12}{'Ho':>9}")
for lo, hi_ in [(6, 11), (12, 19), (20, 34), (35, 59), (60, 10**9)]:
    m = (D >= lo) & (D <= hi_) & ~np.isnan(G)
    if m.sum() > 1000:
        lab = f"{lo}-{hi_}" if hi_ < 10**9 else f"{lo}+"
        print(f"  {lab:>12}{int(m.sum()):>12}{float(np.mean(G[m] == 1)):>9.4f}")

print("\nD. is the deficit concentrated in some loci (null alleles) or spread (error)?")
n = np.sum(~np.isnan(hi), axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    p = np.nansum(hi, axis=1) / (2 * n)
    ho = np.nansum(hi == 1, axis=1) / n
    he = 2 * p * (1 - p)
    ok = (n >= 20) & (he > 0.1)
    fl = 1 - ho[ok] / he[ok]
print(f"  loci assessed: {ok.sum()}")
print(f"  FIS per locus: median {np.median(fl):+.3f}, "
      f"IQR {np.percentile(fl,25):+.3f} to {np.percentile(fl,75):+.3f}")
print(f"  loci with FIS > 0.9 (near-complete deficit): {100*np.mean(fl > 0.9):.1f}%")
print(f"  loci with FIS < 0.1 (behaving normally)   : {100*np.mean(fl < 0.1):.1f}%")

print("\nINTERPRETATION")
if not np.isnan(mw):
    if mw < 0.15 and f_all - mw > 0.10:
        print("  >>> WAHLUND. The deficit is largely an artefact of pooling structured stands.")
        print("  >>> Within stands the data behave normally; no correction needed for")
        print("  >>> within-stand analyses, and the low pooled Ho is expected, not damage.")
    elif mw > 0.3:
        print("  >>> NOT Wahlund. A large deficit persists WITHIN stands, so inbreeding, null")
        print("  >>> alleles or calling error dominate -- use panels C and D to choose.")
    else:
        print("  >>> MIXED. Wahlund contributes but a within-stand deficit remains.")
