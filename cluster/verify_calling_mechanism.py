#!/usr/bin/env python3
"""What actually produces the replicate disagreements: the reads, or the caller?

The ipyrad VCF carries raw per-sample allele counts in the CATG field, so the reads underlying
every genotype call are available. That lets the mechanism be observed directly instead of
inferred from a signature -- which is where the last three diagnoses went wrong.

THREE MECHANISMS, THREE DISTINCT READ SIGNATURES, stated before running:

  1. THRESHOLD SUPPRESSION (caller's fault, FIXABLE)
     A homozygous call sits on reads that clearly show two alleles, e.g. 7 REF / 3 ALT. ipyrad
     calls a heterozygote only when the minor allele passes a likelihood threshold that depends
     on its per-sample heterozygosity estimate; below it, the majority allele wins. Two libraries
     splitting opposite ways then yield an opposite homozygote.
     PREDICTS: at disagreement sites, homozygous calls carry substantial minor-allele counts
     (>=2 reads, minor fraction well above 0), and opposite-homozygote pairs show BOTH libraries
     with real minority alleles.

  2. TRUE DROPOUT (sampling, PARTLY fixable by depth)
     The reads genuinely contain only one allele.
     PREDICTS: homozygous calls at disagreement sites have minor-allele count ~0.

  3. NON-HOMOLOGY (different loci, NOT fixable)
     Each library's reads are internally clean but for different sequences.
     PREDICTS: both libraries clean (minor count ~0) AND fixed for different alleles, with no
     depth dependence.

Mechanism 1 means the heterozygotes are still in the reads and genotype likelihoods computed
from CATG would recover them. Mechanism 3 means they were never there.
"""
import csv, re, sys
import numpy as np

PROJ = "/home/ecastron/araucaria"
V = sys.argv[1] if len(sys.argv) > 1 else (
    f"{PROJ}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf")

pairs = [(r["run_a"], r["run_b"]) for r in
         csv.DictReader(open(f"{PROJ}/metadata/replicate_pairs_validated.tsv"), delimiter="\t")
         if r["status"] == "validated"]

BASES = "CATG"          # order of the CATG format field

samples = []
agree_minor, dis_minor = [], []          # minor-allele fraction at homozygous calls
opp_pairs = []                           # (minor frac lib1, minor frac lib2) at opposite homs
hethom_pairs = []
n_sites = 0

with open(V) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            idx = {s: i for i, s in enumerate(samples)}
            continue
        ref, alt = f[3], f[4].split(",")[0]
        if ref not in BASES or alt not in BASES:
            continue
        ri, ai = BASES.index(ref), BASES.index(alt)
        cells = f[9:]
        n_sites += 1

        def counts(cell):
            p = cell.split(":")
            if len(p) < 3:
                return None, None, None
            g = p[0].replace("|", "/")
            try:
                catg = [int(x) for x in p[2].split(",")]
            except ValueError:
                return None, None, None
            return g, catg[ri], catg[ai]

        for a, b in pairs:
            if a not in idx or b not in idx:
                continue
            ga, ra, aa = counts(cells[idx[a]])
            gb, rb, ab = counts(cells[idx[b]])
            if ga is None or gb is None:
                continue
            if ga in ("./.", ".") or gb in ("./.", "."):
                continue
            ta, tb = ra + aa, rb + ab
            if ta < 6 or tb < 6:
                continue

            def hom(g):
                x, _, y = g.partition("/")
                return x == y

            def minor_frac(r_, a_):
                t = r_ + a_
                return min(r_, a_) / t if t else 0.0

            if ga == gb:
                if hom(ga):
                    agree_minor.append(minor_frac(ra, aa))
            else:
                if hom(ga):
                    dis_minor.append(minor_frac(ra, aa))
                if hom(gb):
                    dis_minor.append(minor_frac(rb, ab))
                if hom(ga) and hom(gb) and ga != gb:
                    opp_pairs.append((minor_frac(ra, aa), minor_frac(rb, ab),
                                      min(ra, aa), min(rb, ab)))
                elif hom(ga) != hom(gb):
                    hh = (ra, aa) if hom(ga) else (rb, ab)
                    hethom_pairs.append(minor_frac(*hh))

agree_minor = np.array(agree_minor)
dis_minor = np.array(dis_minor)
print(f"{V.split('/')[-1]}: {n_sites} sites, {len(pairs)} validated pairs\n")

print("MINOR-ALLELE FRACTION at HOMOZYGOUS calls (reads that disagree with the call)")
print(f"  {'context':<34}{'n':>9}{'mean':>9}{'median':>9}{'% with >=2 minor reads':>24}")
for lab, arr, cnts in (("sites where replicates AGREE", agree_minor, None),
                       ("sites where replicates DISAGREE", dis_minor, None)):
    if len(arr):
        print(f"  {lab:<34}{len(arr):>9}{arr.mean():>9.4f}{np.median(arr):>9.4f}"
              f"{100*np.mean(arr > 0):>23.1f}%")

if opp_pairs:
    o = np.array(opp_pairs)
    print(f"\nOPPOSITE-HOMOZYGOTE disagreements: {len(o)}")
    print(f"  minor-allele fraction, library A : mean {o[:,0].mean():.4f}, "
          f"median {np.median(o[:,0]):.4f}")
    print(f"  minor-allele fraction, library B : mean {o[:,1].mean():.4f}, "
          f"median {np.median(o[:,1]):.4f}")
    print(f"  minor-allele READ COUNT >= 2 in library A: {100*np.mean(o[:,2] >= 2):.1f}%")
    print(f"  minor-allele READ COUNT >= 2 in library B: {100*np.mean(o[:,3] >= 2):.1f}%")
    print(f"  BOTH libraries carry >= 2 minority reads : "
          f"{100*np.mean((o[:,2] >= 2) & (o[:,3] >= 2)):.1f}%")
    print(f"  BOTH libraries have ZERO minority reads  : "
          f"{100*np.mean((o[:,2] == 0) & (o[:,3] == 0)):.1f}%")

if hethom_pairs:
    h = np.array(hethom_pairs)
    print(f"\nHET-vs-HOM disagreements: {len(h)}")
    print(f"  minor fraction at the HOM call: mean {h.mean():.4f}, median {np.median(h):.4f}, "
          f"% with >=2 minor reads {100*np.mean(h > 0):.1f}%")

print("\nINTERPRETATION")
if len(opp_pairs):
    both_zero = np.mean((np.array(opp_pairs)[:, 2] == 0) & (np.array(opp_pairs)[:, 3] == 0))
    both_two = np.mean((np.array(opp_pairs)[:, 2] >= 2) & (np.array(opp_pairs)[:, 3] >= 2))
    if both_two > 0.3:
        print("  >>> THRESHOLD SUPPRESSION. The heterozygotes are present in the reads and the")
        print("  >>> caller is discarding them. Genotype likelihoods from CATG would recover")
        print("  >>> them; this is an assembly artefact, not a limit of the data.")
    elif both_zero > 0.7:
        print("  >>> TRUE DROPOUT or NON-HOMOLOGY. The minority alleles are absent from the")
        print("  >>> reads, so no re-calling can recover them.")
    else:
        print("  >>> MIXED. Both mechanisms contribute; quantify each before choosing a fix.")
