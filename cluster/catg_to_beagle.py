#!/usr/bin/env python3
"""Convert ipyrad CATG read counts into genotype likelihoods (Beagle format) for ngsRelate.

Why likelihoods rather than called genotypes: this dataset has NO technical replicates -- the
fifteen paired tree_ids turned out to be related trees, not replicates -- so the genotyping error
rate cannot be measured directly. Low-coverage sites are also demonstrably under-called
(heterozygosity 0.013 at 6-11x rising to 0.115 at 60x+). Likelihoods propagate that uncertainty
into the kinship estimate instead of forcing a hard call and discarding the evidence.

The ipyrad VCF carries GT:DP:CATG, where CATG holds raw per-base read counts. That is the
sufficient statistic for a genotype likelihood, so no re-assembly or BAM handling is needed.

MODEL. For a biallelic site with reference allele R and alternate A, and observed counts
n_R and n_A, under a per-read error rate e:

    P(read is A | genotype RR) = e
    P(read is A | genotype RA) = 0.5
    P(read is A | genotype AA) = 1 - e

so the likelihood of each genotype is binomial in n_A out of (n_R + n_A). Reads supporting
neither R nor A are dropped: at a biallelic site they are sequencing error or a third allele, and
counting them would inflate the apparent error rate.

e = 0.01 is used as a fixed per-read error. This is the one assumption in the conversion. It is
deliberately conservative -- higher than a typical Illumina per-base error after quality
filtering, so genotypes are called less confidently rather than more.

Beagle format: marker allele1 allele2 then three columns per individual (P(RR), P(RA), P(AA)),
normalised to sum to 1 per individual.
"""
import math, re, sys

VCF = sys.argv[1]
OUT = sys.argv[2]
ERR = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
BASES = "CATG"
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}     # Beagle allele coding

def binom_logpmf(k, n, p):
    if p <= 0:
        return 0.0 if k == 0 else -math.inf
    if p >= 1:
        return 0.0 if k == n else -math.inf
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log(1 - p))

n_out = n_skip = 0
with open(VCF) as fh, open(OUT, "w") as out:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = [re.sub(r"\.[12]$", "", s) for s in f[9:]]
            hdr = ["marker", "allele1", "allele2"]
            for s in samples:
                hdr += [s, s, s]
            out.write("\t".join(hdr) + "\n")
            continue
        ref, alt = f[3], f[4].split(",")[0]
        if ref not in BASES or alt not in BASES or ref == alt:
            n_skip += 1
            continue
        ri, ai = BASES.index(ref), BASES.index(alt)
        row = [f"{f[0]}_{f[1]}", str(CODE[ref]), str(CODE[alt])]
        for c in f[9:]:
            p = c.split(":")
            try:
                catg = [int(x) for x in p[2].split(",")]
            except (IndexError, ValueError):
                row += ["0.333333", "0.333333", "0.333333"]
                continue
            nr, na = catg[ri], catg[ai]
            n = nr + na
            if n == 0:
                row += ["0.333333", "0.333333", "0.333333"]      # no data: flat likelihood
                continue
            ls = [binom_logpmf(na, n, ERR),
                  binom_logpmf(na, n, 0.5),
                  binom_logpmf(na, n, 1 - ERR)]
            m = max(ls)
            ex = [math.exp(x - m) for x in ls]
            t = sum(ex)
            row += [f"{v / t:.6f}" for v in ex]
        out.write("\t".join(row) + "\n")
        n_out += 1

print(f"{VCF} -> {OUT}")
print(f"  {n_out} biallelic sites written, {n_skip} skipped (non-biallelic or indel)")
print(f"  {len(samples)} individuals, per-read error rate e = {ERR}")
