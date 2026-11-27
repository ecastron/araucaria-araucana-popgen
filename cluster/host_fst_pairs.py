#!/usr/bin/env python
"""Pairwise Hudson F_ST between the 12 stands (66 pairs), with great-circle distances.

PRODUCES  results/derived/host_fst_pairs.tsv
          columns: stand_a stand_b km fst between nloci
          (copied to figdata/host_fst_pairs.tsv;
          read by figures/build_F2_journal.py, figures/qst_fst.py, figures/build_supp_journal.py,
          and referred to by scripts/abc/mc1_fst_ci.py and scripts/ref/fb_fst_pairs.py)
          and results/derived/stand_fst_matrix.tsv (12 x 12 symmetric matrix of the same
          F_ST values, diagonal 0, read by figures/build_neighbornet.R)

INPUTS    rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf   (via scripts/abc/load_observed.py)
          results/Q1/feems/keep.txt                                 (the 163 relatedness-pruned trees)
          metadata/popmap_177.tsv, metadata/manifest rad runs.csv   (stand and coordinates)

METHOD    Genotype dosages from GT (biallelic SNPs only, 163 trees, 78,026 SNPs), stand allele
          frequencies from called genotypes, Hudson F_ST per locus combined as a ratio of averages
          (Bhatia et al. 2013), loci used for a pair only where both stands have >= 6 called
          chromosomes. Distance is the haversine distance between stand centroids (mean lon/lat of
          the stand's retained trees). 'between' = 1 when one stand is Coastal (PN, TR, VA) and the
          other Andean. These are exactly the functions the ABC code uses for the observed vector
          (scripts/abc/sumstats.py: freqs_by_deme, hudson_fst; scripts/abc/chain.py: _km), so this
          is project code, not a published package; Hudson's estimator is implemented directly
          (no scikit-allel call).

PRINTS    163 trees x 78,026 SNPs; IBD r = 0.647; within-range n = 39, 0.079-0.130;
          between-range n = 27, 0.136-0.206 (all as quoted in the Results).
"""
import sys, itertools as it, numpy as np
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
from load_observed import load            # parses the cr30 VCF restricted to keep.txt
from sumstats import freqs_by_deme, hudson_fst
from chain import _km

OUTDIR = "/home/ecastron/araucaria/results/derived"
MINCHR = 6        # called chromosomes required in each stand of a pair

G, DEME, COORDS, COASTAL, stands = load()
print(f"genotype matrix: {G.shape[0]} trees x {G.shape[1]} SNPs | {len(stands)} stands")
P, N = freqs_by_deme(G, DEME, len(stands))

rows = []
for a, b in it.combinations(range(len(stands)), 2):
    m = (N[a] >= MINCHR) & (N[b] >= MINCHR)
    f = hudson_fst(P[a][m], P[b][m], N[a][m], N[b][m])
    rows.append((stands[a], stands[b], _km(COORDS[a], COORDS[b]), f,
                 int(COASTAL[a] != COASTAL[b]), int(m.sum())))

with open(f"{OUTDIR}/host_fst_pairs.tsv", "w") as fh:
    fh.write("stand_a\tstand_b\tkm\tfst\tbetween\tnloci\n")
    for r in rows:
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.6f}\t{r[4]}\t{r[5]}\n")

# 12 x 12 matrix for the NeighborNet; values are the 6-decimal strings of the pairs table so the
# two files can never disagree by rounding
M = {(r[0], r[1]): f"{r[3]:.6f}" for r in rows}
with open(f"{OUTDIR}/stand_fst_matrix.tsv", "w") as fh:
    fh.write("\t".join(stands) + "\n")
    for s in stands:
        fh.write("\t".join("0.000000" if s == t else M.get((s, t), M.get((t, s)))
                           for t in stands) + "\n")

# quantities quoted in the Results, for the record
f = np.array([r[3] for r in rows]); km = np.array([r[2] for r in rows])
bw = np.array([r[4] for r in rows], bool)
print(f"IBD r (F_ST vs log km) = {np.corrcoef(np.log(km), f)[0, 1]:.3f}")
print(f"within-range  n={(~bw).sum()}  {f[~bw].min():.3f}-{f[~bw].max():.3f}")
print(f"between-range n={bw.sum()}  {f[bw].min():.3f}-{f[bw].max():.3f}")
print(f"written -> {OUTDIR}/host_fst_pairs.tsv, stand_fst_matrix.tsv")
