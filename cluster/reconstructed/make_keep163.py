#!/usr/bin/env python
"""The 163-tree working set: prune one member of each close-kin pair within a stand.

PRODUCES  results/reconstructed/keep.txt            two columns (FID IID = run accession), plink --keep
          results/reconstructed/keep_vcfnames.txt   same trees under their VCF sample names (SRRxxx.1)
          results/reconstructed/coords.txt          lon lat per kept tree, in keep order (feems input)
          (originals: results/Q1/feems/keep.txt, keep_vcfnames.txt, coords.txt, 24 Aug 2026; keep.txt
          is read by clock_da.py, f3_corrected.py, make_migrate_infile.py, abc/load_observed.py and
          others)
          It also prints the SNP count of the 163-tree matrix used downstream (78,026).

INPUTS    results/Q2/ngsrelate.tsv                     ngsRelate on the 177 libraries (run_q2.sh)
          metadata/popmap_177.tsv                      run -> stand
          rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf   called genotypes per library
          metadata/manifest rad runs.csv               coordinates

RULE      A pair is close kin when ngsRelate theta >= 0.25, both trees are in the same stand, and
          the pair shares at least MIN_SITES sites (nSites). From each such pair the tree with
          FEWER called genotypes in the cr30 VCF is dropped. 14 pairs qualify, no tree is in more
          than one pair, so the order of removal does not matter: 177 -> 163.

          The shared-sites floor is NOT stated in ANALYSIS_SPEC.md (which says only "14 within-
          stand pairs at ngsRelate theta >= 0.25"), but it is required: 19 within-stand pairs reach
          theta >= 0.25, and the 5 that were not pruned all rest on 102-203 shared sites, while the
          14 pruned ones rest on 1,592-7,368. Any floor from 204 to 1,592 gives the same set; 1,000
          is used here. Which member to drop is likewise inferred: in all 14 pairs the dropped tree
          has fewer called genotypes (it also has the lower summed depth in all 14, so the two
          criteria cannot be told apart from this output).

SNP COUNT  78,026 is the number of biallelic SNPs in full_ct091.cr30.vcf (78,802 SNP records minus
          776 multiallelic records), which is the filter scripts/abc/load_observed.py applies. No further
          locus filter is applied for the 163 trees (loci with no calls in them are kept as
          missing columns).

PROVENANCE  Reconstructed 18-19 Nov 2026 to replace lost inline code from 24 Aug 2026.
VERIFICATION  (run 18 Nov 2026 on McCartney, srun -p basic) keep.txt, keep_vcfnames.txt and
          coords.txt are each BYTE-IDENTICAL (cmp) to the originals in results/Q1/feems/. The run
          printed 14 kin pairs, 177 -> 163 trees, 78,802 SNP records, 78,026 biallelic.
"""
import csv, re, collections

P = "/home/ecastron/araucaria"
VCF = f"{P}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf"
OUT = f"{P}/results/reconstructed"
THETA, MIN_SITES = 0.25, 1000

stand = {r[0]: r[1] for r in csv.reader(open(f"{P}/metadata/popmap_177.tsv"), delimiter="\t")
         if len(r) >= 2}

# called genotypes per library, and the biallelic SNP count
calls = collections.Counter(); vcfname = {}; order = []; nrec = nbi = 0
for line in open(VCF):
    if line.startswith("##"):
        continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"):
        for s in f[9:]:
            base = re.sub(r"\.[12]$", "", s); vcfname[base] = s; order.append(base)
        continue
    nrec += 1
    if len(f[3]) == 1 and len(f[4]) == 1 and "," not in f[4]:
        nbi += 1
    for s, c in zip(order, f[9:]):
        if c.split(":", 1)[0] not in ("./.", "."):
            calls[s] += 1

pairs = []
for r in csv.DictReader(open(f"{P}/results/Q2/ngsrelate.tsv"), delimiter="\t"):
    a, b = r["ida"], r["idb"]
    if (float(r["theta"]) >= THETA and int(r["nSites"]) >= MIN_SITES
            and stand.get(a) is not None and stand.get(a) == stand.get(b)):
        pairs.append((a, b, float(r["theta"]), int(r["nSites"])))

drop = set()
for a, b, th, ns in sorted(pairs, key=lambda p: -p[2]):
    d = a if calls[a] < calls[b] else b
    print(f"{stand[a]}  {a} ({calls[a]:,} calls)  {b} ({calls[b]:,} calls)  "
          f"theta {th:.3f}  nSites {ns:,}  -> drop {d}")
    drop.add(d)
assert len(drop) == len(pairs), "a tree sits in two kin pairs; removal order would matter"

keep = [s for s in order if s in stand and s not in drop]
print(f"\n{len(pairs)} kin pairs | {len(order)} libraries -> {len(keep)} kept")
print(f"SNP records {nrec:,} | biallelic {nbi:,}")

lat, lon = {}, {}
for r in csv.DictReader(open(f"{P}/metadata/manifest rad runs.csv")):
    lat[r["run_accession"]] = r["latitude"]; lon[r["run_accession"]] = r["longitude"]

with open(f"{OUT}/keep.txt", "w") as fh:
    fh.writelines(f"{s}\t{s}\n" for s in keep)
with open(f"{OUT}/keep_vcfnames.txt", "w") as fh:
    fh.writelines(f"{vcfname[s]}\t{vcfname[s]}\n" for s in keep)
with open(f"{OUT}/coords.txt", "w") as fh:
    fh.writelines(f"{lon[s]} {lat[s]}\n" for s in keep)
print(f"written -> {OUT}/keep.txt, keep_vcfnames.txt, coords.txt")
