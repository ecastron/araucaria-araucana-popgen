#!/usr/bin/env python3
"""Assemble figdata/het_diagnosis.json (Figure S2 B, C, E and the heterozygote-deficit record).

PRODUCES  results/derived/het_diagnosis.json
          (copied to figdata/het_diagnosis.json; read by figures/build_supp_journal.py for S2B, S2C)

METHOD    Each number is read from the output file of the script that produced it. Every value is
          tagged in SOURCES below; values with no producer output are fixed values and are
          marked FIXED.

INPUTS (all on McCartney, under /home/ecastron/araucaria)
  results/derived/het_diag/*.txt   stdout of scripts/ref/{two_phenomena, rarefy_and_cv,
                                         locus_depth_filter, cnv_filter, hwe_per_stand, hwe_by_maf}.py,
                                         run on the cr30 VCF by het_diagnosis_inputs.sbatch
  results/Q1/ct_compare/ct_compare.txt   median ngsRelate F per clustering threshold (job 40007)
  results/polyrad/polyrad_F_by_depth_fixed.csv   polyRAD de novo F by depth (= figdata copy, md5)
  results/polyrad_fb/fb_F_by_depth.csv   FreeBayes + polyRAD F by depth (= figdata copy, md5)
  logs/hobsfilt_40316.out                H_obs ceilings (hobs_paralog_filter.R stdout)
  results/derived/het_diag/figdata_locus_filter_F.tsv   copy of the local figdata/locus_filter_F.tsv
                                         (structural filters)
  results/refmap/hdplot.tsv              HDplot per-SNP H and D
  logs/simF2_39642.out, logs/simW50_39644.out   allele-dropout and Wahlund+dropout ceilings
  logs/dosetest_41088.out                hemizygosity test
  results/sv/fb_fst_indiv.tsv, results/figdata/host_fst_pairs.tsv   de novo vs reference F_ST

"""
import csv, json, re, math

P = "/home/ecastron/araucaria"
H = f"{P}/results/derived/het_diag"
OUT = f"{P}/results/derived/het_diagnosis.json"
R3 = lambda x: round(x + 0.0, 3)
txt = lambda p: open(p).read()
SOURCES = {}


def table_rows(text, header_regex):
    """lines after a header line, until a blank line"""
    lines = text.splitlines()
    for i, l in enumerate(lines):
        if re.search(header_regex, l):
            out = []
            for l2 in lines[i + 1:]:
                if not l2.strip():
                    break
                out.append(l2)
            return out
    raise KeyError(header_regex)


def spread_after(text, header_regex):
    lines = text.splitlines()
    i = next(k for k, l in enumerate(lines) if re.search(header_regex, l))
    for l in lines[i:]:
        m = re.search(r"spread:?\s+([-\d.]+)", l)
        if m:
            return float(m.group(1))


# ---- observed F: median ngsRelate F at ct 0.91; clustering test at 0.85 --------------------------
ct = {r.split()[0]: float(r.split()[4]) for r in txt(f"{P}/results/Q1/ct_compare/ct_compare.txt").splitlines()
      if re.match(r"0\.(85|88|91)\s+0\.\d+\s+\d+\s+[\d.]+\s+[\d.]+$", r)}   # the ipyradH/SNPs/TsTv/medianF table
observed_F = ct["0.91"]
SOURCES["observed_F"] = "results/Q1/ct_compare/ct_compare.txt (ct 0.91 median F)"

# ---- F by genotype depth, within stand (two_phenomena.py) -----------------------------------------
tp = [l.split() for l in txt(f"{H}/two_phenomena.txt").splitlines()
      if re.match(r"^(6-10|10-15|15-20|20-30|30-50|50-100|100\+)\s", l)]
F_by_depth = {"bins": [r[0] for r in tp], "n": [int(r[1]) for r in tp],
              "F_within_stand": [float(r[2]) for r in tp],
              "balance_spread": [float(r[4]) for r in tp]}
SOURCES["F_by_genotype_depth"] = "scripts/ref/two_phenomena.py stdout"

# ---- rarefaction to 30 reads (rarefy_and_cv.py) ---------------------------------------------------
rc = txt(f"{H}/rarefy_and_cv.txt")
nat = [l.split() for l in table_rows(rc, r"\(0\) NO FILTER")[1:] if l.startswith("orig")]
rar = [l.split() for l in table_rows(rc, r"\(1\) RAREFIED")[1:] if l.startswith("orig")]
rarefaction = {"bins": [r[1] for r in nat], "native": [R3(float(r[4])) for r in nat],
               "rarefied": [R3(float(r[4])) for r in rar]}
spread_rarefied = spread_after(rc, r"\(1\) RAREFIED")
spread_cv_kept = spread_after(rc, r"\(2\) DEPTH-CV FILTER: loci KEPT")
m = re.search(r"loci evaluated: (\d+) \| high depth-CV \(>1.0\) dropped: (\d+)", rc)
n_cv_loci, n_cv_drop = int(m.group(1)), int(m.group(2))
SOURCES["rarefaction_to_30_reads"] = "scripts/ref/rarefy_and_cv.py stdout, sections (0) and (1)"

# ---- pipelines: polyRAD de novo, FreeBayes + polyRAD -----------------------------------------------
def fcsv(path):
    return {r["depth"]: float(r["F_IS"]) for r in csv.DictReader(open(path))}
pr = fcsv(f"{P}/results/polyrad/polyrad_F_by_depth_fixed.csv")
fb = fcsv(f"{P}/results/polyrad_fb/fb_F_by_depth.csv")
PB = ["6-10", "20-30", "50-100", "100+"]
ALL7 = ["6-10", "10-15", "15-20", "20-30", "30-50", "50-100", "100+"]
# FIXED: the ipyrad binomial curve on the polyRAD locus set is the list in
# figures/build_polyrad_F.py (Fi = 0.8694 ... -0.3363); cluster/s2e_polyrad_route.R recomputes it.
IPYRAD_BINOMIAL = dict(zip(ALL7, [0.8694, 0.8153, 0.7847, 0.6512, 0.3061, -0.1960, -0.3363]))
pipelines = {"bins": PB,
             "ipyrad binomial (de novo)": [R3(IPYRAD_BINOMIAL[b]) for b in PB],
             "polyRAD beta-binomial (de novo)": [R3(pr[b]) for b in PB],
             "FreeBayes + polyRAD (reference)": [R3(fb[b]) for b in PB]}
spr = lambda d: max(d[b] for b in ALL7) - min(d[b] for b in ALL7)
SOURCES["pipelines"] = ("polyrad_F_by_depth_fixed.csv, fb_F_by_depth.csv; ipyrad binomial FIXED "
                        "(constant from figures/build_polyrad_F.py)")

# ---- the intervention table -----------------------------------------------------------------------
def tsv_spread(path, label_prefix):
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["filter"].startswith(label_prefix):
            return float(r["spread"]), r
# H_obs ceilings: the producer's own stdout (hobs_paralog_filter.R, job 40316)
HO = {}
for l in txt(f"{P}/logs/hobsfilt_40316.out").splitlines():
    m = re.match(r"^(no Hobs filter|Hobs <= [\d.]+)\s+(\d+)\s.*spread ([-\d.]+)", l)
    if m:
        HO[m.group(1)] = (int(m.group(2)), float(m.group(3)))
hobs_base, hobs20 = HO["no Hobs filter"][1], HO["Hobs <= 0.20"][1]
hobs20_drop = round(100 * (1 - HO["Hobs <= 0.20"][0] / HO["no Hobs filter"][0]))
# structural filters: figdata/locus_filter_F.tsv (written by cluster/s2e_polyrad_route.R), copied in as an input
struct_mod, _ = tsv_spread(f"{H}/figdata_locus_filter_F.tsv", "cv<=0.80 depth<=95")

ld = txt(f"{H}/locus_depth_filter.txt")
m = re.search(r"mean depth <= 1\.25x median.*?loci kept\s+\d+ \((\d+)%\).*?spread ([-\d.]+)", ld, re.S)
ld_kept_pct, ld_spread = int(m.group(1)), float(m.group(2))

cn = txt(f"{H}/cnv_filter.txt")
cnv_15 = spread_after(cn, r"keeping only ratio <= 1\.5")
cnv_ratio = spread_after(cn, r"A: F by DEPTH RATIO")
cnv_abs = spread_after(cn, r"reference: F by ABSOLUTE depth")

spread = [
    ["none - unfiltered baseline", R3(hobs_base), "baseline"],
    ["rarefy every genotype to 30 reads", R3(spread_rarefied), "test"],
    [f"locus mean depth <= 1.25x median ({100 - ld_kept_pct}% dropped)", R3(ld_spread), "test"],
    [f"H_obs ceiling 0.20 ({hobs20_drop}% of loci dropped)", R3(hobs20), "test"],
    ["beta-binomial overdispersion, polyRAD d = 5", R3(spr(pr)), "test"],
    [f"depth-CV filter ({n_cv_drop} of {n_cv_loci:,} loci)", R3(spread_cv_kept), "test"],
    ["independent pipeline: bowtie2 + FreeBayes", R3(spr(fb)), "test"],
    ["locus-relative depth (CNV) filter <= 1.5", R3(cnv_15), "test"],
    ["structural locus filter, moderate setting", R3(struct_mod), "worse"],
]
SOURCES["spread_after_intervention"] = ("hobsfilt_40316.out, rarefy_and_cv.py, locus_depth_filter.py, "
                                        "polyRAD/FB csv, cnv_filter.py, locus_filter_F.tsv")

# ---- mechanisms -----------------------------------------------------------------------------------
pooled = float(re.search(r"^ALL\s+\S+\s+\S+\s+\S+\s+([-\d.]+)", txt(f"{H}/hwe_by_maf.txt"), re.M).group(1))
within = float(re.search(r"^MEAN\s+\S+\s+\S+\s+\S+\s+([-\d.]+)", txt(f"{H}/hwe_per_stand.txt"), re.M).group(1))
drop50 = float(re.search(r"drop50: .*median F = ([\d.]+)", txt(f"{P}/logs/simF2_39642.out")).group(1))
wd50 = float(re.search(r"wahlund_drop50: .*median F = ([\d.]+)", txt(f"{P}/logs/simW50_39644.out")).group(1))
hd = list(csv.DictReader(open(f"{P}/results/refmap/hdplot.tsv"), delimiter="\t"))
hd_flag = sum(1 for r in hd if float(r["H"]) > 0.55 and abs(float(r["D"])) > 7)
hemi = float(re.search(r"H_obs ratio high/low = ([\d.]+)", txt(f"{P}/logs/dosetest_41088.out")).group(1))
single = sum(1 for l in open(f"{P}/results/sv/loci_copies.tsv") if l.rstrip().endswith("\t1"))
ib = spr(IPYRAD_BINOMIAL)

mechanisms = [
    ["Wahlund effect", f"measured directly: pooled {pooled:.3f} - within-stand {within:.3f}",
     f"contributes F = {pooled - within:.3f} of {observed_F:.3f}"],
    ["Allele dropout", "arithmetic ceiling",
     f"cannot exceed {drop50:.3f}; with Wahlund maxed, {wd50:.3f} < {observed_F:.3f} observed"],
    ["Clustering oversplitting", "clust_threshold 0.91 -> 0.85",
     f"F ROSE, {ct['0.91']:.3f} -> {ct['0.85']:.3f}"],
    ["Paralogous loci", "HDplot", f"{hd_flag} of {len(hd):,} loci flagged"],
    # FIXED text rows: these tests print no number used here
    ["PCR duplication", "random-barcode analysis", "refuted"],
    ["Nuclear plastid inserts", "organellar mapping", "refuted"],
    ["Sequencing depth itself", "rarefy all genotypes to 30 reads",
     f"F still spans {rarefaction['rarefied'][0]:+.3f} to {rarefaction['rarefied'][-1]:+.3f}"],
    ["Aligner choice", "minimap2 / bwa-mem2 / bowtie2", "identical depth curves"],
    ["De novo assembly artifact", "bowtie2 + FreeBayes + polyRAD on the reference genome",
     f"reproduces, spread {spr(fb):.3f} vs {spr(pr):.3f} and {ib:.3f}"],
    ["Genotype-model misfit", "beta-binomial overdispersion d = 5",
     f"spread {spr(pr):.3f} vs {ib:.3f} binomial"],
    ["Segregating copy number", "locus-relative depth",
     f"predicts F WORSE than raw depth ({cnv_ratio:.3f} vs {cnv_abs:.3f})"],
    # FIXED: copy-class spreads from cluster/ref/F_by_copynumber.py (run by s2e_mechanisms.sbatch)
    ["Copy number, measured", f"{single:,} RAD loci mapped to the 19.3 Gb genome",
     "F identical across 1, 2, 3-5 and 6+ copy classes (spread 1.404-1.448)"],
    ["Hemizygosity / dosage", "per-individual bedcov over single-copy loci, 177 BAMs",
     f"H_obs ratio high/low = {hemi:.3f}; hemizygosity predicts ~4"],
    # SD numbers: selfaln + F_by_sd (fsd_41134.out holds the per-class F; the pair count and % are
    # FIXED values from the self-alignment)
    ["Segmental duplication", "genome self-alignment, 83,424 SD pairs, 8.8% of genome",
     "only 1.1% of RAD loci in SDs; SD tracks non-SD to within 0.02"],
]

# ---- consequences ----------------------------------------------------------------------------------
den = {(r["stand_a"], r["stand_b"]): float(r["fst"])
       for r in csv.DictReader(open(f"{P}/results/figdata/host_fst_pairs.tsv"), delimiter="\t")}
ref = {(r["stand_a"], r["stand_b"]): float(r["fst"])
       for r in csv.DictReader(open(f"{P}/results/sv/fb_fst_indiv.tsv"), delimiter="\t")}
k = sorted(set(den) & set(ref)); x = [den[i] for i in k]; y = [ref[i] for i in k]
mx, my = sum(x) / len(x), sum(y) / len(y)
r_struct = (sum((a - mx) * (b - my) for a, b in zip(x, y))
            / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)))
consequences = {
    "structure_robust_r": round(r_struct, 4),
    "homozygous_concordance": "98.3-99.5%",      # FIXED value
    "note": ("Homozygous calls are concordant with their read counts, so allele-frequency and "
             "structure-based statistics remain informative; individual heterozygosity statistics "
             "are not reported."),
}

out = {
    "_source": ("Assembled by cluster/het_diagnosis_assemble.py from "
                "the producer outputs named in its SOURCES; fixed values marked FIXED there."),
    "observed_F": observed_F,
    "F_by_genotype_depth": F_by_depth,
    "rarefaction_to_30_reads": rarefaction,
    "pipelines": pipelines,
    "spread_after_intervention": spread,
    "mechanisms": mechanisms,
    "consequences": consequences,
}
json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)
print(json.dumps(SOURCES, indent=1))
print(f"written -> {OUT}")
