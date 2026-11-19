#!/bin/bash
#SBATCH --job-name=arau-ipy-full
#SBATCH --partition=basic
#SBATCH --cpus-per-task=24   # n001 shares with other users; 24 fits alongside them
#SBATCH --mem=900G           # ample: the 30-library sweep peaked well under 96G
#SBATCH --nodelist=n001
#SBATCH --time=120:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/ipyfull_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/ipyfull_%j.err

# Full ipyrad assembly, run as a CO-PRIMARY alongside the Stacks assembly rather than as a
# diagnostic. The original reason for ranking Stacks first -- its de Bruijn paired-end contig
# assembly (Rochette et al. 2019) -- does not apply here, because R2 in this library is anchored
# rather than randomly sheared and the resulting contig is a fixed ~190 bp instead of 400-800 bp.
# With that advantage gone there is no basis to rank the two assemblers in advance, and Shafer
# et al. (2017) recommend running more than one regardless. Both are reported.
#
# Usage: sbatch ipyrad_full.sh <clust_threshold>
# The threshold comes from the replicate-concordance sweep, applied under the same pre-declared
# rule used for Stacks M/n. Do not choose it by inspecting downstream results.

set -uo pipefail
CT=${1:?usage: ipyrad_full.sh <clust_threshold>}
PROJ=/home/ecastron/araucaria
E=/opt/miniconda/miniconda-py312_24.3.0-0/envs/ipyrad
export PATH=$E/bin:$PATH

NAME="full_ct${CT/./}"
WORK=$PROJ/rad/ipyrad_full
LINKS=$WORK/r1_links
mkdir -p "$LINKS"; cd "$WORK"

# Own IPython dir: ipyrad's ipyparallel controller collides with any concurrent ipyrad run
# sharing ~/.ipython/profile_default, and the loser blocks forever at "loading reads" with no
# error. A sweep task lost 3.5 h to this before it was isolated.
export IPYTHONDIR="$WORK/.ipython_${NAME}"
mkdir -p "$IPYTHONDIR"

# Libraries entering the assembly: every cleaned library except the ones with too few reads to
# genotype. The floor is applied here, once, from the spot counts recorded at download, so the
# same set goes into both assemblers.
MIN_SPOTS=50000
awk -v m=$MIN_SPOTS 'NR>1 && $2>=m {print $1}' "$PROJ/metadata/library_spots.tsv" \
  > "$WORK/keep_libs.txt" 2>/dev/null \
  || ls "$PROJ/rad/cleaned/"*.1.fq.gz | xargs -n1 basename | sed 's/\.1\.fq\.gz$//' > "$WORK/keep_libs.txt"

while read -r s; do
  [[ -s "$PROJ/rad/cleaned/${s}.1.fq.gz" ]] && ln -sf "$PROJ/rad/cleaned/${s}.1.fq.gz" "$LINKS/${s}.fq.gz"
done < "$WORK/keep_libs.txt"
echo "libraries entering assembly: $(ls "$LINKS" | wc -l)"

P=$WORK/params-${NAME}.txt
cat > "$P" <<EOF
------- ipyrad params file (v.0.9.95)-------------------------------------
${NAME}                        ## [0] assembly_name
${WORK}                        ## [1] project_dir
                               ## [2] raw_fastq_path
                               ## [3] barcodes_path
${LINKS}/*.fq.gz               ## [4] sorted_fastq_path
denovo                         ## [5] assembly_method
                               ## [6] reference_sequence
rad                            ## [7] datatype
TGCAG,                         ## [8] restriction_overhang
5                              ## [9] max_low_qual_bases
33                             ## [10] phred_Qscore_offset
6                              ## [11] mindepth_statistical
6                              ## [12] mindepth_majrule
10000                          ## [13] maxdepth
${CT}                          ## [14] clust_threshold
0                              ## [15] max_barcode_mismatch
2                              ## [16] filter_adapters
35                             ## [17] filter_min_trim_len
2                              ## [18] max_alleles_consens
0.05                           ## [19] max_Ns_consens
0.05                           ## [20] max_Hs_consens
4                              ## [21] min_samples_locus
0.2                            ## [22] max_SNPs_locus
8                              ## [23] max_Indels_locus
0.5                            ## [24] max_shared_Hs_locus
0, 0, 0, 0                     ## [25] trim_reads
0, 0, 0, 0                     ## [26] trim_loci
p, s, v, k, g                  ## [27] output_formats
                               ## [28] pop_assign_file
                               ## [29] reference_as_filter
EOF

# min_samples_locus is left permissive at 4 here and the missingness threshold is applied later
# in the per-analysis SNP filtering, so that structure, Ne and relatedness can each take the
# filter they need from one assembly rather than forcing one compromise at assembly time.

echo "=== ipyrad full assembly, clust_threshold=${CT} ==="; date
if ! ipyrad -p "$P" -s 1234567 -c "$SLURM_CPUS_PER_TASK" -f > "$WORK/${NAME}.stdout" 2>&1; then
  echo "FAILED" >&2; tail -40 "$WORK/${NAME}.stdout" >&2; exit 1
fi
date

VCF="$WORK/${NAME}_outfiles/${NAME}.vcf"
[[ -s "$VCF" ]] || { echo "no VCF produced" >&2; exit 1; }
echo "DONE snps=$(grep -vc '^#' "$VCF") samples=$(grep -m1 '^#CHROM' "$VCF" | awk '{print NF-9}')"
cat "$WORK/${NAME}_outfiles/${NAME}_stats.txt" | head -40
