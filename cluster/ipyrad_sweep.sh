#!/bin/bash
#SBATCH --job-name=arau-ipyrad
#SBATCH --partition=basic
#SBATCH --array=1-4%2
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=36:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/ipyrad_%A_%a.out
#SBATCH --error=/home/ecastron/araucaria/logs/ipyrad_%A_%a.err

# ipyrad sweep over clustering threshold, on the same 30 replicate-pair libraries used for
# the Stacks sweep, scored by the same criterion: genotyping error between the 13 usable
# replicate library pairs.
#
# Purpose is diagnostic. Stacks gives 21-31% discordance between two libraries of the SAME
# tree under every treatment tried (depth, HDplot paralog filtering, PCR-duplicate removal,
# M = 1-8). ipyrad clusters within and between samples rather than matching samples to a
# catalog, so if the discordance is a catalog-matching artefact it should largely disappear.
# If it does not, the data cannot support individual genotyping and Q6 needs genotype
# likelihoods rather than called genotypes.
#
# R1 only. R2 in this library is anchored, not randomly sheared (1.38 distinct R2 per distinct
# R1; 77% of R1 have exactly one R2 partner), so it provides no shear diversity and only 18%
# of SNPs. Excluding it keeps the assembler comparison clean; whether it earns inclusion is a
# separate question to settle after this one.
#
# 0.85 reproduces Varas-Myrik et al. (2022) as a reference point. min_samples_locus is 8 of 30
# (~27%), matching the r = 0.30 working point the Stacks locus-recovery curve identified;
# their 4-of-190 (2%) was tuned for structure, not for relatedness.

set -uo pipefail
PROJ=/home/ecastron/araucaria
E=/opt/miniconda/miniconda-py312_24.3.0-0/envs/ipyrad
export PATH=$E/bin:$PATH

THRESHOLDS=(0.85 0.88 0.91 0.94)
CT=${THRESHOLDS[$((SLURM_ARRAY_TASK_ID - 1))]}
NAME="ct${CT/./}"
WORK=$PROJ/rad/ipyrad
mkdir -p "$WORK"
cd "$WORK"

# Each array task gets its OWN IPython directory. ipyrad drives an ipyparallel cluster, and
# concurrent tasks sharing ~/.ipython/profile_default write competing ipcontroller files: one
# task captures the cluster and the others block forever at "loading reads" with no error.
# Observed here -- ct085 sat at 0% for 3.5 h alongside ct088 before this was isolated.
export IPYTHONDIR="$WORK/.ipython_${NAME}"
mkdir -p "$IPYTHONDIR"

# ipyrad names samples from the filename; our cleaned files carry a .1/.2 mate suffix that
# would become part of the sample name, so expose R1 under a clean name via symlinks.
LINKS=$WORK/r1_links
if [[ ! -d "$LINKS" ]]; then
  mkdir -p "$LINKS"
  while read -r s pop; do
    ln -sf "$PROJ/rad/cleaned/${s}.1.fq.gz" "$LINKS/${s}.fq.gz"
  done < "$PROJ/metadata/popmap_replicates.tsv"
fi
echo "libraries linked: $(ls $LINKS | wc -l)"

P=$WORK/params-${NAME}.txt
cat > "$P" <<EOF
------- ipyrad params file (v.0.9)ateway ---------------------------------------
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
8                              ## [21] min_samples_locus
0.2                            ## [22] max_SNPs_locus
8                              ## [23] max_Indels_locus
0.5                            ## [24] max_shared_Hs_locus
0, 0, 0, 0                     ## [25] trim_reads
0, 0, 0, 0                     ## [26] trim_loci
p, s, v                        ## [27] output_formats
                               ## [28] pop_assign_file
                               ## [29] reference_as_filter
EOF
# The header line above must start with the literal marker ipyrad expects.
sed -i '1s/.*/------- ipyrad params file (v.0.9.95)-------------------------------------/' "$P"

echo "=== ipyrad clust_threshold=${CT} ==="
date
if ! ipyrad -p "$P" -s 1234567 -c "$SLURM_CPUS_PER_TASK" -f > "$WORK/${NAME}.stdout" 2>&1; then
  echo "FAIL ct=$CT" >&2
  tail -30 "$WORK/${NAME}.stdout" >&2
  exit 1
fi
date

VCF="$WORK/${NAME}_outfiles/${NAME}.vcf"
[[ -s "$VCF" ]] || { echo "FAIL ct=$CT no VCF at $VCF" >&2; ls -R "$WORK/${NAME}_outfiles" >&2; exit 1; }
echo "OK ct=$CT snps=$(grep -vc '^#' "$VCF")"
