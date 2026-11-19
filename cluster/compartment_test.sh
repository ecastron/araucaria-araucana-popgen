#!/bin/bash
#SBATCH --job-name=arau-comptest
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/comptest_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/comptest_%j.err

# Is the disputed 2018 compartment root or trunk? Decide from sequence content.
#
# Araucaria araucana forms arbuscular mycorrhizas (Glomeromycota) and forms NO
# ectomycorrhizas (Diehl et al. 2010; Andrade et al. 2000; Chavez et al. 2024 report 17
# Glomeromycota species from A. araucana rhizosphere at Nahuelbuta, Conguillio and
# Tolhuaca, which are three of our own sites). Arbuscular mycorrhizal fungi colonise
# ROOTS and do not occur in trunk or stem tissue. So:
#
#   disputed samples are ROOT  -> Glomeromycota present, strongly enriched over foliage
#   disputed samples are TRUNK -> Glomeromycota ~absent, indistinguishable from foliage
#
# The inference is differential (disputed vs branch/needle), not absolute abundance, so it
# is robust to ITS1F primer bias against Glomeromycota and to reference incompleteness.
#
# Foliar controls are depth-matched so that a difference cannot be an artefact of the
# disputed samples being the deepest ITS libraries in the dataset.

set -uo pipefail

PROJ=/home/ecastron/araucaria
REF=$PROJ/refs/sh_general_release_dynamic_19.02.2025.fasta
OUT=$PROJ/compartment_test
SUB=20000      # reads subsampled per library
mkdir -p "$OUT"

module load VSEARCH/vsearch-2.30.0 2>/dev/null || export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/vsearch/bin:$PATH
command -v vsearch >/dev/null || { echo "FATAL: vsearch not found" >&2; exit 1; }
command -v seqtk   >/dev/null || { echo "FATAL: seqtk not found" >&2; exit 1; }

classify () {          # $1 = run accession, $2 = group label
  local run="$1" grp="$2"
  local fq=$PROJ/raw/microbiome/${run}_1.fastq.gz
  [[ -s "$fq" ]] || { echo "  skip $run (no file)"; return; }
  local wk=$OUT/tmp_$run; mkdir -p "$wk"
  seqtk sample -s 42 "$fq" $SUB | seqtk seq -A - > "$wk/reads.fa" 2>/dev/null
  local n; n=$(grep -c "^>" "$wk/reads.fa" 2>/dev/null || echo 0)
  if [[ "$n" -lt 200 ]]; then echo "  skip $run (only $n reads)"; rm -rf "$wk"; return; fi
  vsearch --usearch_global "$wk/reads.fa" --db "$REF" \
          --id 0.80 --maxaccepts 1 --maxrejects 32 --top_hits_only \
          --blast6out "$wk/hits.tsv" --threads "$SLURM_CPUS_PER_TASK" \
          --quiet 2>/dev/null
  # Best-hit phylum per read, then per-sample phylum fractions.
  awk -v run="$run" -v grp="$grp" -v tot="$n" '
    { if (match($2, /p__[^;]*/)) { p=substr($2, RSTART+3, RLENGTH-3); cnt[p]++; hit++ } }
    END { for (p in cnt) printf "%s\t%s\t%s\t%d\t%d\t%d\n", run, grp, p, cnt[p], hit, tot }
  ' "$wk/hits.tsv" >> "$OUT/phylum_counts.tsv"
  echo "  done $run ($grp): $n reads, $(wc -l < "$wk/hits.tsv") hits"
  rm -rf "$wk"
}

printf "run\tgroup\tphylum\tcount\ttotal_hits\ttotal_reads\n" > "$OUT/phylum_counts.tsv"
while IFS=$'\t' read -r run grp; do
  [[ -z "${run:-}" || "$run" == "run" ]] && continue
  classify "$run" "$grp"
done < "$PROJ/metadata/compartment_test_samples.tsv"

echo "COMPARTMENT_TEST_DONE"
wc -l "$OUT/phylum_counts.tsv"
