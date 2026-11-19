#!/bin/bash
#SBATCH --job-name=arau-defline
#SBATCH --partition=basic
#SBATCH --array=1-188%12
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/defline_%A_%a.out
#SBATCH --error=/home/ecastron/araucaria/logs/defline_%A_%a.err

# Make read names Stacks-compatible. Two requirements, and both must hold:
#
#   1. The name before the mate suffix must be IDENTICAL across mates.
#      fasterq-dump writes "@SRR11841691.2 7_1101_4868_3776 length=91" and the length field
#      differs between mates (91 vs 101 bp), so the full deflines never match. tsv2bam then
#      reports "Failed to find any matching paired-end reads" - two hours into a run, after
#      ustacks, cstacks and sstacks have completed.
#
#   2. A mate suffix must still be PRESENT. Truncating to the first token alone removes it,
#      because process_radtags appended /1 and /2 to the last field, not the first. tsv2bam
#      then reports "Unrecognized paired-end read name format".
#
# So: keep the first token, strip any existing mate suffix, and re-append the correct one.
# Writing it that way makes the script idempotent - rerunning it is a no-op.
#
# Target form:   @SRR11841691.2/1   and   @SRR11841691.2/2

set -uo pipefail

PROJ=/home/ecastron/araucaria
POPMAP=$PROJ/metadata/popmap_all_libraries.tsv
DIR=$PROJ/rad/cleaned

run=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i {print $1}' "$POPMAP" | tr -d '\r')
[[ -n "${run:-}" ]] || { echo "FATAL: no run for index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

# Already in target form? (single token, ends in /1)
first=$(zcat "$DIR/${run}.1.fq.gz" 2>/dev/null | head -1)
if [[ "$first" != *" "* && "$first" == */1 ]]; then
  echo "SKIP  $run (already @name/1)"
  exit 0
fi

WORK=/scratch/$USER/defline_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

for m in 1 2; do
  src=$DIR/${run}.${m}.fq.gz
  [[ -s "$src" ]] || { echo "FAIL  $run missing $src" >&2; exit 1; }
  zcat "$src" \
    | awk -v m="$m" 'NR%4==1 {n=$1; sub(/\/[12]$/, "", n); print n "/" m; next} {print}' \
    | pigz -p "$SLURM_CPUS_PER_TASK" > "$WORK/${run}.${m}.fq.gz"
  a=$(zcat "$src" | wc -l); b=$(zcat "$WORK/${run}.${m}.fq.gz" | wc -l)
  if [[ "$a" -ne "$b" ]]; then
    echo "FAIL  $run mate $m line count changed: $a -> $b" >&2
    exit 1
  fi
done

# Confirm the two mates agree once the suffix is removed, before publishing.
n1=$(zcat "$WORK/${run}.1.fq.gz" | awk 'NR%4==1{print $1; exit}')
n2=$(zcat "$WORK/${run}.2.fq.gz" | awk 'NR%4==1{print $1; exit}')
if [[ "${n1%/1}" != "${n2%/2}" ]]; then
  echo "FAIL  $run mate names disagree: $n1 vs $n2" >&2
  exit 1
fi

mv "$WORK/${run}.1.fq.gz" "$WORK/${run}.2.fq.gz" "$DIR/"
echo "OK    $run  R1=$n1  R2=$n2"
