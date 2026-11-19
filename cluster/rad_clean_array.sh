#!/bin/bash
#SBATCH --job-name=arau-radclean
#SBATCH --partition=basic
#SBATCH --array=1-188%20
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/clean_%A_%a.out
#SBATCH --error=/home/ecastron/araucaria/logs/clean_%A_%a.err

# Step 2, stage 1: clean one RAD library per array task.
#
# The runs are already demultiplexed (one run per library, barcodes stripped by the
# submitter), so -b is omitted and process_radtags is used purely to quality-filter and
# to verify the PstI cut site. R1 begins with the TGCAG remnant, confirmed by inspection.
#
#   -e pstI  check the cut site on the single-end read
#   -c       drop reads containing uncalled bases
#   -q       drop low-quality reads (sliding window)
#   -r       rescue RAD-tag cut sites with recoverable errors
#   -D       keep discards so the loss can be audited rather than assumed
#
# Idempotent: a library whose cleaned mates already exist is skipped.

set -uo pipefail

PROJ=/home/ecastron/araucaria
POPMAP=$PROJ/metadata/popmap_all_libraries.tsv
DEST=$PROJ/rad/cleaned

module load Stacks/stacks-2.68

run=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i {print $1}' "$POPMAP" | tr -d '\r')
if [[ -z "${run:-}" ]]; then
  echo "FATAL: no run for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

mkdir -p "$DEST"
if [[ -s "$DEST/${run}.1.fq.gz" && -s "$DEST/${run}.2.fq.gz" ]]; then
  echo "SKIP  $run (already cleaned)"
  exit 0
fi

R1=$PROJ/raw/rad/${run}_1.fastq.gz
R2=$PROJ/raw/rad/${run}_2.fastq.gz
for f in "$R1" "$R2"; do
  [[ -s "$f" ]] || { echo "FAIL  $run missing input $f" >&2; exit 1; }
done

# Work on node-local NVMe, then publish, so a killed task cannot leave a half-written
# pair that the skip check above would treat as complete.
WORK=/scratch/$USER/radclean_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

if ! process_radtags -1 "$R1" -2 "$R2" -o "$WORK" -e pstI -c -q -r -D \
       --basename "$run" --threads "$SLURM_CPUS_PER_TASK" > "$WORK/prt.log" 2>&1; then
  echo "FAIL  $run process_radtags failed" >&2
  tail -15 "$WORK/prt.log" >&2
  exit 1
fi

if [[ ! -s "$WORK/${run}.1.fq.gz" || ! -s "$WORK/${run}.2.fq.gz" ]]; then
  echo "FAIL  $run expected cleaned pair missing" >&2
  ls -la "$WORK" >&2
  exit 1
fi

mv "$WORK/${run}".*.fq.gz "$DEST/"
# Keep the per-library log; the retained/discarded counts feed the stage-1 QC summary.
cp "$WORK/prt.log" "$PROJ/logs/prt_${run}.log" 2>/dev/null
[[ -f "$WORK/process_radtags.log" ]] && cp "$WORK/process_radtags.log" "$PROJ/logs/prt_stats_${run}.log"

kept=$(zcat "$DEST/${run}.1.fq.gz" | awk 'NR%4==1' | wc -l)
echo "OK    $run kept_read_pairs=$kept"
