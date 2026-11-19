#!/bin/bash
#SBATCH --job-name=arau-verify
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/verify_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/verify_%j.err

# Independent check that the downloaded FASTQ are complete and uncorrupted.
#
# Three tests per run:
#   1. gzip integrity  -- catches truncated or partially written files
#   2. mate concordance -- _1 and _2 must contain the same number of reads
#   3. spot concordance -- read count must equal the spot count NCBI reports
#
# Test 3 is the strongest: it ties the local file back to the archive's own record of
# how many reads the run contains, so a silently short file cannot pass.

set -uo pipefail

PROJ=/home/ecastron/araucaria
EXPECTED=$PROJ/metadata/expected_spots.tsv
OUT=$PROJ/logs/verify_report.tsv

check_one() {
  local ds="$1" run="$2" exp="$3"
  local f1="$PROJ/raw/$ds/${run}_1.fastq.gz"
  local f2="$PROJ/raw/$ds/${run}_2.fastq.gz"

  if ! pigz -t "$f1" 2>/dev/null; then echo -e "$run\t$ds\tGZIP_FAIL_R1\t-\t-\t$exp"; return; fi
  if ! pigz -t "$f2" 2>/dev/null; then echo -e "$run\t$ds\tGZIP_FAIL_R2\t-\t-\t$exp"; return; fi

  # 4 lines per read; count lines once per mate.
  local n1 n2
  n1=$(( $(pigz -dc "$f1" | wc -l) / 4 ))
  n2=$(( $(pigz -dc "$f2" | wc -l) / 4 ))

  local status="OK"
  [[ "$n1" -ne "$n2" ]] && status="MATE_MISMATCH"
  [[ "$status" == "OK" && "$exp" != "NA" && "$n1" -ne "$exp" ]] && status="SPOT_MISMATCH"
  echo -e "$run\t$ds\t$status\t$n1\t$n2\t$exp"
}
export -f check_one
export PROJ

printf "run\tdataset\tstatus\treads_R1\treads_R2\texpected_spots\n" > "$OUT"
# GNU parallel is installed on the login node but NOT on the compute nodes, so use
# xargs -P instead. export -f above makes check_one visible to the bash children.
tail -n +2 "$EXPECTED" | tr '\t' ' ' \
  | xargs -P "$SLURM_CPUS_PER_TASK" -n 3 bash -c 'check_one "$0" "$1" "$2"' \
  >> "$OUT"

echo "=== verification summary ==="
tail -n +2 "$OUT" | awk -F'\t' '{print $3}' | sort | uniq -c
echo "=== any non-OK rows ==="
tail -n +2 "$OUT" | awk -F'\t' '$3!="OK"' | head -20
