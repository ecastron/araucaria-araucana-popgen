#!/bin/bash
#SBATCH --job-name=lanes
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/lanes_%j.out
set -uo pipefail
R=/home/ecastron/araucaria/raw/rad
O=/home/ecastron/araucaria/results/refmap
: > $O/lane_assignment.tsv
for f in $R/*_1.fastq.gz; do
  S=$(basename $f _1.fastq.gz)
  # defline: @SRR11841691.2 7_1101_4868_3776 length=91  -> field2 = lane_tile_x_y
  L=$(zcat "$f" 2>/dev/null | head -20000 \
      | awk 'NR%4==1{split($2,a,"_"); if(a[1]!="") print a[1]}' | sort | uniq -c | sort -rn | head -1)
  N=$(echo "$L" | awk '{print $2}'); C=$(echo "$L" | awk '{print $1}')
  printf "%s\t%s\t%s\n" "$S" "$N" "$C" >> $O/lane_assignment.tsv
done
echo "libraries processed: $(wc -l < $O/lane_assignment.tsv)"
echo "=== lane distribution ==="
awk -F'\t' '{n[$2]++} END{for(k in n) printf "  lane %-6s %d libraries\n",k,n[k]}' $O/lane_assignment.tsv
echo "=== any library spanning lanes (purity of top lane) ==="
awk -F'\t' '$3<5000{print "  MIXED: "$0}' $O/lane_assignment.tsv | head
