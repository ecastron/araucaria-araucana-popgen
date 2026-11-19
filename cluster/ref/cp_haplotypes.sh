#!/bin/bash
#SBATCH --job-name=cphap
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/cphap_%j.out
set -eo pipefail
module load Minimap2/minimap2
module load SAMTOOLS/samtools.1.20
module load Bcftools/1.20
D=/home/ecastron/araucaria/refs/araucana
S=/home/ecastron/araucaria/results/sv/cp
R=/home/ecastron/araucaria/raw/rad
cd "$S"
cp "$D/v1_chloroplast.fasta" cp_ref.fa
samtools faidx cp_ref.fa
minimap2 -x sr -d cp.mmi cp_ref.fa
echo "chloroplast: $(awk '/^>/{next}{n+=length($0)}END{print n}' cp_ref.fa) bp"
: > cp_coverage.tsv
# all 177 libraries in the analysis set -- the reference is 146 kb so this is cheap
for SRR in $(awk '{print $1}' /home/ecastron/araucaria/metadata/popmap_177.tsv); do
  [ -f "$R/${SRR}_1.fastq.gz" ] || continue
  minimap2 -ax sr -t 16 cp.mmi "$R/${SRR}_1.fastq.gz" "$R/${SRR}_2.fastq.gz" 2>/dev/null \
    | samtools sort -@ 4 -o "${SRR}.cp.bam" -
  samtools index "${SRR}.cp.bam"
  mp=$(samtools view -c -F 4 "${SRR}.cp.bam")
  read cvg mdep <<< $(samtools depth -a "${SRR}.cp.bam" | awk '{if($3>0)c++; s+=$3} END{printf "%.2f %.1f", 100*c/NR, s/NR}')
  printf "%s\t%s\t%s\t%s\n" "$SRR" "$mp" "$cvg" "$mdep" >> cp_coverage.tsv
done
echo "=== chloroplast yield across libraries ==="
sort -k2,2nr cp_coverage.tsv | head -5
echo "..."
sort -k2,2nr cp_coverage.tsv | tail -3
awk '{n++; r+=$2; c+=$3; d+=$4} END{printf "\nlibraries %d | mean mapped reads %.0f | mean %% sites covered %.1f%% | mean depth %.1fx\n", n, r/n, c/n, d/n}' cp_coverage.tsv
