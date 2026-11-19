#!/bin/bash
#SBATCH --job-name=simctrl
#SBATCH --partition=basic
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/simctrl_%j.out
set -euo pipefail
module load Minimap2/minimap2
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
S=/home/ecastron/araucaria/results/simcontrol
cd "$S"
# small reference subset so indexing is seconds, not minutes
if [ ! -s simref.fa ]; then
  # NB: 'awk | head' triggers SIGPIPE under `set -o pipefail`; self-limit instead
  awk '$2>500000{print $1; n++} n>=120{exit}' "$D/v1.fasta.fai" > simcontigs.txt
  samtools faidx "$D/v1.fasta" $(tr '\n' ' ' < simcontigs.txt) > simref.fa
  samtools faidx simref.fa
fi
awk '{s+=$2} END{printf "sim reference: %d contigs, %.1f Mb\n", NR, s/1e6}' simref.fa.fai
minimap2 -x sr -d simref.mmi simref.fa 2>/dev/null

for MODE in diploid mixture; do
  echo
  echo "################################################################"
  echo "###  $MODE"
  [ "$MODE" = diploid ] && echo "###  TRUE heterozygotes -> statistic MUST return ~1.0x" \
                        || echo "###  3 haplotypes, random ratios -> MUST return high"
  echo "################################################################"
  python3 /home/ecastron/araucaria/scripts/ref/simulate_reads.py \
     --ref simref.fa --out sim_$MODE --mode $MODE --loci 20000 --seed 42
  minimap2 -ax sr -t 12 simref.mmi sim_${MODE}_1.fq.gz sim_${MODE}_2.fq.gz 2>/dev/null \
    | samtools sort -@ 4 -m 2G -o sim_$MODE.bam -
  samtools index sim_$MODE.bam
  samtools flagstat sim_$MODE.bam | grep -E "in total|[^y] mapped \("
  # IDENTICAL pileup flags and IDENTICAL script as every real-data run today
  samtools mpileup -Q 20 -q 30 -d 2000 -f simref.fa sim_$MODE.bam \
    | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
        --label "SIMULATED $MODE (known truth)" > sim_$MODE.balance.txt
  sed -n '1,14p' sim_$MODE.balance.txt
done
echo
echo "############ REAL DATA, for comparison ############"
awk '/^[0-9]/&&NF>=9{printf "  %-8s %s\n",$1,$8}' \
  /home/ecastron/araucaria/results/refmap/SRR11841615.q30.balance.txt
