#!/bin/bash
#SBATCH --job-name=prepref
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/prepref_%j.out
set -euo pipefail
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
cd "$D"
echo "node $(hostname) start $(date)"
gzip -t v1.fasta.gz && echo "gzip integrity OK"
# samtools faidx needs plain or bgzip'd fasta; the download is plain gzip.
zcat v1.fasta.gz > v1.fasta
samtools faidx v1.fasta
echo "=== assembly stats ==="
awk '{b+=$2; n++; if($2>m)m=$2} END{printf "contigs: %d\ntotal  : %d bp (%.2f Gb)\nlongest: %d bp\n",n,b,b/1e9,m}' v1.fasta.fai
sort -k2,2nr v1.fasta.fai | awk '{b+=$2; L[NR]=$2; T=T+0} END{ }' 
awk '{print $2}' v1.fasta.fai | sort -nr | awk '{a[NR]=$1; s+=$1} END{h=s/2; c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=h){printf "N50    : %d bp\n",a[i]; break}}}'
ls -lh v1.fasta v1.fasta.fai
echo "done $(date)"
