#!/bin/bash
#SBATCH --job-name=dlgenome
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dlgenome_%j.out
set -euo pipefail
B=http://plantaeviz.tomsbiolab.com:8001/easygdb/araucaria_db/downloads
D=/home/ecastron/araucaria/refs/araucana
cd "$D"
get() {  # -C - resumes a partial file rather than restarting 5.3 GB
  echo ">>> $1"
  curl -sS -L -C - --retry 20 --retry-delay 15 --connect-timeout 60 -o "$(basename "$1")" "$B/$1"
}
get v1/assembly/v1.fasta.gz
get v1/annotation/v1_annotation.gff3
get v1/assembly/v1_chloroplast.fasta
get v1/annotation/v1_annotation_proteins.fasta
echo "=== sizes ==="; ls -l
echo "=== gzip integrity ==="
gzip -t v1.fasta.gz && echo "v1.fasta.gz OK"
echo "=== contig/base count ==="
zcat v1.fasta.gz | awk '/^>/{n++;next}{b+=length($0)} END{printf "contigs: %d\nbases  : %d (%.2f Gb)\n",n,b,b/1e9}'
echo "=== chloroplast ==="
awk '/^>/{n++;next}{b+=length($0)} END{printf "seqs: %d  bases: %d\n",n,b}' v1_chloroplast.fasta
