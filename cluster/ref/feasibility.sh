#!/bin/bash
#SBATCH --job-name=feas
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/feas_%j.out
set -euo pipefail
module load SAMTOOLS/samtools.1.20
O=/home/ecastron/araucaria/results/refmap
cd "$O"
# Per-position depth across all 8 libraries at once. samtools depth emits only
# positions with coverage somewhere, so this stays tractable on a 19.3 Gb genome.
# Same MAPQ/BQ filters as the balance analysis so the numbers are comparable.
for tag in nodup dup; do
  echo ">>> $tag"
  samtools depth -Q 30 -q 20 -f "bamlist_${tag}.txt" > "depthmat_${tag}.tsv"
  wc -l < "depthmat_${tag}.tsv" | xargs echo "   positions with any coverage:"
done
python3 - <<'PY'
def summarise(path,label):
    import sys
    tot=0
    # count libraries at >=6x and >=10x per position
    hist6=[0]*9; hist10=[0]*9
    for ln in open(path):
        f=ln.rstrip('\n').split('\t')
        if len(f)<3: continue
        d=[int(x) for x in f[2:]]
        tot+=1
        hist6[sum(1 for x in d if x>=6)]+=1
        hist10[sum(1 for x in d if x>=10)]+=1
    print('\n=== %s ==='%label)
    print('positions with any coverage: %d'%tot)
    print('%-28s %12s %12s'%('libraries meeting depth','>=6x','>=10x'))
    for k in range(8,-1,-1):
        print('  %-26d %12d %12d'%(k,hist6[k],hist10[k]))
    a6=sum(hist6[k] for k in range(6,9)); a10=sum(hist10[k] for k in range(6,9))
    a8_6=hist6[8]; a8_10=hist10[8]
    print('  in >=6 of 8 libraries      %12d %12d'%(a6,a10))
    print('  in ALL 8 libraries         %12d %12d'%(a8_6,a8_10))
summarise('/home/ecastron/araucaria/results/refmap/depthmat_nodup.tsv','DEDUPLICATED reads')
summarise('/home/ecastron/araucaria/results/refmap/depthmat_dup.tsv','ALL reads (with duplicates)')
print('\nBenchmark: the current analysis uses 78,802 SNPs across 177 libraries.')
print('A site is only useful if enough individuals retain depth there.')
PY
