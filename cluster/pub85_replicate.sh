#!/bin/bash
#SBATCH -J pub85
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 8:00:00
#SBATCH -o /home/ecastron/araucaria/logs/pub85_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/pub85_%j.err
# Approximate the published 2021/2022 filter chain on OUR assemblies and see what survives.
# Published: ipyrad ct 85%, then VCFtools excluding het>50%, MAF<5%, missing>50%,
# individuals >50% missing -> 134 individuals, 49,122 biallelic SNPs.
set -uo pipefail
P=/home/ecastron/araucaria
R=$P/results/Q1/pub85
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH
module load Popgen/popgen-1.0 2>/dev/null
PY=/opt/miniconda/miniconda-py312_24.3.0-0/envs/ipyrad/bin/python
exec > "$R/pub85.txt" 2>&1
echo "===== approximating the published filter chain on our assemblies ====="
echo "generated $(date)"

# rank libraries by call rate (from the cr30 set, cheap and monotone with overall quality)
$PY - <<'PY'
import re
V="/home/ecastron/araucaria/rad/ipyrad_full/full_ct085_outfiles/full_ct085.cr30.vcf"
per=None; s=[]
for line in open(V):
    if line.startswith("##"): continue
    f=line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"): s=f[9:]; per=[0]*len(s); n=0; continue
    n+=1
    for k,c in enumerate(f[9:]):
        if c.split(":")[0] not in ("./.","."): per[k]+=1
rank=sorted(range(len(s)), key=lambda i:-per[i])
with open("/home/ecastron/araucaria/results/Q1/pub85/best134.txt","w") as o:
    for i in rank[:134]: o.write(s[i]+"\n")
print(f"library call rate over {n} cr30 sites: best {per[rank[0]]/n:.3f}, "
      f"134th {per[rank[133]]/n:.3f}, worst {per[rank[-1]]/n:.3f}")
PY

for CT in 085 091; do
  V=$P/rad/ipyrad_full/full_ct${CT}_outfiles/full_ct${CT}.vcf
  echo
  echo "########## ct 0.${CT#0} ##########"
  for MODE in all134 all177; do
    if [ "$MODE" = all134 ]; then KEEP="--keep $R/best134.txt"; else KEEP=""; fi
    out=$R/ct${CT}_${MODE}
    vcftools --vcf "$V" $KEEP --max-missing 0.5 --maf 0.05 --min-alleles 2 --max-alleles 2 \
             --recode --stdout 2>"$R/log_${CT}_${MODE}.txt" | gzip > "$out.vcf.gz"
    n=$(zcat "$out.vcf.gz" | grep -vc '^#')
    ind=$(zcat "$out.vcf.gz" | grep -m1 '^#CHROM' | awk '{print NF-9}')
    printf "  %-8s individuals=%3s  SNPs=%7s\n" "$MODE" "${ind:-0}" "$n"
  done
done
echo
echo "published target: 134 individuals, 49,122 SNPs"
