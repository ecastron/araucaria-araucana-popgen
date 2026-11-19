#!/bin/bash
#SBATCH --job-name=callsnp
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/callsnp_%j.out
set -euo pipefail
module load SAMTOOLS/samtools.1.20
module load Bcftools/1.20
D=/home/ecastron/araucaria/refs/araucana
O=/home/ecastron/araucaria/results/refmap
cd "$O"
# restrict to positions with >=6x in >=3 of 8 libraries (the cr30-equivalent filter)
if [ ! -s callable.bed ]; then
  awk 'BEGIN{OFS="\t"} {n=0; for(i=3;i<=NF;i++) if($i>=6) n++; if(n>=3) print $1,$2-1,$2}' \
    depthmat_dup.tsv > callable.bed
fi
echo "callable positions: $(wc -l < callable.bed)"
bcftools mpileup -f "$D/v1.fasta" -b bamlist_dup.txt -R callable.bed \
    -q 30 -Q 20 -d 2000 -a AD,DP -Ou \
  | bcftools call -mv -Ov > snps_ref.vcf
echo "############ REFERENCE-BASED SNP CALLING, 8 LIBRARIES ############"
echo "variant sites  : $(grep -vc '^#' snps_ref.vcf)"
echo "biallelic SNPs : $(awk '!/^#/ && length($4)==1 && length($5)==1' snps_ref.vcf | wc -l)"
awk '!/^#/ && length($4)==1 && length($5)==1 {
   n=0; for(i=10;i<=NF;i++){split($i,g,":"); if(g[1]!="./." && g[1]!=".|.") n++}
   c[n]++ } END{ print "genotyped in N libraries:"; for(k=8;k>=1;k--) printf "  %d libs: %d SNPs\n",k,c[k] }' snps_ref.vcf
echo
echo "Benchmark: ipyrad de novo = 78,802 SNPs at cr30 across 177 libraries"
echo "(8 libraries here, so expect fewer; the ratio is what matters)"
