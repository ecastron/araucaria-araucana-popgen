#!/bin/bash
#SBATCH -J dp20chk
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH -o /home/ecastron/araucaria/logs/dp20chk_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/dp20chk_%j.err
# Sensitivity check: re-run Q1/Q2/Q4 with genotypes restricted to >=20x, methods unchanged.
set -uo pipefail
P=/home/ecastron/araucaria
V=$P/results/dp20/full_ct091.dp20.cr30.vcf
D=$P/results/dp20
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH
LOG=$D/run_dp20.log; : > "$LOG"
exec 2>&1

echo "########## per-library calls WITHIN the final 4145-SNP matrix ##########"
python3 - "$V" <<'PY'
import sys, statistics as st
per=None; n=0
for line in open(sys.argv[1]):
    if line.startswith("##"): continue
    f=line.rstrip("\n").split("\t")
    if line.startswith("#CHROM"): per=[0]*(len(f)-9); continue
    n+=1
    for k,c in enumerate(f[9:]):
        if c.split(":")[0] not in ("./.","."): per[k]+=1
print(f"{n} SNPs x {len(per)} libraries")
print(f"per-library calls: median {int(st.median(per))} min {min(per)} max {max(per)}")
for t in (200,500,1000,2000):
    print(f"  libraries with >={t} calls: {sum(1 for x in per if x>=t)}/{len(per)}")
PY

echo; echo "########## Q1: popgen (distance, FST matrix, PCA, diversity) ##########"
python3 $P/scripts/dp20/q1_popgen.py "$V" 2>&1 | tail -30

echo; echo "########## Q1: AMOVA ##########"
module load R 2>/dev/null
Rscript $P/scripts/dp20/q1_amova.R 2>&1 | tail -25

echo; echo "########## Q1: FST via vcftools (input to Q4) ##########"
module load Popgen/popgen-1.0 2>/dev/null
mkdir -p "$D/Q1/fst"
while read A; do while read B; do
  [ "$A" \< "$B" ] || continue
  r=$(vcftools --vcf "$V" --weir-fst-pop "$D/Q1/fst/pop_$A.txt" \
        --weir-fst-pop "$D/Q1/fst/pop_$B.txt" --stdout 2>&1 \
        | grep -i "mean Fst" | awk '{print $NF}')
  printf "%s\t%s\t%s\n" "$A" "$B" "${r:-NA}"
done < "$D/Q1/fst/stands.txt"; done < "$D/Q1/fst/stands.txt" > "$D/Q1/q1_fst_vcftools.tsv"
echo "pairs written: $(wc -l < "$D/Q1/q1_fst_vcftools.tsv")"
head -3 "$D/Q1/q1_fst_vcftools.tsv"

echo; echo "########## Q4: IBD + Coastal residual ##########"
Rscript $P/scripts/dp20/q4_ibd_residual.R 2>&1 | tail -30

echo; echo "########## Q2: ngsRelate at 20x, then SGS ##########"
python3 $P/scripts/catg_to_beagle_mindp.py \
  "$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf" \
  "$D/Q2/gl.beagle" 0.01 20 2>&1 | tail -3
gzip -f "$D/Q2/gl.beagle"
NS=$(zcat "$D/Q2/gl.beagle.gz" | head -1 | awk '{print (NF-3)/3}')
NL=$(( $(zcat "$D/Q2/gl.beagle.gz" | wc -l) - 1 ))
zcat "$D/Q2/gl.beagle.gz" | head -1 | cut -f4- | tr '\t' '\n' | awk 'NR%3==1' > "$D/Q2/samples.txt"
ngsRelate -G "$D/Q2/gl.beagle.gz" -n "$NS" -L "$NL" -z "$D/Q2/samples.txt" \
  -O "$D/Q2/ngsrelate.tsv" > /dev/null 2>&1
echo "--- pairwise overlap (nSites) at 20x vs baseline ---"
python3 - "$D/Q2/ngsrelate.tsv" "$P/results/Q2/ngsrelate.tsv" <<'PY'
import sys, statistics as st
for lab,path in (("dp20",sys.argv[1]),("baseline",sys.argv[2])):
    try:
        hdr=None; ns=[]
        for line in open(path):
            p=line.rstrip("\n").split("\t")
            if hdr is None: hdr=p; continue
            ns.append(int(dict(zip(hdr,p))["nSites"]))
        print(f"  {lab:<9} pairs={len(ns)}  median nSites={int(st.median(ns))}  "
              f">=200: {sum(1 for x in ns if x>=200)}  >=2000: {sum(1 for x in ns if x>=2000)}")
    except Exception as e: print(f"  {lab}: {e}")
PY
Rscript $P/scripts/dp20/q2_sgs.R 2>&1 | tail -25
echo; echo "########## DONE $(date) ##########"
